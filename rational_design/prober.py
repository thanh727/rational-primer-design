import pandas as pd
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Blast import NCBIWWW, NCBIXML
import os
import time
import socket
import urllib


def _positive_product_mask(df):
    if df.empty or "PCR Product Sequence" not in df.columns:
        return pd.Series(False, index=df.index)
    series = df["PCR Product Sequence"]
    text = series.astype("string").str.strip().str.lower()
    return series.notna() & ~text.isin(["", "nan", "none", "n/a", "0", "absent"])


def _gc_percent(seq):
    seq = str(seq).upper()
    if not seq:
        return 0.0
    return (seq.count("G") + seq.count("C")) / len(seq) * 100


def _has_homopolymer(seq, limit=4):
    seq = str(seq).upper()
    run = 1
    last = ""
    for base in seq:
        if base == last:
            run += 1
            if run >= limit:
                return True
        else:
            run = 1
            last = base
    return False

def _gc_percent_degenerate(seq):
    seq_str = str(seq).upper()
    if not seq_str:
        return 0.0
    from .validator import IUPAC as v_iupac
    gc_fraction = 0.0
    for c in seq_str:
        if c in ('G', 'C'):
            gc_fraction += 1.0
        elif c in ('A', 'T'):
            gc_fraction += 0.0
        elif c in v_iupac:
            bases = v_iupac[c]
            gc_fraction += sum(1 for b in bases if b in ('G', 'C')) / len(bases)
        else:
            gc_fraction += 0.5
    return (gc_fraction / len(seq_str)) * 100

class ProbeSelector:
    def design(self, detail_csv, output_csv, summary_csv, params):
        print("--- 🌍 PROBER: Full Genomic Signature Atlas Generation ---")
        detail_dir = os.path.dirname(detail_csv)
        bg_candidates = [
            params.get("background_results_csv") if params else None,
            os.path.join(detail_dir, "master_results_background.csv"),
            os.path.join(detail_dir, "results_background.csv"),
        ]
        bg_csv = next((p for p in bg_candidates if p and os.path.exists(p)), None)
        df_sum = pd.read_csv(summary_csv)
        df_detail = pd.read_csv(detail_csv)
        df_bg = pd.read_csv(bg_csv) if bg_csv else pd.DataFrame()

        final_assays = []
        for _, row in df_sum.iterrows():
            pair_id = row['Primer Pair']
            relevant_amps = df_detail[df_detail['Primer Pair'] == pair_id]['PCR Product Sequence'].dropna().unique().tolist()
            if not relevant_amps: continue

            spec_val = 100.0
            if not df_bg.empty:
                bg_total = df_bg['Strain'].nunique()
                bg_hits = df_bg[(df_bg['Primer Pair'] == pair_id) & _positive_product_mask(df_bg)]['Strain'].nunique()
                spec_val = round((1 - (bg_hits / max(1, bg_total))) * 100, 2)

            p_tm_base = (self._calc_tm(row['Forward Primer']) + self._calc_tm(row['Reverse Primer'])) / 2
            best_p = self._find_best_probe(
                relevant_amps,
                p_tm_base,
                row['Forward Primer'],
                row['Reverse Primer'],
                params
            )

            final_assays.append({
                "Set_ID": pair_id,
                "Sensitivity": f"{row['Sensitivity_Percent']}%",
                "Spec": f"{spec_val}%",
                "Max_Copy": row['Max_Copy_Number'],
                "Fwd_Primer": row['Forward Primer'], "Rev_Primer": row['Reverse Primer'],
                "Amplicon_Size": len(relevant_amps[0]),
                "Amplicon_Sequence": relevant_amps[0],
                "Probe_Seq": best_p['seq'] if best_p else "N/A",
                "Probe_Tm": round(best_p['tm'], 1) if best_p else 0
            })

        if final_assays:
            pd.DataFrame(final_assays).to_csv(output_csv, index=False)
            return True
        return False

    def run_blast_annotation(self, final_csv_path):
        """Translate Amplicon and run BLASTX Protein for gene annotation."""
        if not os.path.exists(final_csv_path): return
        df = pd.read_csv(final_csv_path)
        print(f"🌐 Launching BLASTX Protein for {len(df)} candidates (Timeout: 30s)...")

        # Set strict global socket timeout for NCBI server to prevent hanging
        socket.setdefaulttimeout(30.0)

        df['Target_Gene'] = "N/A"

        for index, row in df.iterrows():
            amp_seq = row['Amplicon_Sequence']
            if pd.isna(amp_seq) or amp_seq == "N/A": continue

            print(f"   🔎 [{index+1}/{len(df)}] BLASTX: {row['Set_ID']}...")
            try:
                # blastx mode: translate nucleotide to protein and search
                result_handle = NCBIWWW.qblast("blastx", "nr", amp_seq)
                blast_record = NCBIXML.read(result_handle)

                if blast_record.alignments:
                    # Get first protein title
                    hit_title = blast_record.alignments[0].title
                    gene_name = hit_title.split('[')[0].split('>', 1)[-1].strip()
                    df.at[index, 'Target_Gene'] = gene_name
                    print(f"      🎯 Gene: {gene_name[:60]}...")

                # NCBI limits request frequency, sleep 5s to avoid being blocked
                time.sleep(5)

            except socket.timeout:
                print(f"   ⚠️ Timeout (30s) at {row['Set_ID']}. NCBI server is too slow, skipping annotation.")
                df.at[index, 'Target_Gene'] = "Timeout"
            except urllib.error.URLError:
                print(f"   ⚠️ Network Error/Timeout at {row['Set_ID']}. Skipping annotation.")
                df.at[index, 'Target_Gene'] = "Network Error"
            except Exception as e:
                print(f"   ⚠️ Error at {row['Set_ID']}: {e}")
                df.at[index, 'Target_Gene'] = "Error"

            # Save temporarily every 5 items to avoid data loss if network drops
            if index % 5 == 0: df.to_csv(final_csv_path, index=False)

        df.to_csv(final_csv_path, index=False)
        print(f"🎉 Complete annotation for all candidates!")

    def _calc_tm(self, seq):
        """IDT OligoAnalyzer / ThermoFisher-compatible Tm (SantaLucia 1998).
        Conditions: Na=50mM, Mg=3.0mM, dNTPs=0.8mM, [Oligo]=250nM."""
        clean = "".join([c for c in str(seq).upper() if c in 'ATGC'])
        return mt.Tm_NN(
            Seq(clean),
            nn_table=mt.DNA_NN3,
            Na=50,
            Mg=3.0,
            dNTPs=0.8,
            dnac1=250,
            dnac2=0
        )

    def _find_best_probe(self, amplicons, p_tm, fwd_primer="", rev_primer="", params=None):
        ref_amp = min(amplicons, key=len)
        left = max(15, len(str(fwd_primer)))
        right = max(15, len(str(rev_primer)))
        if len(ref_amp) <= left + right + 20:
            return None
        search_zone = ref_amp[left:len(ref_amp)-right]
        candidates = []
        from .multiplex import check_dimer, check_hairpin

        degenerate_probes = params.get("degenerate_primers", False) if params else False
        max_iupac = int(params.get("max_iupac_per_primer", 2)) if params else 2

        for length in range(20, 26):
            for i in range(len(search_zone) - length + 1):
                sub = search_zone[i : i+length]
                probe_seq = None

                if any(sub not in amp for amp in amplicons):
                    if degenerate_probes:
                        # Extract matching segments from amplicons and build consensus
                        candidate_seqs = []
                        fwd_len = len(fwd_primer)
                        rev_len = len(rev_primer)
                        import Levenshtein
                        for amp in amplicons:
                            amp_inner = amp[fwd_len : len(amp) - rev_len]
                            ref_inner = ref_amp[fwd_len : len(ref_amp) - rev_len]
                            if len(amp_inner) == len(ref_inner):
                                rel_pos = left + i - fwd_len
                                candidate_seqs.append(amp_inner[rel_pos : rel_pos + length])
                            else:
                                best_sub = None
                                best_dist = 999
                                for k in range(len(amp_inner) - length + 1):
                                    candidate = amp_inner[k : k + length]
                                    dist = Levenshtein.distance(sub, candidate)
                                    if dist < best_dist:
                                        best_dist = dist
                                        best_sub = candidate
                                if best_sub:
                                    candidate_seqs.append(best_sub)

                        if len(candidate_seqs) == len(amplicons):
                            from .utils import generate_iupac_consensus
                            probe_seq = generate_iupac_consensus(candidate_seqs, max_iupac=max_iupac)
                        else:
                            continue
                    else:
                        continue
                else:
                    probe_seq = sub

                if not probe_seq:
                    continue

                tm = self._calc_tm(probe_seq)
                gc = _gc_percent_degenerate(probe_seq)
                if not (p_tm + 5 <= tm <= p_tm + 12):
                    continue
                if not (35 <= gc <= 65):
                    continue
                if probe_seq.startswith("G") or _has_homopolymer(probe_seq):
                    continue
                if check_hairpin(probe_seq) >= 4:
                    continue
                if fwd_primer and check_dimer(probe_seq, str(fwd_primer))[1] >= 4:
                    continue
                if rev_primer and check_dimer(probe_seq, str(rev_primer))[1] >= 4:
                    continue
                candidates.append({'seq': probe_seq, 'tm': tm, 'gc': gc})
        return sorted(candidates, key=lambda x: (abs(x['tm'] - (p_tm + 8)), abs(x['gc'] - 50)))[0] if candidates else None
