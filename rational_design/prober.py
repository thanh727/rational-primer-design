import pandas as pd
import numpy as np
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqUtils import gc_fraction
from Bio.Blast import NCBIWWW, NCBIXML
import primer3
import warnings
import time
import os
import sys

try:
    import Levenshtein
    HAS_LEVENSHTEIN = True
except ImportError:
    HAS_LEVENSHTEIN = False

warnings.filterwarnings("ignore", category=UserWarning, module='primer3')

MIN_GC, MAX_GC = 20.0, 80.0
MIN_LEN, MAX_LEN = 15, 36
AVOID_5_PRIME_G = True  

class ProbeSelector:
    def __init__(self):
        self.enable_blast = True

    def design(self, validation_csv, output_csv, stats_csv=None, params=None):
        print("--- STARTING PROBE DESIGNER (Instrumented) ---")
        t_start = time.time()
        
        if params is None: params = {}
        min_sens = params.get("min_sensitivity", 0.0)
        self.probe_max_mm = params.get("probe_max_mismatch", 3)

        try: df = pd.read_csv(validation_csv)
        except: return False
        if df.empty: return False

        stats_map = {}
        if stats_csv:
            try:
                df_stats = pd.read_csv(stats_csv)
                for _, row in df_stats.iterrows():
                    pair = row["Primer Pair"]
                    stats_map[pair] = {"sens": row.get("Sensitivity_Percent", 0), "spec": 100.0 - row.get("Non-Target_Percent", 0)}
            except: pass

        df = df[df['PCR Product Sequence'].notna() & (df['PCR Product Sequence'] != "")]
        grouped = df.groupby('Primer Pair')
        final_rows = []
        
        print(f"   [Filter] Applying Min Sensitivity Threshold: {min_sens}%")
        
        for set_id, group in grouped:
            stats = stats_map.get(set_id, {"sens": 0, "spec": 0})
            if stats["sens"] < min_sens: continue

            amplicons = group['PCR Product Sequence'].tolist()
            fwd = str(group.iloc[0]['Forward Primer'])
            rev_input = str(group.iloc[0]['Reverse Primer'])
            
            tm_f = self._calc_tm_robust(fwd)
            tm_r = self._calc_tm_robust(rev_input)
            avg_tm = (tm_f + tm_r) / 2
            
            fwd_stats = self.get_structure_stats(fwd)
            rev_stats = self.get_structure_stats(rev_input)
            ht_tm, ht_dg = self.get_heterodimer_stats(fwd, rev_input)
            
            top3 = self.find_probes_priority(amplicons, fwd, rev_input, set_id, avg_tm)
            
            if top3:
                print(f"   [+] {set_id}: Found {len(top3)} probes.")
                for i, p in enumerate(top3):
                    rev_order = str(Seq(rev_input).reverse_complement())
                    final_rows.append({
                        "Set_ID": set_id, "Option": i+1, "Target_Gene": "Pending_BLAST",
                        "Sensitivity": f"{stats['sens']:.1f}%", "Specificity": f"{stats['spec']:.1f}%",
                        "Probe_Coverage": f"{p['coverage_pct']:.1f}%", "Amplicon_Size": len(amplicons[0]),
                        "Note": p['type'],
                        "Fwd_Seq": fwd, "Fwd_Tm": round(tm_f, 1), "Fwd_GC": round(gc_fraction(fwd)*100,1),
                        "Fwd_Hairpin_dG": fwd_stats['dG'], "Fwd_HomoDimer_dG": fwd_stats['Homo_dG'],
                        "Rev_Seq_Order": rev_order, "Rev_Tm": round(tm_r, 1), "Rev_GC": round(gc_fraction(rev_input)*100,1),
                        "Rev_Hairpin_dG": rev_stats['dG'], "Rev_HomoDimer_dG": rev_stats['Homo_dG'],
                        "HeteroDimer_dG": ht_dg,
                        "Probe_Seq": p['seq'], "Probe_Tm": round(p['tm'], 1), "Probe_GC": round(p['gc'], 1),
                        "Tm_Gap": round(p['tm_gap'], 1), "Full_Amplicon": amplicons[0]
                    })

        if final_rows:
            cols = ["Set_ID", "Option", "Target_Gene", "Sensitivity", "Specificity", "Probe_Coverage", "Amplicon_Size", "Note", "Fwd_Seq", "Fwd_Tm", "Rev_Seq_Order", "Rev_Tm", "Probe_Seq", "Probe_Tm", "Probe_GC", "Tm_Gap", "Full_Amplicon"]
            df_out = pd.DataFrame(final_rows)
            df_out = df_out[[c for c in cols if c in df_out.columns] + [c for c in df_out.columns if c not in cols]]
            self._safe_save(df_out, output_csv)
            self.last_saved_path = output_csv
            return True
        else:
            print("\n❌ No valid probes found.")
            return False

    def run_blast_annotation(self, csv_path):
        if hasattr(self, 'last_saved_path'): csv_path = self.last_saved_path
        print("\n--- STARTING SMART BLAST ANNOTATION (Optimized) ---")
        try: df = pd.read_csv(csv_path)
        except: return
        
        unique_sets = df['Set_ID'].unique()
        print(f"🔎 BLASTing {len(unique_sets)} unique primer sets...")
        
        for set_id in unique_sets:
            mask = df['Set_ID'] == set_id
            first_val = str(df.loc[mask, 'Target_Gene'].iloc[0])
            if first_val != "Pending_BLAST" and "Hypothetical" not in first_val:
                continue

            amplicon = df.loc[mask, 'Full_Amplicon'].iloc[0]
            print(f"   [Querying] {set_id}...")
            gene_name = self.identify_target_gene(amplicon)
            df.loc[mask, 'Target_Gene'] = gene_name
            print(f"      -> {gene_name}")
            self._safe_save(df, csv_path)
        print(f"✅ Annotation Complete!")

    def _safe_save(self, df, path):
        counter = 1; save_path = path
        while True:
            try: df.to_csv(save_path, index=False); break
            except PermissionError:
                counter += 1
                save_path = f"{os.path.splitext(path)[0]}_{counter}{os.path.splitext(path)[1]}"

    def _calc_tm_robust(self, seq_str):
        try: return mt.Tm_NN(Seq(seq_str), nn_table=mt.DNA_NN3, Na=50, Mg=1.5, dNTPs=0.2, dnac1=500, dnac2=500)
        except: return 0.0

    def get_structure_stats(self, seq_str):
        try: return {'dG': round(primer3.calcHairpin(seq_str).dg, 1), 'Homo_dG': round(primer3.calcHomodimer(seq_str).dg, 1)}
        except: return {'dG': 0.0, 'Homo_dG': 0.0}

    def get_heterodimer_stats(self, fwd, rev):
        try:
            ht = primer3.calcHeterodimer(fwd, rev)
            return round(ht.tm, 1), round(ht.dg, 1)
        except: return 0.0, 0.0

    def identify_target_gene(self, seq_str):
        if not self.enable_blast: return "BLAST_Disabled"
        try:
            result = NCBIWWW.qblast("blastx", "nr", seq_str, hitlist_size=1)
            record = NCBIXML.read(result)
            if record.alignments:
                title = record.alignments[0].title
                if "|" in title: return title.split("|")[-1].strip()
                return title[:50]
        except: pass
        return "Hypothetical/Unknown"

    def find_probes_priority(self, amplicons, fwd_primer, rev_primer_input, set_id, avg_primer_tm):
        ref_amp = min(amplicons, key=len)
        start_buf, end_buf = 3, 3
        if len(ref_amp) < (start_buf + end_buf + 15): return []
        search_space = ref_amp[start_buf : -end_buf]
        candidates = []
        seen = set()
        
        for length in range(MIN_LEN, MAX_LEN + 1):
            for i in range(len(search_space) - length + 1):
                seq = search_space[i : i+length]
                if seq in seen: continue
                
                valid, reason = self.check_composition(seq)
                if not valid: continue
                
                tm = self._calc_tm_robust(seq)
                tm_gap = tm - avg_primer_tm
                
                probe_type = "Unknown"; base_score = 0
                if tm_gap >= 5.0:
                    probe_type = "Standard TaqMan"; base_score = 1000 
                elif tm_gap >= -5.0:
                    probe_type = "Relaxed/MGB"; base_score = 100 
                else: continue 

                n_ex, n_mis = self.check_coverage(seq, amplicons)
                cov_pct = ((n_ex + n_mis) / len(amplicons)) * 100.0
                
                len_bonus = 0
                if probe_type == "Standard TaqMan":
                    if 20 <= length <= 30: len_bonus = 50
                else:
                    if 13 <= length <= 20: len_bonus = 50
                
                final_score = base_score + (cov_pct * 5) + len_bonus + tm_gap
                
                candidates.append({
                    "seq": seq, "tm": tm, "gc": gc_fraction(seq)*100,
                    "type": probe_type, "tm_gap": tm_gap,
                    "coverage_pct": cov_pct, "score": final_score
                })
                seen.add(seq)
                
        candidates.sort(key=lambda x: x["score"], reverse=True)
        return candidates[:3]

    def fast_hamming(self, s1, s2):
        if len(s1) != len(s2): return 99
        if HAS_LEVENSHTEIN: return Levenshtein.hamming(s1, s2)
        dist = 0
        for i in range(len(s1)):
            if s1[i] != s2[i]: dist += 1
        return dist

    def check_coverage(self, probe_seq, amplicons):
        exact = 0; mismatch = 0; p_len = len(probe_seq)
        for amp in amplicons:
            if probe_seq in amp: exact += 1; continue
            found = False
            scan_limit = len(amp) - p_len + 1
            if scan_limit > 0:
                for k in range(scan_limit):
                    sub_seq = amp[k : k+p_len]
                    if self.fast_hamming(probe_seq, sub_seq) <= self.probe_max_mm:
                        found = True; break
            if found: mismatch += 1
        return exact, mismatch

    def check_composition(self, seq):
        if "N" in seq: return False, "Has N"
        if AVOID_5_PRIME_G and seq.startswith("G"): return False, "5' G"
        if "GGGG" in seq or "CCCC" in seq: return False, "Poly Run"
        gc = gc_fraction(seq) * 100
        if not (MIN_GC <= gc <= MAX_GC): return False, "GC"
        return True, "OK"