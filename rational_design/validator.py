import csv
import multiprocessing
import os
from collections import defaultdict

import Levenshtein
import pandas as pd
import psutil
from Bio import SeqIO
from Bio.Seq import Seq

from .utils import nuclear_ram_flush

DEFAULT_PRODUCT_MIN, DEFAULT_PRODUCT_MAX = 70, 450

IUPAC = {
    "A": {"A"}, "C": {"C"}, "G": {"G"}, "T": {"T"},
    "R": {"A", "G"}, "Y": {"C", "T"}, "S": {"G", "C"},
    "W": {"A", "T"}, "K": {"G", "T"}, "M": {"A", "C"},
    "B": {"C", "G", "T"}, "D": {"A", "G", "T"},
    "H": {"A", "C", "T"}, "V": {"A", "C", "G"},
    "N": {"A", "C", "G", "T"},
}


def _is_degenerate(seq):
    return any(base not in "ATGC" for base in str(seq).upper())


def _iupac_mismatch_count(primer_s, target_s):
    return sum(1 for p, t in zip(primer_s, target_s) if t not in IUPAC.get(p, {p}))


def _primer_segments(primer_s, max_mm):
    """Pigeonhole segments: any <=max_mm alignment must contain one exact segment."""
    n = len(primer_s)
    segment_count = max(1, min(max_mm + 1, n))
    base_len, extra = divmod(n, segment_count)
    offset = 0
    for idx in range(segment_count):
        seg_len = base_len + (1 if idx < extra else 0)
        if seg_len > 0:
            yield offset, primer_s[offset:offset + seg_len]
        offset += seg_len


def positive_product_mask(df):
    if df.empty or "PCR Product Sequence" not in df.columns:
        return pd.Series(False, index=df.index)
    series = df["PCR Product Sequence"]
    text = series.astype("string").str.strip().str.lower()
    return series.notna() & ~text.isin(["", "nan", "none", "n/a", "0", "absent"])


def find_primer_hits(sequence_s, primer_s, max_mm):
    hits = []
    primer_s = str(primer_s).upper()
    sequence_s = str(sequence_s).upper()
    n = len(primer_s)
    if n == 0 or len(sequence_s) < n:
        return hits

    degenerate = _is_degenerate(primer_s)
    candidate_starts = set()

    if degenerate:
        candidate_starts.update(range(0, len(sequence_s) - n + 1))
    else:
        for offset, segment in _primer_segments(primer_s, max_mm):
            start = 0
            while True:
                pos = sequence_s.find(segment, start)
                if pos == -1:
                    break
                start_pos = pos - offset
                if 0 <= start_pos <= len(sequence_s) - n:
                    candidate_starts.add(start_pos)
                start = pos + 1

    for start_pos in sorted(candidate_starts):
        target_sub = sequence_s[start_pos:start_pos + n]
        if degenerate:
            mm_count = _iupac_mismatch_count(primer_s, target_sub)
        else:
            mm_count = Levenshtein.hamming(primer_s, target_sub)
        if mm_count <= max_mm:
            mapping = "".join([t if p != t else "-" for p, t in zip(primer_s, target_sub)])
            hits.append({"pos": start_pos, "map": mapping, "mm": mm_count, "seq": target_sub})
    return hits


ALL_PRIMERS_CACHED = []


def init_worker(primer_list):
    global ALL_PRIMERS_CACHED
    ALL_PRIMERS_CACHED = primer_list


def run_pcr_on_single_strain(strain_task):
    global ALL_PRIMERS_CACHED
    sid, sequences = strain_task
    results = []

    for p_name, fwd, rev_rc, rev_orig, mm, product_min, product_max in ALL_PRIMERS_CACHED:
        f_s, r_rc_s = fwd, rev_rc
        r_s, f_rc_s = rev_orig, str(Seq(fwd).reverse_complement()).upper()

        best_hit, min_mm, total_bands_in_strain = None, 999, 0

        for seq_s in sequences:
            f_hits = find_primer_hits(seq_s, f_s, mm)
            r_hits = find_primer_hits(seq_s, r_rc_s, mm)

            r_len = len(r_hits)
            r_idx = 0
            for f in f_hits:
                while r_idx < r_len and r_hits[r_idx]["pos"] <= f["pos"]:
                    r_idx += 1
                for j in range(r_idx, r_len):
                    r = r_hits[j]
                    p_len = r["pos"] - f["pos"] + len(r_rc_s)
                    if p_len > product_max:
                        break
                    if product_min <= p_len:
                        total_bands_in_strain += 1
                        if (f["mm"] + r["mm"]) < min_mm:
                            min_mm = f["mm"] + r["mm"]
                            best_hit = {
                                "seq": seq_s[f["pos"]:r["pos"] + len(r_rc_s)],
                                "map": f"F:[{f['map']}]; R:[{r['map'][::-1]}]",
                                "f_bind": f["seq"],
                                "r_bind": str(Seq(r["seq"]).reverse_complement()).upper(),
                            }

            r_hits_alt = find_primer_hits(seq_s, r_s, mm)
            f_hits_alt = find_primer_hits(seq_s, f_rc_s, mm)

            f_len = len(f_hits_alt)
            f_idx = 0
            for r in r_hits_alt:
                while f_idx < f_len and f_hits_alt[f_idx]["pos"] <= r["pos"]:
                    f_idx += 1
                for j in range(f_idx, f_len):
                    f = f_hits_alt[j]
                    p_len = f["pos"] - r["pos"] + len(f_rc_s)
                    if p_len > product_max:
                        break
                    if product_min <= p_len:
                        total_bands_in_strain += 1
                        if (r["mm"] + f["mm"]) < min_mm:
                            min_mm = r["mm"] + f["mm"]
                            best_hit = {
                                "seq": seq_s[r["pos"]:f["pos"] + len(f_rc_s)],
                                "map": f"F:[{f['map'][::-1]}]; R:[{r['map']}]",
                                "f_bind": str(Seq(f["seq"]).reverse_complement()).upper(),
                                "r_bind": r["seq"],
                            }

        results.append([
            sid,
            p_name,
            best_hit["seq"] if best_hit else None,
            total_bands_in_strain,
            best_hit["map"] if best_hit else "N/A",
            fwd,
            rev_orig,
            best_hit["f_bind"] if best_hit else None,
            best_hit["r_bind"] if best_hit else None,
        ])
    return results


class InSilicoValidator:
    def __init__(self, params=None):
        self.params = params or {}

    def validate(self, target_path, bg_path, cand_fa, out_dir, params):
        self.params = params
        cpu = max(1, int(self.params.get("cpu_cores", 4)))
        max_mm = int(self.params.get("max_mismatch", 2))
        all_pairs = self._prepare_pairs(cand_fa, max_mm)

        for label, path in [("Background", bg_path), ("Target", target_path)]:
            if not path or not os.path.exists(path):
                continue
            out_csv = os.path.join(out_dir, f"results_{label.lower()}.csv")
            self._process_strictly(path, all_pairs, out_csv, cpu, label)
            if label == "Background":
                valid_ids = self._identify_specific(out_csv)
                all_pairs = [p for p in all_pairs if p[0] in valid_ids]

        self._generate_summary(os.path.join(out_dir, "results_target.csv"), out_dir)
        if self.params.get("degenerate_primers", False) and bg_path and os.path.exists(bg_path):
            self._revalidate_degenerate_background(bg_path, out_dir, cpu, max_mm)

    def _process_strictly(self, input_path, pairs, out_csv, cpu, label):
        strain_data = defaultdict(list)
        if os.path.isfile(input_path):
            for r in SeqIO.parse(input_path, "fasta"):
                sid = r.id.split("|")[0] if "|" in r.id else r.id
                strain_data[sid].append(str(r.seq).upper())
        else:
            for f in sorted(os.listdir(input_path)):
                if f.lower().endswith((".fasta", ".fa", ".fna", ".fas")):
                    sid = f
                    for r in SeqIO.parse(os.path.join(input_path, f), "fasta"):
                        strain_data[sid].append(str(r.seq).upper())

        all_sids = list(strain_data.keys())
        print(f"      > Detected {len(all_sids)} independent strains in {label}.")

        with open(out_csv, "w", newline="") as f:
            csv.writer(f).writerow([
                "Strain", "Primer Pair", "PCR Product Sequence", "Total_Bands",
                "Mismatch_Map", "Forward Primer", "Reverse Primer", "F_Bind_Seq", "R_Bind_Seq",
            ])

        if not pairs or not all_sids:
            return

        avail_gb = psutil.virtual_memory().available / (1024 ** 3)
        if avail_gb > 16:
            batch_size = 100
        elif avail_gb > 8:
            batch_size = 50
        elif avail_gb > 4:
            batch_size = 20
        else:
            batch_size = 5
        print(f"      > Free RAM: {avail_gb:.1f}GB. Auto-assigned BATCH_SIZE = {batch_size}")

        for i in range(0, len(all_sids), batch_size):
            batch_sids = all_sids[i:i + batch_size]
            batch_tasks = [(sid, strain_data[sid]) for sid in batch_sids]
            ctx = multiprocessing.get_context("spawn")
            with ctx.Pool(cpu, initializer=init_worker, initargs=(pairs,), maxtasksperchild=1) as pool:
                results = pool.map(run_pcr_on_single_strain, batch_tasks)
            with open(out_csv, "a", newline="") as f:
                writer = csv.writer(f)
                for res_list in results:
                    writer.writerows(res_list)
            del batch_tasks
            nuclear_ram_flush()

    def _prepare_pairs(self, p_fa, mm):
        recs = list(SeqIO.parse(p_fa, "fasta"))
        if len(recs) % 2 != 0:
            raise ValueError(f"Primer FASTA must contain forward/reverse record pairs: {p_fa}")

        product_min = int(self.params.get("product_size_min", DEFAULT_PRODUCT_MIN))
        product_max = int(self.params.get("product_size_max", DEFAULT_PRODUCT_MAX))
        pairs = []
        for i in range(0, len(recs), 2):
            pairs.append((
                recs[i].id.replace("_Fwd", ""),
                str(recs[i].seq).upper(),
                str(Seq(recs[i + 1].seq).reverse_complement()).upper(),
                str(recs[i + 1].seq).upper(),
                mm,
                product_min,
                product_max,
            ))
        return pairs

    def _identify_specific(self, csv_nt):
        df = pd.read_csv(csv_nt)
        if df.empty:
            return set()
        hits = df[positive_product_mask(df)].groupby("Primer Pair")["Strain"].nunique()
        total = df["Strain"].nunique()
        if total == 0:
            return set(df["Primer Pair"].unique())
        max_xr = float(self.params.get("validation_max_cross_reactivity", 1.0))
        return {p for p in df["Primer Pair"].unique() if (hits.get(p, 0) / total * 100) <= max_xr}

    def _generate_summary(self, t_csv, out_dir):
        df_raw = pd.read_csv(t_csv)
        df_pos = df_raw[positive_product_mask(df_raw)]
        total_strains = df_raw["Strain"].nunique()
        summary = []

        is_degenerate = self.params.get("degenerate_primers", False)
        if is_degenerate:
            from .utils import generate_iupac_consensus

        for p in df_raw["Primer Pair"].unique():
            p_df = df_pos[df_pos["Primer Pair"] == p]
            if p_df.empty:
                continue

            fwd_seq = p_df.iloc[0]["Forward Primer"]
            rev_seq = p_df.iloc[0]["Reverse Primer"]

            if is_degenerate and "F_Bind_Seq" in p_df.columns and "R_Bind_Seq" in p_df.columns:
                f_binds = p_df["F_Bind_Seq"].dropna().tolist()
                r_binds = p_df["R_Bind_Seq"].dropna().tolist()
                if f_binds and r_binds:
                    fwd_seq = generate_iupac_consensus(f_binds)
                    rev_seq = generate_iupac_consensus(r_binds)

            summary.append({
                "Primer Pair": p,
                "Sensitivity_Percent": round((p_df["Strain"].nunique() / max(1, total_strains)) * 100, 2),
                "Max_Copy_Number": int(p_df["Total_Bands"].max()),
                "PCR Product Sequence": p_df.iloc[0]["PCR Product Sequence"],
                "Forward Primer": fwd_seq,
                "Reverse Primer": rev_seq,
            })

        output_csv = os.path.join(out_dir, "pcr_results_summary.csv")
        columns = [
            "Primer Pair", "Sensitivity_Percent", "Max_Copy_Number",
            "PCR Product Sequence", "Forward Primer", "Reverse Primer",
        ]
        if summary:
            pd.DataFrame(summary).sort_values("Max_Copy_Number", ascending=False).to_csv(output_csv, index=False)
        else:
            pd.DataFrame(columns=columns).to_csv(output_csv, index=False)

    def _revalidate_degenerate_background(self, bg_path, out_dir, cpu, max_mm):
        summary_csv = os.path.join(out_dir, "pcr_results_summary.csv")
        if not os.path.exists(summary_csv):
            return
        df_summary = pd.read_csv(summary_csv)
        if df_summary.empty:
            return

        tmp_fasta = os.path.join(out_dir, "degenerate_candidates.fasta")
        with open(tmp_fasta, "w") as handle:
            for _, row in df_summary.iterrows():
                pair_id = row["Primer Pair"]
                handle.write(f">{pair_id}_Fwd\n{row['Forward Primer']}\n")
                handle.write(f">{pair_id}_Rev\n{row['Reverse Primer']}\n")

        degenerate_pairs = self._prepare_pairs(tmp_fasta, max_mm)
        bg_csv = os.path.join(out_dir, "results_background.csv")
        self._process_strictly(bg_path, degenerate_pairs, bg_csv, cpu, "Background")
        valid_ids = self._identify_specific(bg_csv)
        df_summary = df_summary[df_summary["Primer Pair"].isin(valid_ids)]
        df_summary.to_csv(summary_csv, index=False)
