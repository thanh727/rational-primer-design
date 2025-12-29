import os
import csv
import multiprocessing
from Bio import SeqIO
from Bio.Seq import Seq
from tqdm import tqdm
import pandas as pd
import time
try:
    import Levenshtein
    HAS_LEVENSHTEIN = True
except ImportError:
    HAS_LEVENSHTEIN = False

PRODUCT_LENGTH_LIMIT = 4000

def fast_hamming(s1, s2, max_mm):
    if len(s1) != len(s2): return -1
    if HAS_LEVENSHTEIN:
        d = Levenshtein.hamming(s1, s2)
        return d if d <= max_mm else -1
    mm = 0
    for i in range(len(s1)):
        if s1[i] != s2[i]:
            mm += 1
            if mm > max_mm: return -1
    return mm

def find_fuzzy_locations_pigeonhole(sequence, primer, max_mm):
    candidates = set()
    n_parts = max_mm + 1
    part_len = len(primer) // n_parts
    if part_len < 4:
        locs = []
        plen = len(primer)
        for i in range(len(sequence) - plen + 1):
            sub = sequence[i : i+plen]
            if fast_hamming(sub, primer, max_mm) != -1: locs.append(i)
        return locs
    for i in range(n_parts):
        start = i * part_len
        end = (i + 1) * part_len if i < n_parts - 1 else len(primer)
        chunk = primer[start:end]
        search_pos = 0
        while True:
            found_at = sequence.find(chunk, search_pos)
            if found_at == -1: break
            full_start = found_at - start
            if 0 <= full_start <= len(sequence) - len(primer): candidates.add(full_start)
            search_pos = found_at + 1
    valid_locs = []
    for loc in candidates:
        sub = sequence[loc : loc + len(primer)]
        if fast_hamming(sub, primer, max_mm) != -1: valid_locs.append(loc)
    return sorted(valid_locs)

def check_strand(sequence, fwd_seq, rev_seq, max_mm):
    fwd_locs = find_fuzzy_locations_pigeonhole(sequence, fwd_seq, max_mm)
    if not fwd_locs: return None
    rev_locs = find_fuzzy_locations_pigeonhole(sequence, rev_seq, max_mm)
    if not rev_locs: return None
    for f_start in fwd_locs:
        for r_start in rev_locs:
            if f_start < r_start:
                dist = r_start - f_start
                full_len = dist + len(rev_seq)
                if full_len < PRODUCT_LENGTH_LIMIT:
                    product_seq = sequence[f_start : f_start + full_len]
                    f_sub = sequence[f_start : f_start + len(fwd_seq)]
                    r_sub = sequence[r_start : r_start + len(rev_seq)]
                    f_mm = fast_hamming(f_sub, fwd_seq, 100)
                    r_mm = fast_hamming(r_sub, rev_seq, 100)
                    return product_seq, full_len, f_mm, r_mm
    return None

def process_genome_task(args):
    genome_id, genome_seq, primer_pairs, max_mm = args
    results = {}
    rc_genome = str(Seq(genome_seq).reverse_complement())
    for p_name, _, fwd, _, rev_search, rev_orig in primer_pairs:
        res = check_strand(genome_seq, fwd, rev_search, max_mm)
        if not res: res = check_strand(rc_genome, fwd, rev_search, max_mm)
        prod_seq, prod_len, f_mm, r_mm = res if res else (None, None, -1, -1)
        results[p_name] = {
            "Sequence": prod_seq if prod_seq else "",
            "Length": prod_len if prod_len else "No PCR product",
            "Fwd Primer Seq": fwd,
            "Rev Primer Seq": rev_orig,
            "Fwd Mismatches": f_mm if f_mm != -1 else 'N/A',
            "Rev Mismatches": r_mm if r_mm != -1 else 'N/A'
        }
    return (genome_id, results)

class InSilicoValidator:
    def __init__(self):
        pass

    def validate(self, target_fasta, bg_fasta, candidate_fasta, output_dir, params):
        print(f"--- STARTING VALIDATOR (Turbo Pigeonhole Mode) ---")
        out_path = output_dir
        if not os.path.exists(out_path): os.makedirs(out_path)
        
        max_mm = params.get("max_mismatch", 3)
        cpu = params.get("cpu_cores", 0)
        if cpu <= 0: cpu = max(1, multiprocessing.cpu_count() - 2)

        # 1. Background
        t0 = time.time()
        print("\n--- [1/4] Checking Non-Target (Background) ---")
        csv_nontarget = os.path.join(output_dir, "results_non_target.csv")
        self._run_pcr_sim(bg_fasta, candidate_fasta, csv_nontarget, max_mm, cpu)
        print(f"   [Time] {time.time()-t0:.1f}s")
        
        # 2. Filter
        t0 = time.time()
        print("\n--- [2/4] Filtering Specific Primers ---")
        max_xr = params.get("validation_max_cross_reactivity", 5.0)
        valid_pairs = self._identify_specific_primers(csv_nontarget, max_xr)
        print(f"   [Time] {time.time()-t0:.1f}s")
        
        if not valid_pairs:
            print("❌ Pipeline Stopped: No primers passed specificity check.")
            return

        # 3. Target
        t0 = time.time()
        print("\n--- [3/4] Checking Target Sensitivity ---")
        csv_target = os.path.join(output_dir, "results_target.csv")
        self._run_pcr_sim(target_fasta, candidate_fasta, csv_target, max_mm, cpu)
        print(f"   [Time] {time.time()-t0:.1f}s")
        
        # 4. Stats
        print("\n--- [4/4] Generating Summary Report ---")
        self._analyze_results(csv_target, csv_nontarget, valid_pairs, output_dir)

    def _run_pcr_sim(self, genome_file, primer_file, output_csv, max_mm, num_workers):
        print(f"🧬 Loading Genomes: {os.path.basename(genome_file)}")
        try: genome_dict = {r.id: str(r.seq).upper() for r in SeqIO.parse(genome_file, "fasta")}
        except: return
        if not genome_dict: return

        primers_list = []
        try:
            for r in SeqIO.parse(primer_file, "fasta"): primers_list.append((r.id, str(r.seq).upper()))
        except: return
        pairs = self._convert_primers_to_pairs(primers_list)
        if not pairs: return
        
        print(f"🖥️ Simulating PCR on {len(genome_dict)} genomes using {num_workers} cores...")
        tasks = [(gid, seq, pairs, max_mm) for gid, seq in genome_dict.items()]
        
        all_rows = []
        with multiprocessing.Pool(num_workers) as pool:
            results = list(tqdm(pool.imap_unordered(process_genome_task, tasks, chunksize=1), total=len(tasks)))
            
        for gid, res_dict in results:
            for p_name, data in res_dict.items():
                all_rows.append({
                    "Genome": gid, "Primer Pair": p_name,
                    "PCR Product Length": data['Length'], "PCR Product Sequence": data['Sequence'],
                    "Forward Primer": data['Fwd Primer Seq'], "Reverse Primer": data['Rev Primer Seq'],
                    "Fwd Mismatches": data['Fwd Mismatches'], "Rev Mismatches": data['Rev Mismatches']
                })
        
        all_rows.sort(key=lambda x: (x['Primer Pair'], x['Genome']))

        keys = ["Genome", "Primer Pair", "PCR Product Length", "PCR Product Sequence", 
                "Forward Primer", "Reverse Primer", "Fwd Mismatches", "Rev Mismatches"]
        with open(output_csv, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=keys)
            writer.writeheader()
            writer.writerows(all_rows)

    def _convert_primers_to_pairs(self, primers_list):
        pairs = []
        for i in range(0, len(primers_list) - 1, 2):
            fwd_id, fwd_seq = primers_list[i]
            rev_id, rev_seq_original = primers_list[i+1]
            rev_seq_search = str(Seq(rev_seq_original).reverse_complement())
            pair_name = fwd_id.split("_Fwd")[0]
            if "Set" not in pair_name: pair_name = f"Set_{i//2 + 1}"
            pairs.append((pair_name, fwd_id, fwd_seq, rev_id, rev_seq_search, rev_seq_original))
        return pairs

    def _identify_specific_primers(self, csv_file, max_xr):
        df = pd.read_csv(csv_file)
        if df.empty: return set()
        total = df['Genome'].nunique()
        valid = set()
        success = df[df['PCR Product Length'] != "No PCR product"]
        counts = success['Primer Pair'].value_counts()
        all_pairs = df['Primer Pair'].unique()
        for p in all_pairs:
            hits = counts.get(p, 0)
            xr = (hits / total) * 100.0
            if xr <= max_xr: valid.add(p)
            else: print(f"      [Discarded] {p}: {xr:.2f}% cross-reactivity")
        print(f"   [Filter] Retained {len(valid)} pairs.")
        return valid

    def _analyze_results(self, target_csv, non_target_csv, valid_pairs, output_dir):
        df_t = pd.read_csv(target_csv)
        df_nt = pd.read_csv(non_target_csv)
        t_total = df_t['Genome'].nunique()
        nt_total = df_nt['Genome'].nunique()
        summary = []
        for p in valid_pairs:
            t_hits = df_t[(df_t['Primer Pair'] == p) & (df_t['PCR Product Length'] != "No PCR product")]
            sens = (len(t_hits['Genome'].unique()) / t_total) * 100.0
            nt_hits = df_nt[(df_nt['Primer Pair'] == p) & (df_nt['PCR Product Length'] != "No PCR product")]
            xr = (len(nt_hits['Genome'].unique()) / nt_total) * 100.0
            summary.append({"Primer Pair": p, "Sensitivity_Percent": round(sens, 2), "Non-Target_Percent": round(xr, 2)})
        df_sum = pd.DataFrame(summary).sort_values("Sensitivity_Percent", ascending=False)
        df_sum.to_csv(os.path.join(output_dir, "pcr_results_summary.csv"), index=False)
        print(f"   📊 Statistics saved.")