import multiprocessing
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
from collections import defaultdict, Counter
import csv
import primer3
from pathlib import Path
import re
import gc
import time # <--- Added

# ... (Include all the helper functions: init_worker, get_canonical_2bit, mine_binary_worker, check_bg_binary_worker, fast_hamming, find_fuzzy_occurrence, validate_candidate_fuzzy from previous answer) ...
# (I am omitting them here for brevity, assume they are the same as the "2-Bit Turbo" version)

# --- GLOBAL SHARED MEMORY ---
WORKER_GENOMES = None

def init_worker(genomes_list):
    global WORKER_GENOMES
    WORKER_GENOMES = genomes_list

TRANS_TABLE = str.maketrans("ACGT", "0123")
RC_TABLE = str.maketrans("ACGT", "TGCA")

def get_canonical_2bit(kmer_str):
    fwd = kmer_str
    rev = fwd.translate(RC_TABLE)[::-1]
    target = fwd if fwd < rev else rev
    try: return int(target.translate(TRANS_TABLE), 4)
    except ValueError: return -1

def mine_binary_worker(args):
    seq_indices, k = args
    global WORKER_GENOMES
    counts = defaultdict(int)
    if not WORKER_GENOMES: return counts
    for idx in seq_indices:
        seq = WORKER_GENOMES[idx]
        n = len(seq)
        if n < k: continue
        unique_in_genome = set()
        for i in range(n - k + 1):
            k_str = seq[i : i+k]
            val = get_canonical_2bit(k_str)
            if val != -1: unique_in_genome.add(val)
        for val in unique_in_genome:
            counts[val] += 1
    return counts

def check_bg_binary_worker(args):
    bg_seqs_chunk, k, candidates_set = args
    bad_counts = defaultdict(int)
    for seq in bg_seqs_chunk:
        n = len(seq)
        if n < k: continue
        unique_bg = set()
        for i in range(n - k + 1):
            k_str = seq[i : i+k]
            val = get_canonical_2bit(k_str)
            if val in candidates_set: unique_bg.add(val)
        for val in unique_bg: bad_counts[val] += 1
    return bad_counts

def fast_hamming(s1, s2, max_mm):
    mm = 0
    for i in range(len(s1)):
        if s1[i] != s2[i]:
            mm += 1
            if mm > max_mm: return False
    return True

def find_fuzzy_occurrence(genome_seq, primer, max_mm):
    n = len(primer)
    half = n // 2
    part1 = primer[:half]
    part2 = primer[half:]
    start = 0
    while True:
        idx = genome_seq.find(part1, start)
        if idx == -1: break
        sub = genome_seq[idx : idx+n]
        if len(sub) == n and fast_hamming(sub, primer, max_mm): return True
        start = idx + 1
    start = 0
    while True:
        idx = genome_seq.find(part2, start)
        if idx == -1: break
        full_start = idx - half
        if full_start >= 0:
            sub = genome_seq[full_start : full_start+n]
            if len(sub) == n and fast_hamming(sub, primer, max_mm): return True
        start = idx + 1
    return False

def validate_candidate_fuzzy(pair_data):
    global WORKER_GENOMES
    if not WORKER_GENOMES: return 0.0
    fwd = pair_data['Fwd']
    rev = pair_data['Rev']
    rev_rc = str(Seq(rev).reverse_complement())
    MAX_MM = 2
    hits = 0
    total = len(WORKER_GENOMES)
    for t_seq in WORKER_GENOMES:
        if find_fuzzy_occurrence(t_seq, fwd, MAX_MM):
            if find_fuzzy_occurrence(t_seq, rev_rc, MAX_MM):
                hits += 1
    return (hits / total) * 100.0

class PrimerDesigner:
    def __init__(self, params):
        self.params = params
        self.k = params.get("primer_length", 20)
        self.max_candidates = params.get("design_max_candidates", 20)
        self.cpu = params.get("cpu_cores", 1)

    def design(self, target_fasta, bg_fasta, output_csv):
        print(f"--- STARTING PRIMER DESIGN (2-Bit Integer Mode) ---")
        
        t_start = time.time()
        print("   [IO] Loading Target Genomes into RAM...")
        targets = [str(r.seq).upper() for r in SeqIO.parse(target_fasta, "fasta")]
        if not targets: return
        print(f"   [IO] Loaded in {time.time()-t_start:.1f}s")
        
        # 2. Mine
        t0 = time.time()
        min_prev = self.params.get("design_min_conservation", 0.90)
        conserved_ints = self._mine_integers(targets, min_prev)
        print(f"   [Step 1 Time] {time.time()-t0:.1f}s")
        
        if not conserved_ints: return

        # 3. Filter
        t0 = time.time()
        max_bg = self.params.get("design_max_bg_prevalence", 0.20)
        valid_ints = self._filter_integers(conserved_ints, bg_fasta, max_bg)
        print(f"   [Step 2 Time] {time.time()-t0:.1f}s")
        
        if not valid_ints: return

        # 4. Recover
        t0 = time.time()
        print("   [Step 3] Converting Integers back to DNA...")
        valid_seqs = self._recover_sequences(valid_ints, targets)
        print(f"   [Step 3 Time] {time.time()-t0:.1f}s")

        # 5. Pair
        t0 = time.time()
        self._pair_and_validate(valid_seqs, targets, output_csv)
        print(f"   [Step 4 Time] {time.time()-t0:.1f}s")

    # ... (Keep methods _mine_integers, _filter_integers, _recover_sequences, _pair_and_validate EXACTLY as in previous answer) ...
    
    def _mine_integers(self, targets, min_prevalence):
        print(f"   [Step 1] Mining Conserved Integers (Length: {self.k}bp)...")
        total_records = len(targets)
        indices = list(range(total_records))
        chunk_size = max(1, total_records // (self.cpu * 4))
        batches = [indices[i:i + chunk_size] for i in range(0, len(indices), chunk_size)]
        final_counts = defaultdict(int)
        with multiprocessing.Pool(self.cpu, initializer=init_worker, initargs=(targets,)) as pool:
            tasks = [(b, self.k) for b in batches]
            results = pool.map(mine_binary_worker, tasks)
            for res in results:
                for kmer_int, count in res.items():
                    final_counts[kmer_int] += count
        threshold = int(total_records * min_prevalence)
        candidates = {k for k, v in final_counts.items() if v >= threshold}
        print(f"   [Result] {len(candidates)} unique integer candidates found.")
        return candidates

    def _filter_integers(self, candidates_set, bg_file, max_bg_prev):
        print(f"   [Step 2] Filtering against Background (Integer Mode)...")
        bg_chunks = []
        chunk = []
        for r in SeqIO.parse(bg_file, "fasta"):
            chunk.append(str(r.seq).upper())
            if len(chunk) >= 50:
                bg_chunks.append(chunk); chunk = []
        if chunk: bg_chunks.append(chunk)
        total_bg = sum(len(c) for c in bg_chunks)
        if total_bg == 0: return candidates_set
        bad_counts = defaultdict(int)
        tasks = [(c, self.k, candidates_set) for c in bg_chunks]
        with multiprocessing.Pool(self.cpu) as pool:
            results = pool.map(check_bg_binary_worker, tasks)
            for res in results:
                for kmer_int, count in res.items(): bad_counts[kmer_int] += count
        limit = int(total_bg * max_bg_prev)
        final_set = {k for k in candidates_set if bad_counts[k] <= limit}
        print(f"   [Result] {len(final_set)} integers survived background filter.")
        return final_set

    def _recover_sequences(self, valid_ints, targets):
        recovered = []
        needed = set(valid_ints)
        found_map = {}
        for seq in targets[:50]:
            if len(found_map) == len(needed): break
            n = len(seq)
            if n < self.k: continue
            for i in range(n - self.k + 1):
                k_str = seq[i : i+self.k]
                val = get_canonical_2bit(k_str)
                if val in needed and val not in found_map:
                    found_map[val] = k_str
        return list(found_map.values())

    def _calculate_tm(self, seq_str):
        try: return mt.Tm_NN(Seq(seq_str), nn_table=mt.DNA_NN3, Na=50, Mg=1.5, dNTPs=0.2, dnac1=500, dnac2=500)
        except: return mt.Tm_NN(Seq(seq_str))

    def _is_good_sequence(self, seq_str, tm):
        gc_pct = (seq_str.count('G') + seq_str.count('C')) / len(seq_str) * 100
        if gc_pct < 20.0 or gc_pct > 80.0: return False
        if re.search(r'A{6,}|T{6,}|G{6,}|C{6,}', seq_str): return False
        return True

    def _pair_and_validate(self, valid_oligos, target_seqs, output_csv):
        valid_set = set(valid_oligos)
        ref_seq = target_seqs[0]
        ref_len = len(ref_seq)
        
        print("   [QC] Checking Thermodynamics...")
        oligo_stats = {}
        opt_tm = self.params.get("primer_opt_tm", 60.0)
        
        for oligo in valid_oligos:
            tm = self._calculate_tm(oligo)
            if abs(tm - opt_tm) <= 6.0:
                if self._is_good_sequence(oligo, tm):
                    oligo_stats[oligo] = tm
        
        if not oligo_stats: return

        print("   [Pairing] Mapping seeds to reference...")
        fwd_candidates = []
        for i in range(ref_len - self.k + 1):
            sub = ref_seq[i : i+self.k]
            if sub in oligo_stats:
                fwd_candidates.append({'start': i, 'seq': sub, 'tm': oligo_stats[sub]})
                
        raw_pairs = []
        min_amp = self.params.get("product_size_min", 70)
        max_amp = self.params.get("product_size_max", 250)
        
        for f in fwd_candidates:
            f_start = f['start']
            for r_end in range(f_start + min_amp, f_start + max_amp):
                if r_end > ref_len: break
                target_site = ref_seq[r_end-self.k : r_end]
                rev_seq = str(Seq(target_site).reverse_complement())
                if rev_seq in oligo_stats:
                    r_tm = oligo_stats[rev_seq]
                    if abs(f['tm'] - r_tm) <= 7.0:
                        raw_pairs.append({
                            'Fwd': f['seq'], 'Rev': rev_seq, 
                            'AmpLen': r_end - f_start, 'FwdStart': f_start,
                            'FwdTm': f['tm'], 'RevTm': r_tm
                        })
        
        opt_len = (min_amp + max_amp) // 2
        raw_pairs.sort(key=lambda x: (abs(x['AmpLen'] - opt_len), abs(x['FwdTm'] - x['RevTm'])))

        print("   [Filtering] Removing redundant pairs (overlap < 5bp)...")
        unique_candidates = []
        seen_positions = []
        for p in raw_pairs:
            p1 = p['FwdStart']
            is_redundant = False
            for accepted_p1 in seen_positions:
                if abs(p1 - accepted_p1) < 5:
                    is_redundant = True; break
            if not is_redundant:
                seen_positions.append(p1)
                unique_candidates.append(p)
                if len(unique_candidates) >= 1000: break
        
        print(f"   [Filtering] Retained {len(unique_candidates)} unique candidates.")

        print(f"   [Turbo] Validating using {self.cpu} cores (Fuzzy)...")
        final_list = []
        design_threshold = self.params.get("design_min_conservation", 0.95) * 100
        
        pool = multiprocessing.Pool(self.cpu, initializer=init_worker, initargs=(target_seqs,))
        try:
            results = pool.map(validate_candidate_fuzzy, unique_candidates)
            for i, prev in enumerate(results):
                if len(final_list) >= self.max_candidates: break
                if prev >= design_threshold:
                    pair = unique_candidates[i]
                    final_list.append({
                        'Set_ID': f"Set_{len(final_list)+1}",
                        'Forward Primer': pair['Fwd'],
                        'Reverse Primer': pair['Rev'],
                        'Prevalence': f"{prev:.1f}%"
                    })
                    print(f"      [+] Pair Validated! Prev: {prev:.1f}% (Fuzzy)")
        finally:
            pool.close(); pool.join()
        
        if final_list:
            with open(output_csv, 'w', newline='') as f:
                dict_writer = csv.DictWriter(f, final_list[0].keys())
                dict_writer.writeheader()
                dict_writer.writerows(final_list)
            print(f"   ✅ [SUCCESS] Saved {len(final_list)} pairs.")
        else:
            print("   ❌ No pairs passed validation.")