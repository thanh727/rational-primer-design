import multiprocessing
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
from collections import defaultdict, Counter, deque
import pandas as pd
import time
import gc
import array
import os
import math
import sys

# --- SYSTEM OPTIMIZATION ---
gc.set_threshold(100, 10, 10)
WORKER_GENOMES = None

def init_worker(genome_tuples):
    """Initialize shared data for worker processes."""
    global WORKER_GENOMES
    WORKER_GENOMES = genome_tuples

DNA_MAP = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'a': 0, 'c': 1, 'g': 2, 't': 3}


def is_low_complexity_sequence(seq):
    """Reject short tandem repeats such as ATAT... or ATGCATGC...."""
    seq = str(seq).upper()
    if not seq:
        return True
    for unit in range(1, min(5, len(seq) // 2 + 1)):
        if len(seq) % unit == 0 and seq == seq[:unit] * (len(seq) // unit):
            return True
    return False

def is_bio_safe(h, k):
    """Check PCR primer standards on bit-packed h."""
    gc_count = 0
    consecutive = 1
    last_base = -1
    temp_h = h

    for i in range(k):
        base = temp_h & 3
        if base == 1 or base == 2:
            gc_count += 1
        if base == last_base:
            consecutive += 1
            if consecutive >= 5: return False
        else:
            consecutive = 1
        last_base = base
        if i == 0 and (base == 0 or base == 3):
            return False
        temp_h >>= 2

    if not (0.4 * k <= gc_count <= 0.65 * k):
        return False
    return True

def calculate_entropy(h, k):
    """Calculate Shannon Entropy."""
    temp_h = h
    bases = [0, 0, 0, 0] # A, C, G, T
    for _ in range(k):
        bases[temp_h & 3] += 1
        temp_h >>= 2

    entropy = 0
    for c in bases:
        if c > 0:
            p = c / k
            entropy -= p * math.log2(p)
    return entropy

# --- MAIN WORKER: INTRA-STRAIN PAIRING (STRIDE OPTIMIZED) ---
def intra_strain_pairing_worker(args):
    """
    Intra-strain process (Stage 2 - Phase 1):
    Focuses on counting and collecting valid K-mers, applying per-strain filters.
    """
    sid, sequences, k, step, min_copy = args
    valid_hashes = set()
    hash_counts = Counter()
    mask = (1 << (2 * k)) - 1
    step = max(1, int(step)) # Ensure step is at least 1

    for seq in sequences:
        seq = seq.upper()
        n = len(seq)

        # Stride processing
        for i in range(0, n - k + 1, step):
            sub = seq[i : i + k]

            h = 0
            for char in sub:
                val = DNA_MAP.get(char, -1)
                if val == -1:
                    h = -1
                    break
                h = (h << 2) | val
            if h == -1: continue
            h &= mask

            # Filter valid primers & Per-strain filter
            if is_bio_safe(h, k) and calculate_entropy(h, k) >= 1.6 and not is_low_complexity_sequence(sub):
                if min_copy > 1:
                    hash_counts[h] += 1
                    if hash_counts[h] >= min_copy:
                        valid_hashes.add(h)
                else:
                    valid_hashes.add(h)

    return array.array('Q', sorted(list(valid_hashes)))

class PrimerDesigner:
    def __init__(self, params):
        self.params = params
        self.k_min = int(params.get("primer_length_min", 18))
        self.k_max = int(params.get("primer_length_max", 22))
        self.max_candidates = int(params.get("design_max_candidates", 50))
        self.min_p = int(params.get("product_size_min", 70))
        self.max_p = int(params.get("product_size_max", 350))
        self.cpu = int(params.get("cpu_cores", 0)) or max(1, multiprocessing.cpu_count())
        # IDT/ThermoFisher-compatible Tm bounds (user-configurable from GUI)
        self.tm_min = float(params.get("primer_tm_min", 55.0))
        self.tm_max = float(params.get("primer_tm_max", 68.0))
        self.tm_delta_max = 3.0  # IDT standard: Fwd/Rev Tm difference <= 3°C

    def design(self, target_input, bg_input, output_csv):
        print(f"--- 🚀 V-EXTREME DESIGNER: Population Genomics Mode ---")
        t_start = time.time()

        # Get step from config, default is 3 (3x speedup)
        k_step = int(self.params.get("kmer_step", 3))

        target_strains = self._load_fastas(target_input)
        bg_strains = self._load_fastas(bg_input)

        # Auto-adjust k_step based on data size to prevent freezing
        total_bases = sum(len(seq) for seqs in target_strains.values() for seq in seqs)

        target_max_p = 0

        all_raw_pairs = []
        for current_k in range(self.k_min, self.k_max + 1):
            self.k = current_k
            print(f"   --- Primer size k={self.k} ---")

            # PHASE 1: MINING TARGET PAIRS
            print(f"   [Phase 1] Mining Target Pairs ({len(target_strains)} strains) with Step={k_step}...")
            target_pool = Counter()
            min_copy = int(self.params.get("min_intra_strain_copy", 1))

            with multiprocessing.Pool(self.cpu) as pool:
                tasks = [(sid, seqs, self.k, k_step, min_copy) for sid, seqs in target_strains.items()]
                try:
                    for arr in pool.imap_unordered(intra_strain_pairing_worker, tasks):
                        target_pool.update(arr)
                except KeyboardInterrupt:
                    pool.terminate()
                    print("\n   ⚠️ Design interrupted by user.")
                    sys.exit(1)

            base_cons_ratio = self.params.get("design_min_conservation", 0.75)
            min_cons = int(len(target_strains) * base_cons_ratio)

            # Filter k-mers based on initial threshold
            valid_kmers = {h: count for h, count in target_pool.items() if count >= min_cons}

            # Dynamic Filtering: Limit to max 50,000 k-mers to save RAM & CPU
            MAX_ELITE = 50000
            if len(valid_kmers) > MAX_ELITE:
                # Sort by frequency descending
                sorted_kmers = sorted(valid_kmers.items(), key=lambda x: x[1], reverse=True)
                elite_targets = set(h for h, count in sorted_kmers[:MAX_ELITE])
            else:
                elite_targets = set(valid_kmers.keys())

            print(f"      > Pushing {len(elite_targets)} elite K-mers to Phase 2.")

            # Explicit memory cleanup
            del target_pool
            del valid_kmers
            gc.collect()

            # PHASE 2: BACKGROUND SCRUBBING
            if bg_strains and elite_targets:
                bg_blacklist = set()
                elite_frozen = frozenset(elite_targets)

                with multiprocessing.Pool(self.cpu) as pool:
                    bg_step = max(1, k_step // 2)
                    tasks = [(sid, seqs, self.k, elite_frozen, bg_step) for sid, seqs in bg_strains.items()]
                    try:
                        for arr in pool.imap_unordered(self._targeted_kmer_worker, tasks):
                            bg_blacklist.update(arr)
                    except KeyboardInterrupt:
                        pool.terminate()
                        print("\n   ⚠️ Design interrupted by user.")
                        sys.exit(1)

                elite_targets = elite_targets - bg_blacklist
                print(f"      > Remaining {len(elite_targets)} clean primers.")

            # PHASE 3: FINAL SELECTION
            if elite_targets:
                raw_pairs = self._extract_raw_pairs(elite_targets, target_strains)
                all_raw_pairs.extend(raw_pairs)

        if all_raw_pairs:
            print(f"   [Phase 3] Scoring & Ranking {len(all_raw_pairs)} candidate pairs...")
            self._score_and_save(all_raw_pairs, target_strains, output_csv)
            print(f"    ✅ Stage 2 Complete: {time.time()-t_start:.1f}s")
        else:
            print("    ⚠️ Could not find valid primers for any length.")

    def _targeted_kmer_worker(self, args):
        """Scan Background with stride for speedup."""
        sid, sequences, k, elite_set, step = args
        found_hashes = set()
        mask = (1 << (2 * k)) - 1
        for seq in sequences:
            seq_u = seq.upper()
            n = len(seq_u)
            for i in range(0, n - k + 1, step):
                sub = seq_u[i : i + k]
                h = 0
                for char in sub:
                    val = DNA_MAP.get(char, -1)
                    if val == -1:
                        h = -1
                        break
                    h = (h << 2) | val
                if h != -1 and (h & mask) in elite_set:
                    found_hashes.add(h & mask)
        return array.array('Q', sorted(list(found_hashes)))

    def _load_fastas(self, input_path):
        """Load FASTA data using the project contract: 1 FASTA file = 1 strain.

        Important:
        - When input_path is a directory, every FASTA file inside that directory is
          treated as one strain, regardless of how many contigs/records are inside
          the file. The strain ID is therefore the FASTA filename stem.
        - Record IDs inside each FASTA are treated only as contig IDs and must NOT
          redefine strain identity. This avoids collapsing many files into one
          strain when contig headers contain a shared prefix or a pipe symbol.
        - When input_path is a single FASTA file (e.g., the merged output from
          LibraryConstructor), record IDs follow the "strain|contig" convention
          produced by the constructor. These are split by '|' to recover per-strain
          identity. Files without '|' in record IDs are treated as one strain
          (filename = strain ID).
        """
        data = defaultdict(list)
        fasta_exts = ('.fasta', '.fa', '.fna', '.fas')

        if not input_path or not os.path.exists(input_path):
            print(f"   ⚠️ FASTA input not found: {input_path}")
            return data

        input_path = os.path.abspath(str(input_path))

        if os.path.isdir(input_path):
            files = sorted(
                os.path.join(input_path, f)
                for f in os.listdir(input_path)
                if f.lower().endswith(fasta_exts)
            )

            print(
                f"   📂 Loading FASTA folder: {input_path} "
                f"({len(files)} FASTA file(s); 1 file = 1 strain)"
            )

            for f in files:
                sid = os.path.splitext(os.path.basename(f))[0]
                contigs = 0
                total_bp = 0
                try:
                    for r in SeqIO.parse(f, "fasta"):
                        seq = str(r.seq).upper()
                        if not seq:
                            continue
                        data[sid].append(seq)
                        contigs += 1
                        total_bp += len(seq)
                except Exception as e:
                    print(f"      ❌ Error parsing {os.path.basename(f)}: {e}")
                    continue

                if contigs == 0:
                    data.pop(sid, None)
                    print(f"      ⚠️ Empty FASTA skipped: {os.path.basename(f)}")
                else:
                    print(
                        f"      ✅ {sid}: loaded 1 strain from {os.path.basename(f)} "
                        f"({contigs} contig(s), {total_bp/1_000_000:.2f} Mb)."
                    )

            print(f"   ✅ Detected {len(data)} strain(s) from folder: {input_path}")
            return data

        # Single FASTA file input (e.g., merged output from constructor).
        # The constructor always formats IDs as "strain|contig", so we split
        # by '|' to recover per-strain identity. If no '|' is present, the
        # entire file is treated as one strain (filename = strain ID).
        default_sid = os.path.splitext(os.path.basename(input_path))[0]

        print(f"   📄 Loading single FASTA file: {input_path}")

        try:
            for r in SeqIO.parse(input_path, "fasta"):
                seq = str(r.seq).upper()
                if not seq:
                    continue
                sid = r.id.split('|', 1)[0] if r.id and '|' in r.id else default_sid
                data[sid].append(seq)
        except Exception as e:
            print(f"      ❌ Error parsing {os.path.basename(input_path)}: {e}")
            return data

        print(f"   ✅ Detected {len(data)} strain(s) from single FASTA input: {os.path.basename(input_path)}")
        return data

    def _extract_raw_pairs(self, elite_hashes, target_strains):
        ref_sid = list(target_strains.keys())[0]
        ref_seqs = target_strains[ref_sid]
        raw_pairs, seen_pos, pos_gap = [], set(), 150
        mask = (1 << (2 * self.k)) - 1

        for seq in ref_seqs:
            anchors = []
            for i in range(len(seq) - self.k + 1):
                sub = seq[i : i + self.k]
                h = 0
                valid = True
                for char in sub:
                    if char in DNA_MAP: h = (h << 2) | DNA_MAP[char]
                    else: {valid := False}; break
                if valid and (h & mask) in elite_hashes:
                    anchors.append((i, sub))

            for i in range(len(anchors)):
                pos_a, seq_a = anchors[i]
                if any(abs(pos_a - p) < pos_gap for p in seen_pos): continue
                for j in range(i + 1, len(anchors)):
                    pos_b, seq_b = anchors[j]
                    dist = pos_b - pos_a + self.k
                    if self.min_p <= dist <= self.max_p:
                        f, r = seq_a, str(Seq(seq_b).reverse_complement())
                        tm_f = self._calc_tm(f)
                        tm_r = self._calc_tm(r)
                        # IDT/ThermoFisher-compatible pre-filter:
                        # 1) Both primers must be within configured Tm range
                        # 2) Tm_Delta (|Tm_F - Tm_R|) must be <= 3°C
                        if (self.tm_min <= tm_f <= self.tm_max and
                            self.tm_min <= tm_r <= self.tm_max and
                            abs(tm_f - tm_r) <= self.tm_delta_max and
                            self._passes_structural_qc(f, r)):
                            raw_pairs.append({'Fwd': f, 'Rev': r, 'Pos': pos_a})
                            seen_pos.add(pos_a)
                            break
                    elif dist > self.max_p: break
                if len(raw_pairs) >= self.max_candidates: break
        return raw_pairs

    def _score_and_save(self, all_raw_pairs, target_strains, output_csv):
        if not all_raw_pairs: return
        flat_tg = [(sid, s) for sid, seqs in target_strains.items() for s in seqs]
        with multiprocessing.Pool(self.cpu, initializer=init_worker, initargs=(flat_tg,)) as p:
            try:
                results = p.map(self._final_score_worker, all_raw_pairs)
            except KeyboardInterrupt:
                p.terminate()
                print("\n   ⚠️ Design interrupted by user.")
                sys.exit(1)

        final_list = []
        for i, (p, res) in enumerate(zip(all_raw_pairs, results)):
            final_list.append({
                "Set_ID": f"Set_{i+1}", "Forward Primer": p['Fwd'], "Reverse Primer": p['Rev'],
                "Prevalence": f"{res[0]:.1f}%", "Avg_Copy_Number": f"{res[1]:.2f}"
            })
        df = pd.DataFrame(final_list).sort_values("Avg_Copy_Number", ascending=False)
        df.head(self.max_candidates).to_csv(output_csv, index=False)

    def _final_score_worker(self, pair):
        global WORKER_GENOMES
        f, r_rc = pair['Fwd'], str(Seq(pair['Rev']).reverse_complement())
        hits = defaultdict(int)
        for sid, seq in WORKER_GENOMES:
            hits[sid] += min(seq.count(f), seq.count(r_rc))
        pos_strains = sum(1 for h in hits.values() if h > 0)
        return (pos_strains/len(set(s[0] for s in WORKER_GENOMES)))*100, sum(hits.values())/pos_strains if pos_strains > 0 else 0

    def _passes_structural_qc(self, fwd, rev):
        try:
            from .multiplex import check_dimer, check_hairpin
            for seq in (fwd, rev):
                if is_low_complexity_sequence(seq):
                    return False
                if check_hairpin(seq) >= 4:
                    return False
                any_match, three_prime = check_dimer(seq, seq)
                if three_prime >= 4 or any_match >= 8:
                    return False
            any_match, three_prime = check_dimer(fwd, rev)
            return three_prime < 4 and any_match < 8
        except Exception:
            return True

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
