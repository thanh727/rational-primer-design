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

from .validator import find_primer_hits

# --- SYSTEM OPTIMIZATION ---
gc.set_threshold(100, 10, 10)
WORKER_GENOMES = None

def init_worker(genome_tuples):
    """Initialize shared data for worker processes."""
    global WORKER_GENOMES
    WORKER_GENOMES = genome_tuples

DNA_MAP = bytearray(256)
DNA_MAP[:] = b'\xff' * 256
DNA_MAP[65] = DNA_MAP[97] = 0   # A/a
DNA_MAP[67] = DNA_MAP[99] = 1   # C/c
DNA_MAP[71] = DNA_MAP[103] = 2  # G/g
DNA_MAP[84] = DNA_MAP[116] = 3  # T/t


def _is_degenerate(seq):
    """Check if a sequence contains IUPAC degenerate codes."""
    return any(base not in "ATGC" for base in str(seq).upper())


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
    """Check PCR primer standards on bit-packed h (LSB = 3' end)."""
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
        temp_h >>= 2

    if not (0.4 * k <= gc_count <= 0.75 * k):
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

# --- MAIN WORKER: INTRA-STRAIN PAIRING (EXHAUSTIVE SCAN) ---
def intra_strain_pairing_worker(args):
    """
    Intra-strain process (Stage 2 - Phase 1):
    Focuses on counting and collecting valid K-mers, applying per-strain filters.
    """
    sid, sequences, k, step, min_copy = args
    valid_hashes = set()
    hash_counts = Counter()
    mask = (1 << (2 * k)) - 1

    total_windows = 0
    for seq in sequences:
        seq = seq.upper()
        n = len(seq)

        # Exhaustive scan (stride=1)
        for i in range(0, n - k + 1):
            total_windows += 1
            sub = seq[i : i + k]

            h = 0
            for char in sub:
                val = DNA_MAP[ord(char)]
                if val == 0xFF:
                    h = -1
                    break
                h = (h << 2) | val
            if h == -1: continue
            h &= mask

            # Filter valid primers & Per-strain filter
            if is_bio_safe(h, k) and calculate_entropy(h, k) >= 1.4 and not is_low_complexity_sequence(sub):
                if min_copy > 1:
                    hash_counts[h] += 1
                    if hash_counts[h] >= min_copy:
                        valid_hashes.add(h)
                else:
                    valid_hashes.add(h)

    return (sid, array.array('Q', sorted(list(valid_hashes))), total_windows)

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

        target_strains = self._load_fastas(target_input)
        bg_strains = self._load_fastas(bg_input)

        total_bases = sum(len(seq) for seqs in target_strains.values() for seq in seqs)
        print(f"   📊 Target: {len(target_strains)} strains, {total_bases/1e6:.1f} Mb")
        if bg_strains:
            bg_bases = sum(len(seq) for seqs in bg_strains.values() for seq in seqs)
            print(f"   📊 Background: {len(bg_strains)} strains, {bg_bases/1e6:.1f} Mb")
        print()

        base_cons_ratio = self.params.get("design_min_conservation", 0.75)
        min_copy = int(self.params.get("min_intra_strain_copy", 1))

        all_raw_pairs = []
        for current_k in range(self.k_min, self.k_max + 1):
            self.k = current_k
            k_t = time.time()
            print(f"   ┌─ k={self.k} ─────────────────────────────┐")

            # ── PHASE 1 ──
            print(f"   │ [Phase 1] Mining {len(target_strains)} strains (stride=1, exhaustive)")
            target_pool = Counter()
            total_windows_all = 0

            with multiprocessing.Pool(self.cpu) as pool:
                tasks = [(sid, seqs, self.k, 1, min_copy) for sid, seqs in target_strains.items()]
                try:
                    for result in pool.imap_unordered(intra_strain_pairing_worker, tasks):
                        sid, arr, n_windows = result
                        target_pool.update(arr)
                        total_windows_all += n_windows
                except KeyboardInterrupt:
                    pool.terminate()
                    print("\n   ⚠️ Design interrupted by user.")
                    sys.exit(1)
                except (AssertionError, OSError, EOFError):
                    print("\n   ⚠️ Multiprocessing error (connection issue). Trying to continue...")
                    pool.terminate()

            min_cons = int(len(target_strains) * base_cons_ratio)
            before_cons = len(target_pool)

            degenerate = self.params.get("degenerate_primers", False)
            valid_kmers = {}
            for h, count in target_pool.items():
                if count >= min_cons:
                    valid_kmers[h] = count
                elif degenerate:
                    family = count
                    for pos in range(self.k):
                        shift = 2 * pos
                        mask = 3 << shift
                        base = (h >> shift) & 3
                        for alt in range(4):
                            if alt != base:
                                vh = (h & ~mask) | (alt << shift)
                                if vh in target_pool:
                                    family += target_pool[vh]
                    if family >= min_cons:
                        valid_kmers[h] = count

            # Dynamic Filtering: Limit to max 50,000 k-mers
            MAX_ELITE = 50000
            if len(valid_kmers) > MAX_ELITE:
                sorted_kmers = sorted(valid_kmers.items(), key=lambda x: x[1], reverse=True)
                elite_targets = {h: count for h, count in sorted_kmers[:MAX_ELITE]}
            else:
                elite_targets = dict(valid_kmers)

            print(f"   │   {total_windows_all:,} windows scanned → {before_cons:,} unique valid k-mers")
            print(f"   │   Conservation ≥{base_cons_ratio:.0%} ({min_cons}/{len(target_strains)} strains): {len(elite_targets)} elite k-mers")
            del target_pool, valid_kmers
            gc.collect()

            # ── PHASE 2 ──
            anchors_for_pairing = None
            bg_info = None

            if bg_strains and elite_targets:
                bg_hit_counter = Counter()
                elite_frozen = frozenset(elite_targets.keys())

                with multiprocessing.Pool(self.cpu) as pool:
                    bg_step = 1
                    tasks = [(sid, seqs, self.k, elite_frozen, bg_step) for sid, seqs in bg_strains.items()]
                    try:
                        for arr in pool.imap_unordered(self._targeted_kmer_worker, tasks):
                            bg_hit_counter.update(arr)
                    except KeyboardInterrupt:
                        pool.terminate()
                        print("\n   ⚠️ Design interrupted by user.")
                        sys.exit(1)
                    except (AssertionError, OSError, EOFError):
                        print("\n   ⚠️ Background scan multiprocessing error. Continuing with partial data...")
                        pool.terminate()

                anchors_for_pairing, bg_info = self._select_anchors_for_pairing(
                    elite_dict=elite_targets,
                    bg_hit_counter=bg_hit_counter,
                    n_background=len(bg_strains),
                    n_target=len(target_strains),
                )

                self._log_background_screen(bg_info, bg_hit_counter, elite_targets)
            else:
                print(f"   │   Background: none")
                if elite_targets:
                    anchors_for_pairing = set(elite_targets.keys())
                    bg_info = {"mode_used": "strict_single_kmer_absence"}

            # ── PHASE 3 ──
            if anchors_for_pairing:
                raw_pairs = self._extract_raw_pairs(anchors_for_pairing, target_strains)

                use_rescue = bg_info and bg_info.get("mode_used") == "rescue_amplicon_level"
                if use_rescue:
                    n_before = len(raw_pairs)
                    raw_pairs = self._filter_pairs_by_background_amplicon(raw_pairs, bg_strains)
                    n_rejected = n_before - len(raw_pairs)
                    print(f"   │   Rescue pairing:")
                    print(f"   │     Anchors retained for pairing: {bg_info.get('rescue_anchor_total', len(anchors_for_pairing))}")
                    print(f"   │     Target candidate pairs formed: {n_before}")
                    print(f"   │     Background amplicon-positive pairs rejected: {n_rejected}")
                    print(f"   │     Final assay-level clean pairs: {len(raw_pairs)}")
                else:
                    print(f"   │   Pairs: {len(raw_pairs)} found from {len(anchors_for_pairing)} anchors")

                all_raw_pairs.extend(raw_pairs)
            else:
                print(f"   │   Pairs: none (no clean primers)")

            print(f"   └─ {time.time()-k_t:.1f}s{' ' * 25}┘")
            print()

        if all_raw_pairs:
            print(f"   [Final] Scoring & Ranking {len(all_raw_pairs)} candidate pairs across all k...")
            self._score_and_save(all_raw_pairs, target_strains, output_csv)
            print(f"   ✅ Stage 2 Complete: {time.time()-t_start:.1f}s")
        else:
            print("   ⚠️  FAILED: No valid assays found. Consider lowering conservation or expanding k range.")

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
                    val = DNA_MAP[ord(char)]
                    if val == 0xFF:
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
                if f.lower().endswith(fasta_exts) and not f.startswith(".")
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

    def _select_anchors_for_pairing(self, elite_dict, bg_hit_counter, n_background, n_target):
        config = self.params
        mode = config.get("background_filter_mode", "auto")
        strict_min_clean = int(config.get("strict_min_clean_kmers", 500))
        strict_min_ratio = float(config.get("strict_min_clean_ratio", 0.05))
        rescue_top_n = int(config.get("rescue_top_kmers_for_pairing", 5000))
        rescue_max_freq = float(config.get("rescue_max_single_primer_background_freq", 0.25))
        rescue_penalty = float(config.get("rescue_background_penalty_weight", 300.0))

        clean_count = 0
        for h, cons_count in elite_dict.items():
            if bg_hit_counter.get(h, 0) == 0:
                clean_count += 1

        elite_total = len(elite_dict)
        clean_ratio = clean_count / max(elite_total, 1)

        if mode == "strict":
            use_rescue = False
        elif mode == "amplicon_level":
            use_rescue = True
        else:
            use_rescue = (clean_count < strict_min_clean or clean_ratio < strict_min_ratio)

        if not use_rescue:
            clean_hashes = {h for h in elite_dict if bg_hit_counter.get(h, 0) == 0}
            return clean_hashes, {
                "mode_used": "strict_single_kmer_absence",
                "elite_total": elite_total,
                "clean_total": clean_count,
                "clean_ratio": clean_ratio,
            }

        scored = []
        for h, cons_count in elite_dict.items():
            bg_count = bg_hit_counter.get(h, 0)
            bg_freq = bg_count / max(n_background, 1)
            if bg_freq > rescue_max_freq:
                continue
            target_freq = cons_count / max(n_target, 1)
            score = 100.0 * target_freq - rescue_penalty * bg_freq
            scored.append((h, score))

        scored.sort(key=lambda x: x[1], reverse=True)
        rescue_anchors = {h for h, s in scored[:rescue_top_n]}

        return rescue_anchors, {
            "mode_used": "rescue_amplicon_level",
            "elite_total": elite_total,
            "clean_total": clean_count,
            "clean_ratio": clean_ratio,
            "rescue_anchor_total": len(rescue_anchors),
        }

    def _log_background_screen(self, bg_info, bg_hit_counter, elite_dict=None):
        mode = bg_info.get("mode_used", "strict_single_kmer_absence")
        elite_total = bg_info.get("elite_total", 0)
        clean_total = bg_info.get("clean_total", 0)
        clean_ratio = bg_info.get("clean_ratio", 0)

        print(f"   │   Background screen:")
        print(f"   │     Elite k-mers: {elite_total}")
        print(f"   │     Clean after background: {clean_total}")
        print(f"   │     Clean ratio: {clean_ratio:.1%}")

        if mode == "rescue_amplicon_level":
            strict_min_clean = int(self.params.get("strict_min_clean_kmers", 500))
            print(f"   │")
            print(f"   │   Strict pool too small: {clean_total} < {strict_min_clean}")
            print(f"   │   Auto-rescue activated.")
            print(f"   │")
            print(f"   │   Rescue rule:")
            print(f"   │     Single-primer background hits are annotated/penalized, not discarded.")
            print(f"   │     Final specificity will be evaluated at amplicon/product level.")
        else:
            print(f"   │   Mode selected: strict single-kmer exclusion")

    def _filter_pairs_by_background_amplicon(self, candidate_pairs, bg_strains):
        final_max_hits = int(self.params.get("final_max_background_amplicon_hits", 0))
        max_mm = int(self.params.get("max_mismatch", 2))
        min_p = int(self.params.get("product_size_min", 70))
        max_p = int(self.params.get("product_size_max", 350))

        final_pairs = []
        for pair in candidate_pairs:
            fwd = pair['Fwd']
            rev = pair['Rev']
            rev_rc = str(Seq(rev).reverse_complement()).upper()

            total_bg_amplicons = 0

            for sid, seqs in bg_strains.items():
                for seq_s in seqs:
                    f_hits = find_primer_hits(seq_s, fwd, max_mm)
                    r_hits = find_primer_hits(seq_s, rev_rc, max_mm)

                    r_idx = 0
                    r_len = len(r_hits)
                    for f in f_hits:
                        while r_idx < r_len and r_hits[r_idx]["pos"] <= f["pos"]:
                            r_idx += 1
                        for j in range(r_idx, r_len):
                            r = r_hits[j]
                            p_len = r["pos"] - f["pos"] + len(rev_rc)
                            if p_len > max_p:
                                break
                            if min_p <= p_len:
                                total_bg_amplicons += 1

            if total_bg_amplicons <= final_max_hits:
                final_pairs.append(pair)

        return final_pairs

    def _extract_raw_pairs(self, elite_hashes, target_strains):
        ref_sid = list(target_strains.keys())[0]
        ref_seqs = target_strains[ref_sid]
        raw_pairs, seen_pos, pos_gap = [], set(), 150
        mask = (1 << (2 * self.k)) - 1
        n_anchors = 0
        n_dist_fail = 0
        n_tm_fail = 0
        n_qc_fail = 0

        degenerate = self.params.get("degenerate_primers", False)
        if degenerate:
            variant_hashes = set()
            for h in elite_hashes:
                for pos in range(self.k):
                    shift = 2 * pos
                    base = (h >> shift) & 3
                    for alt in range(4):
                        if alt != base:
                            vh = (h & ~(3 << shift)) | (alt << shift)
                            variant_hashes.add(vh)

        for seq in ref_seqs:
            anchors = []
            n = len(seq)
            if n < self.k:
                continue
            h = 0
            roll_count = 0
            for i in range(n):
                val = DNA_MAP[ord(seq[i])]
                if val != 0xFF:
                    h = ((h << 2) | val) & mask
                    roll_count += 1
                    if roll_count >= self.k:
                        start = i - self.k + 1
                        if h in elite_hashes or (degenerate and h in variant_hashes):
                            anchors.append((start, seq[start:start + self.k]))
                else:
                    roll_count = 0
                    h = 0

            n_anchors += len(anchors)
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
                        if (self.tm_min <= tm_f <= self.tm_max and
                            self.tm_min <= tm_r <= self.tm_max and
                            abs(tm_f - tm_r) <= self.tm_delta_max):
                            if self._passes_structural_qc(f, r):
                                raw_pairs.append({'Fwd': f, 'Rev': r, 'Pos': pos_a})
                                seen_pos.add(pos_a)
                                break
                            else:
                                n_qc_fail += 1
                        else:
                            n_tm_fail += 1
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
            except (AssertionError, OSError, EOFError):
                print("\n   ⚠️ Scoring multiprocessing error. Saving partial results...")
                p.terminate()
                results = []

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
            if _is_degenerate(fwd) or _is_degenerate(rev):
                return True
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
        Conditions: Na=50mM, Mg=3.0mM, dNTPs=0.8mM, [Oligo]=250nM.
        For degenerate sequences, Tm is approximated from non-IUPAC bases."""
        clean = "".join([c for c in str(seq).upper() if c in 'ATGC'])
        if len(clean) < 4:
            return 0.0
        return mt.Tm_NN(
            Seq(clean),
            nn_table=mt.DNA_NN3,
            Na=50,
            Mg=3.0,
            dNTPs=0.8,
            dnac1=250,
            dnac2=0
        )
