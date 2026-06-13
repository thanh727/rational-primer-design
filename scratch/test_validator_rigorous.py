import re
import time
import random
import tracemalloc
import Levenshtein
from Bio.Seq import Seq

# IUPAC matching rules
IUPAC = {
    "A": {"A"}, "C": {"C"}, "G": {"G"}, "T": {"T"},
    "R": {"A", "G"}, "Y": {"C", "T"}, "S": {"G", "C"},
    "W": {"A", "T"}, "K": {"G", "T"}, "M": {"A", "C"},
    "B": {"C", "G", "T"}, "D": {"A", "G", "T"},
    "H": {"A", "C", "T"}, "V": {"A", "C", "G"},
    "N": {"A", "C", "G", "T"},
}

IUPAC_REGEX = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "[AG]", "Y": "[CT]", "S": "[GC]",
    "W": "[AT]", "K": "[GT]", "M": "[AC]",
    "B": "[CGT]", "D": "[AGT]", "H": "[ACT]",
    "V": "[ACG]", "N": "[ACGT]",
}

def _iupac_to_regex(iupac_seq):
    return "".join(IUPAC_REGEX.get(c, c) for c in iupac_seq.upper())

def _iupac_mismatch_count(primer_s, target_s):
    return sum(1 for p, t in zip(primer_s, target_s) if t not in IUPAC.get(p, {p}))

def _primer_segments(primer_s, max_mm):
    n = len(primer_s)
    segment_count = max(1, min(max_mm + 1, n))
    base_len, extra = divmod(n, segment_count)
    offset = 0
    for idx in range(segment_count):
        seg_len = base_len + (1 if idx < extra else 0)
        if seg_len > 0:
            yield offset, primer_s[offset:offset + seg_len]
        offset += seg_len

# OPTIMIZED regex-based pigeonhole search
def find_primer_hits_optimized(sequence_s, primer_s, max_mm):
    hits = []
    primer_s = str(primer_s).upper()
    sequence_s = str(sequence_s).upper()
    n = len(primer_s)
    if n == 0 or len(sequence_s) < n:
        return hits
    candidate_starts = set()
    for offset, segment in _primer_segments(primer_s, max_mm):
        pattern = _iupac_to_regex(segment)
        for match in re.finditer(f"(?=({pattern}))", sequence_s):
            start_pos = match.start() - offset
            if 0 <= start_pos <= len(sequence_s) - n:
                candidate_starts.add(start_pos)
    for start_pos in sorted(candidate_starts):
        target_sub = sequence_s[start_pos:start_pos + n]
        mm_count = _iupac_mismatch_count(primer_s, target_sub)
        if mm_count <= max_mm:
            mapping = "".join([t if p != t else "-" for p, t in zip(primer_s, target_sub)])
            hits.append({"pos": start_pos, "map": mapping, "mm": mm_count, "seq": target_sub})
    return hits

# OLD brute force window scanning search (when degenerate)
def find_primer_hits_old(sequence_s, primer_s, max_mm):
    hits = []
    primer_s = str(primer_s).upper()
    sequence_s = str(sequence_s).upper()
    n = len(primer_s)
    if n == 0 or len(sequence_s) < n:
        return hits
    # Brute-force window scanning
    for i in range(len(sequence_s) - n + 1):
        target_sub = sequence_s[i:i + n]
        mm_count = _iupac_mismatch_count(primer_s, target_sub)
        if mm_count <= max_mm:
            mapping = "".join([t if p != t else "-" for p, t in zip(primer_s, target_sub)])
            hits.append({"pos": i, "map": mapping, "mm": mm_count, "seq": target_sub})
    return hits

# Setup benchmark environment
print("⚙️ Generating mock 200,000 bp genome sequence...")
bases = ["A", "C", "G", "T"]
genome = "".join(random.choice(bases) for _ in range(200000))

# Test definitions
test_cases = [
    # (name, sequence, max_mismatch)
    ("Standard Primer (Exact)", "TACGTATCGATCGATCGATC", 0),
    ("Standard Primer (2 Mismatches)", "TACGTATCGATCGATCGATC", 2),
    ("Standard Primer (4 Mismatches)", "TACGTATCGATCGATCGATC", 4),
    ("IUPAC Primer (1 Degenerate base, 1 Mismatch)", "WGAACCTGTATTCACAAGYG", 1),
    ("IUPAC Primer (Multi degenerate, 3 Mismatches)", "CYYCCATCTAAGGCTAGTTT", 3),
    ("Highly Degenerate Primer (5 Mismatches)", "RYNSRYNSATCGATCGNRYN", 5),
]

print("\n🚀 RUNNING COMPARISON & BENCHMARK:")
print("-" * 105)
print(f"{'Test Case':<45} | {'Old Hits':<8} | {'New Hits':<8} | {'Status':<10} | {'Old Time':<10} | {'New Time':<10}")
print("-" * 105)

for name, primer, max_mm in test_cases:
    # Benchmark Old
    tracemalloc.start()
    t0 = time.perf_counter()
    old_hits = find_primer_hits_old(genome, primer, max_mm)
    t_old = time.perf_counter() - t0
    _, peak_old = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    # Benchmark New
    tracemalloc.start()
    t0 = time.perf_counter()
    new_hits = find_primer_hits_optimized(genome, primer, max_mm)
    t_new = time.perf_counter() - t0
    _, peak_new = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    # Check validity
    status = "✅ MATCH" if old_hits == new_hits else "❌ MISMATCH"
    
    print(f"{name:<45} | {len(old_hits):<8} | {len(new_hits):<8} | {status:<10} | {t_old:.4f}s | {t_new:.4f}s")
    
    if old_hits != new_hits:
        print(f"   ⚠️ WARNING: Discrepancy found in hits details!")

print("\n🔥 RUNNING MEMORY/RAM STRESS TEST (Batch of 15 Primers on 200k bp genome):")

primers_batch = []
for _ in range(15):
    length = random.randint(18, 22)
    # Mix of degenerate and standard
    if random.random() < 0.5:
        primers_batch.append(("".join(random.choice(list(IUPAC.keys())) for _ in range(length)), random.randint(0, 3)))
    else:
        primers_batch.append(("".join(random.choice(bases) for _ in range(length)), random.randint(0, 3)))

# Measure Old Batch RAM
tracemalloc.start()
t0 = time.perf_counter()
for pr, mm in primers_batch:
    _ = find_primer_hits_old(genome, pr, mm)
time_old_batch = time.perf_counter() - t0
_, peak_old_batch = tracemalloc.get_traced_memory()
tracemalloc.stop()

# Measure New Batch RAM
tracemalloc.start()
t0 = time.perf_counter()
for pr, mm in primers_batch:
    _ = find_primer_hits_optimized(genome, pr, mm)
time_new_batch = time.perf_counter() - t0
_, peak_new_batch = tracemalloc.get_traced_memory()
tracemalloc.stop()

print(f"Old Method (Brute force scanning): Time = {time_old_batch:.3f}s, Peak RAM = {peak_old_batch / (1024 * 1024):.2f} MB")
print(f"New Method (Optimized regex-pigeon): Time = {time_new_batch:.3f}s, Peak RAM = {peak_new_batch / (1024 * 1024):.2f} MB")
speedup = time_old_batch / time_new_batch
ram_reduction = peak_old_batch / peak_new_batch if peak_new_batch > 0 else 1
print(f"⚡ Speedup: {speedup:.1f}x faster")
print(f"💾 Peak RAM Reduction: {ram_reduction:.1f}x lower")
