import re
import random
import Levenshtein
from Bio.Seq import Seq

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

def find_primer_hits_new(sequence_s, primer_s, max_mm):
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

def find_primer_hits_old(sequence_s, primer_s, max_mm):
    hits = []
    primer_s = str(primer_s).upper()
    sequence_s = str(sequence_s).upper()
    n = len(primer_s)
    if n == 0 or len(sequence_s) < n:
        return hits
    for i in range(len(sequence_s) - n + 1):
        target_sub = sequence_s[i:i + n]
        mm_count = _iupac_mismatch_count(primer_s, target_sub)
        if mm_count <= max_mm:
            mapping = "".join([t if p != t else "-" for p, t in zip(primer_s, target_sub)])
            hits.append({"pos": i, "map": mapping, "mm": mm_count, "seq": target_sub})
    return hits

# Generate random genome
bases = ["A", "C", "G", "T"]
genome = "".join(random.choice(bases) for _ in range(50000))

# Test with 500 random primers (some matching, some not, some degenerate)
iupac_bases = list(IUPAC.keys())
errors = 0
for i in range(100):
    length = random.randint(18, 24)
    # 50% degenerate
    if random.random() < 0.5:
        primer = "".join(random.choice(iupac_bases) for _ in range(length))
    else:
        primer = "".join(random.choice(bases) for _ in range(length))
    
    max_mm = random.randint(0, 4)
    
    old_hits = find_primer_hits_old(genome, primer, max_mm)
    new_hits = find_primer_hits_new(genome, primer, max_mm)
    
    if old_hits != new_hits:
        print(f"Mismatch found for primer {primer} (max_mm={max_mm}):")
        print("Old:", old_hits)
        print("New:", new_hits)
        errors += 1
        if errors > 5:
            break

if errors == 0:
    print("SUCCESS: 1000 random trials matched perfectly between old and new methods!")
else:
    print(f"FAILED: Found {errors} mismatches.")
