# Technical Guideline: Adaptive Background Filtering for Difficult Close-Taxa Primer Design

## Purpose

Implement a backend-only refinement for Rational Primer Design so the current strict k-mer background exclusion remains unchanged for normal cases, but the system automatically switches to an amplicon-level rescue strategy when strict filtering collapses the candidate pool.

This is needed for difficult close-taxa cases such as **S. parasuis vs S. suis**, where a valid primer pair may be highly specific at the PCR-product level even if one or both single primers appear somewhere in the background genomes.

---

## Core Problem

Current strict logic:

```text
1. Mine elite target k-mers.
2. If a k-mer appears in any background genome, blacklist/remove it.
3. Form primer pairs only from clean k-mers.
```

This works well for distant taxa, but fails for very close taxa.

Observed failure example:

```text
Conservation ≥65%: 2598 elite k-mers
Background: 2590/2598 blacklisted (99.7%) → 8 clean
Pairs: 0 found from 8 anchors
```

The biological issue:

```text
single-primer background hit ≠ valid PCR amplification in background
```

A background false positive should only be called when:

```text
Forward primer + reverse primer
same genome/contig
correct orientation
valid product size
optional probe binds, if probe assay
```

---

## Required Design

Do **not** remove the current strict mode.

Add an adaptive backend mode:

```text
background_filter_mode = auto
```

Allowed values:

```text
strict          use existing single-kmer absence logic only
amplicon_level  always use rescue logic
auto            try strict first, then rescue if candidate pool is too small
```

Default should be:

```text
auto
```

No major UI change is required.

---

## Rescue Trigger

After elite k-mer mining and background screening, compute:

```python
elite_total = len(elite_kmers)
clean_total = len(clean_kmers)
clean_ratio = clean_total / max(elite_total, 1)
```

Trigger rescue mode when either condition is true:

```python
clean_total < strict_min_clean_kmers
clean_ratio < strict_min_clean_ratio
```

Recommended defaults:

```json
{
  "strict_min_clean_kmers": 500,
  "strict_min_clean_ratio": 0.05
}
```

Example:

```text
elite_total = 2598
clean_total = 8
clean_ratio = 0.003
```

Expected decision:

```text
Activate rescue_amplicon_level mode.
```

---

## Exact Backend Change

### Current behavior to preserve in strict mode

Existing logic likely resembles:

```python
clean_kmers = []
blacklisted = []

for kmer in elite_kmers:
    if kmer in background_index:
        blacklisted.append(kmer)
        continue
    clean_kmers.append(kmer)
```

Keep this behavior for strict mode.

### New behavior in auto-rescue mode

If the clean pool is too small, do not use only `clean_kmers`.

Instead:

```text
1. Return to all elite k-mers before blacklist removal.
2. Annotate each k-mer with background hit count/frequency.
3. Penalize background single-primer hits, but do not hard-remove them.
4. Use ranked elite k-mers for primer-pair formation.
5. Reject candidate pairs only if they produce a valid background amplicon.
```

---

## Pseudocode

```python
def select_anchors_for_pairing(
    elite_kmers,
    background_index,
    n_background,
    config,
):
    screened = []

    for kmer_obj in elite_kmers:
        seq = kmer_obj.seq

        bg_count = background_index.count_strains(seq)
        bg_freq = bg_count / max(n_background, 1)

        kmer_obj.background_count = bg_count
        kmer_obj.background_freq = bg_freq
        kmer_obj.is_clean = bg_count == 0

        screened.append(kmer_obj)

    clean_kmers = [x for x in screened if x.is_clean]
    clean_ratio = len(clean_kmers) / max(len(screened), 1)

    if config.background_filter_mode == "strict":
        use_rescue = False
    elif config.background_filter_mode == "amplicon_level":
        use_rescue = True
    else:  # auto
        use_rescue = (
            len(clean_kmers) < config.strict_min_clean_kmers
            or clean_ratio < config.strict_min_clean_ratio
        )

    if not use_rescue:
        return clean_kmers, {
            "mode_used": "strict_single_kmer_absence",
            "elite_total": len(screened),
            "clean_total": len(clean_kmers),
            "clean_ratio": clean_ratio,
        }

    rescue_candidates = []

    for x in screened:
        # Optional safety filter: remove only extremely common background k-mers.
        if x.background_freq > config.rescue_max_single_primer_background_freq:
            continue

        x.rescue_score = (
            100.0 * x.target_freq
            - config.rescue_background_penalty_weight * x.background_freq
            + x.primer_quality_score
        )

        rescue_candidates.append(x)

    rescue_candidates.sort(key=lambda x: x.rescue_score, reverse=True)
    rescue_candidates = rescue_candidates[:config.rescue_top_kmers_for_pairing]

    return rescue_candidates, {
        "mode_used": "rescue_amplicon_level",
        "elite_total": len(screened),
        "clean_total": len(clean_kmers),
        "clean_ratio": clean_ratio,
        "rescue_anchor_total": len(rescue_candidates),
    }
```

Recommended rescue parameters:

```json
{
  "rescue_top_kmers_for_pairing": 5000,
  "rescue_max_single_primer_background_freq": 0.25,
  "rescue_background_penalty_weight": 300.0
}
```

---

## Pair Formation

Do not brute-force all possible k-mer pairs.

Use the existing target-side pair formation logic, but pass `anchors_for_pairing` selected by the function above.

Strict mode:

```text
anchors_for_pairing = clean_kmers
```

Rescue mode:

```text
anchors_for_pairing = ranked elite k-mers with background-hit annotation
```

Pairing should remain constrained by:

```text
same target genome/contig
correct orientation
amplicon_min_size <= product_size <= amplicon_max_size
primer QC score
ranking/top-N limits
```

---

## Background Exclusion in Rescue Mode

In rescue mode, background exclusion must happen after primer-pair formation.

Reject a pair only if it forms a valid background amplicon.

```python
def filter_pairs_by_background_amplicon(
    candidate_pairs,
    background_genomes,
    config,
):
    final_pairs = []

    for pair in candidate_pairs:
        bg_hits = find_valid_background_amplicons(
            forward=pair.forward_seq,
            reverse=pair.reverse_seq,
            genomes=background_genomes,
            min_size=config.amplicon_min_size,
            max_size=config.amplicon_max_size,
            mismatch_policy=config.mismatch_policy,
        )

        pair.background_amplicon_hits = len(bg_hits)

        if len(bg_hits) > config.final_max_background_amplicon_hits:
            continue

        final_pairs.append(pair)

    return final_pairs
```

Recommended default:

```json
{
  "final_max_background_amplicon_hits": 0
}
```

Important:

```text
Do not reject a candidate only because one single primer appears in background.
Reject only if the full F/R pair creates a valid background product.
```

---

## Logging Requirements

### Strict mode selected

```text
Background screen:
  Elite k-mers: 50000
  Clean after background: 41845
  Clean ratio: 83.7%
  Mode selected: strict single-kmer exclusion
```

### Auto-rescue activated

```text
Background screen:
  Elite k-mers: 2598
  Clean after background: 8
  Clean ratio: 0.3%

Strict pool too small: 8 < 500
Auto-rescue activated.

Rescue rule:
  Single-primer background hits are annotated/penalized, not discarded.
  Final specificity will be evaluated at amplicon/product level.
```

### After rescue pairing

```text
Rescue pairing:
  Anchors retained for pairing: 2598
  Target candidate pairs formed: X
  Background amplicon-positive pairs rejected: Y
  Final assay-level clean pairs: Z
```

---

## Final Validation Must Not Change

Final validation still uses the full validation panel.

Use biological strain-level counting:

```text
one genome/sample/strain = one validation unit
contigs are sequence fragments, not independent biological records
```

Final specificity should be based on full assay behavior:

```text
background positive = valid PCR product
not single-primer hit
```

---

## Benchmark Case: Published S. parasuis Primer

Use this published S. parasuis primer pair as a regression benchmark:

```text
Forward: CAACTGCTGGATAGTTTCGG
Reverse: GTCTGGCTGAGCTTAATTGG
```

Observed in-silico validation:

```text
Total genomes: 478
Target: S. parasuis
Background: S. suis
Sensitivity: 98.9%
Specificity: 100.0%
Cross-contamination: 0
Verdict: PASS
```

Expected interpretation:

```text
This primer pair is valid at amplicon level.
It must not be rejected solely because a single primer/k-mer has background hits.
```

Use this as an acceptance test for rescue behavior.

---

## Acceptance Tests

### Test 1: Easy distant taxa

Input condition:

```text
clean_ratio high
clean_kmers above threshold
```

Expected:

```text
mode_used = strict_single_kmer_absence
current behavior preserved
```

### Test 2: Close taxa candidate-pool collapse

Input condition:

```text
elite_kmers = 2598
clean_kmers = 8
clean_ratio = 0.003
```

Expected:

```text
mode_used = rescue_amplicon_level
anchors_for_pairing > 8
background exclusion deferred until pair/product validation
```

### Test 3: Published S. parasuis primer backtrace

Input:

```text
Forward: CAACTGCTGGATAGTTTCGG
Reverse: GTCTGGCTGAGCTTAATTGG
```

Expected:

```text
Target sensitivity approximately 98.9%
Background valid amplicon hits = 0
Verdict = PASS
```

---

## Summary for Implementation

Implement this as a backend-only adaptive fallback.

Do not replace strict mode.

Only activate rescue mode when strict background exclusion leaves too few clean elite k-mers.

Core change:

```text
Before:
  background single-kmer hit → discard

After, in rescue mode:
  background single-kmer hit → annotate/penalize
  background valid amplicon → discard
```

This preserves speed for normal cases and improves biological correctness for very close taxa.
