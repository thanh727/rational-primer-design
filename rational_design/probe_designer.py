import os
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

from .insilico_pcr_advanced import IndustrialEngine, is_match
from .prober import ProbeSelector, _gc_percent_degenerate, _has_homopolymer
from .multiplex import check_dimer, check_hairpin

def design_probes_for_primers(primers_csv, target_dir, bg_dir=None, output_csv="designed_probes.csv", max_error=4, max_len=1500, probe_max_error=2):
    # 1. Load primer pairs
    df_primers = pd.read_csv(primers_csv)
    primers = df_primers.to_dict('records')
    
    # 2. Build task list of genomes
    tasks = []
    if target_dir and os.path.exists(target_dir):
        for f in os.listdir(target_dir):
            if f.lower().endswith((".fasta", ".fna", ".fa")):
                tasks.append((f, os.path.join(target_dir, f), True))
                
    if bg_dir and os.path.exists(bg_dir):
        for f in os.listdir(bg_dir):
            if f.lower().endswith((".fasta", ".fna", ".fa")):
                tasks.append((f, os.path.join(bg_dir, f), False))
                
    print(f"🧬 Loaded {len(primers)} primer pairs.")
    print(f"📊 Scanning {len(tasks)} genomes using Industrial PCR Engine...")
    
    selector = ProbeSelector()
    results = []
    
    for p in primers:
        name = p['name']
        fwd = p['fwd'].strip().upper()
        rev = p['rev'].strip().upper()
        print(f"▶️ Designing probe for marker: {name}...")
        
        engine = IndustrialEngine(name, fwd, rev, max_error=max_error, max_p=max_len, extract_seq=True, probe_max_error=probe_max_error)
        
        target_amps = []
        bg_amps = []
        total_targets = 0
        total_bgs = 0
        target_hits = 0
        bg_hits = 0
        
        # Run process_genome on all tasks
        for task in tasks:
            filename, filepath, is_target = task
            res = engine.process_genome(task)
            status = res[3]
            amp_seq = res[9]
            if is_target:
                total_targets += 1
                if status == "Positive" and amp_seq != "N/A":
                    target_hits += 1
                    target_amps.append(amp_seq)
            else:
                total_bgs += 1
                if status == "Positive" and amp_seq != "N/A":
                    bg_hits += 1
                    bg_amps.append(amp_seq)
                    
        sens_pct = round((target_hits / max(1, total_targets)) * 100, 1)
        spec_pct = round((1 - (bg_hits / max(1, total_bgs))) * 100, 1) if total_bgs > 0 else 100.0
        
        if not target_amps:
            print(f"   ⚠️ No positive target amplicons detected for {name}. Skipping probe design.")
            results.append({
                "name": name,
                "fwd": fwd,
                "rev": rev,
                "probe": "N/A",
                "probe_tm": 0.0,
                "sensitivity": f"{sens_pct}%",
                "specificity": f"{spec_pct}%"
            })
            continue
            
        p_tm_base = (selector._calc_tm(fwd) + selector._calc_tm(rev)) / 2

        # Helper to check if a probe matches a sequence
        def matches_sequence(probe_seq, seq):
            probe_len = len(probe_seq)
            seq_len = len(seq)
            if seq_len < probe_len:
                return False
            probe_rc = str(Seq(probe_seq).reverse_complement()).upper()
            # Slide window
            for i in range(seq_len - probe_len + 1):
                sub = seq[i : i + probe_len]
                err_fwd = sum(1 for j in range(probe_len) if not is_match(probe_seq[j], sub[j]))
                err_rev = sum(1 for j in range(probe_len) if not is_match(probe_rc[j], sub[j]))
                if err_fwd <= 2 or err_rev <= 2:
                    return True
            return False

        best_probe = _find_probe_for_fixed_primers(target_amps, bg_amps, fwd, rev, p_tm_base)
        
        if best_probe:
            print(f"   🎯 Found Probe: {best_probe['seq']} (Tm: {round(best_probe['tm'], 1)}°C, GC: {round(best_probe['gc'], 1)}%)")
            
            # Recalculate qPCR-level sensitivity and specificity
            qpcr_target_hits = sum(1 for amp in target_amps if matches_sequence(best_probe['seq'], amp))
            qpcr_bg_hits = sum(1 for amp in bg_amps if matches_sequence(best_probe['seq'], amp))
            
            sens_pct_qpcr = round((qpcr_target_hits / max(1, total_targets)) * 100, 1)
            spec_pct_qpcr = round((1 - (qpcr_bg_hits / max(1, total_bgs))) * 100, 1) if total_bgs > 0 else 100.0
            
            results.append({
                "name": name,
                "fwd": fwd,
                "rev": rev,
                "probe": best_probe['seq'],
                "probe_tm": round(best_probe['tm'], 2),
                "sensitivity": f"{sens_pct_qpcr}%",
                "specificity": f"{spec_pct_qpcr}%"
            })
        else:
            print(f"   ❌ No probe found satisfying thermodynamic/specificity rules for {name}.")
            results.append({
                "name": name,
                "fwd": fwd,
                "rev": rev,
                "probe": "None Found",
                "probe_tm": 0.0,
                "sensitivity": f"{sens_pct}%",
                "specificity": f"{spec_pct}%"
            })
            
    # Save output CSV
    pd.DataFrame(results).to_csv(output_csv, index=False)
    print(f"✨ Probe design completed. Report saved at: {output_csv}")
    return True


def _find_probe_for_fixed_primers(target_amps, bg_amps, fwd, rev, p_tm_base):
    # Dynamic relaxation tiers: (max_target_err, min_gc, max_gc, min_tm_diff, max_tm_diff)
    tiers = [
        (0, 35, 65, 5.0, 12.0), # Tier 0: Strict matching, standard GC/Tm
        (0, 30, 70, 4.0, 15.0), # Tier 1: Strict matching, relaxed GC/Tm
        (1, 30, 70, 4.0, 15.0), # Tier 2: 1 Target mismatch allowed, relaxed GC/Tm
        (2, 25, 75, 3.0, 18.0)  # Tier 3: 2 Target mismatches allowed (IUPAC), highly relaxed GC/Tm
    ]

    # Diagnostic counters
    diag = {
        "total_windows": 0,
        "failed_target_match": 0,
        "failed_bg_match": 0,
        "failed_tm": 0,
        "failed_gc": 0,
        "failed_5g": 0,
        "failed_homopolymer": 0,
        "failed_hairpin_dimer": 0
    }

    found_candidate = None

    for max_target_err, min_gc, max_gc, min_tm_diff, max_tm_diff in tiers:
        candidates = []
        
        # Helper to check match with error limit
        def matches_sequence(probe_seq, seq, max_err):
            probe_len = len(probe_seq)
            seq_len = len(seq)
            if seq_len < probe_len:
                return False
            probe_rc = str(Seq(probe_seq).reverse_complement()).upper()
            for i in range(seq_len - probe_len + 1):
                sub = seq[i : i + probe_len]
                err_fwd = sum(1 for j in range(probe_len) if not is_match(probe_seq[j], sub[j]))
                err_rev = sum(1 for j in range(probe_len) if not is_match(probe_rc[j], sub[j]))
                if err_fwd <= max_err or err_rev <= max_err:
                    return True
            return False

        def count_degenerate(seq):
            return sum(1 for char in seq if char not in "ATGC")

        def _has_degenerate_homopolymer(seq, limit=4):
            from .insilico_pcr_advanced import IUPAC
            resolved = [IUPAC.get(b, {b}) for b in seq.upper()]
            n = len(resolved)
            for i in range(n - limit + 1):
                common = resolved[i]
                for j in range(1, limit):
                    common = common.intersection(resolved[i + j])
                if common:
                    return True
            return False

        ref_amp = min(target_amps, key=len)
        left = len(fwd)
        right = len(rev)
        if len(ref_amp) <= left + right + 20:
            continue
            
        search_zone = ref_amp[left : len(ref_amp) - right]
        selector = ProbeSelector()

        # Search length from 20 to 35 bp
        for length in range(20, 36):
            for i in range(len(search_zone) - length + 1):
                diag["total_windows"] += 1
                sub = search_zone[i : i + length]
                probe_seq = None
                
                # Candidate must match all target amplicons within max_target_err
                matches_all_targets = True
                for amp in target_amps:
                    if not matches_sequence(sub, amp, max_target_err):
                        matches_all_targets = False
                        break
                
                if not matches_all_targets:
                    # If not matching all targets, build IUPAC consensus (only if max_target_err > 0)
                    if max_target_err == 0:
                        diag["failed_target_match"] += 1
                        continue
                    candidate_seqs = []
                    import Levenshtein
                    for amp in target_amps:
                        amp_inner = amp[left : len(amp) - right]
                        ref_inner = ref_amp[left : len(ref_amp) - right]
                        if len(amp_inner) == len(ref_inner):
                            candidate_seqs.append(amp_inner[i : i + length])
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
                    
                    if len(candidate_seqs) == len(target_amps):
                        from .utils import generate_iupac_consensus
                        probe_seq = generate_iupac_consensus(candidate_seqs, max_iupac=2)
                    else:
                        diag["failed_target_match"] += 1
                        continue
                else:
                    probe_seq = sub
                    
                if not probe_seq:
                    diag["failed_target_match"] += 1
                    continue
                    
                # Candidate MUST NOT match any background amplicons
                matches_any_bg = False
                for bg_amp in bg_amps:
                    # Background matching is strict (up to 2 mismatches)
                    if matches_sequence(probe_seq, bg_amp, max_err=2):
                        matches_any_bg = True
                        break
                if matches_any_bg:
                    diag["failed_bg_match"] += 1
                    continue
                    
                # Thermodynamic rules
                tm = selector._calc_tm(probe_seq)
                gc = _gc_percent_degenerate(probe_seq)
                
                # 1. Tm must be in relaxed/strict range above primer Tm
                if not (p_tm_base + min_tm_diff <= tm <= p_tm_base + max_tm_diff):
                    diag["failed_tm"] += 1
                    continue
                # 2. GC content range
                if not (min_gc <= gc <= max_gc):
                    diag["failed_gc"] += 1
                    continue
                # 3. No 5' G (including degenerate bases resolving to G)
                from .insilico_pcr_advanced import IUPAC
                first_bases = IUPAC.get(probe_seq[0], {probe_seq[0]})
                if "G" in first_bases:
                    diag["failed_5g"] += 1
                    continue
                # 4. No homopolymers (using degenerate-aware check)
                if _has_degenerate_homopolymer(probe_seq):
                    diag["failed_homopolymer"] += 1
                    continue
                # 5. Hairpin and dimer checks
                if check_hairpin(probe_seq) >= 4 or check_dimer(probe_seq, fwd)[1] >= 4 or check_dimer(probe_seq, rev)[1] >= 4:
                    diag["failed_hairpin_dimer"] += 1
                    continue
                    
                degen_count = count_degenerate(probe_seq)
                candidates.append({'seq': probe_seq, 'tm': tm, 'gc': gc, 'degen': degen_count})
                
        if candidates:
            # Sort order:
            # 1. Number of degenerate bases (ascending, prioritizing 0)
            # 2. Tm closeness to p_tm_base + 8
            # 3. GC content closeness to 50%
            found_candidate = sorted(candidates, key=lambda x: (x['degen'], abs(x['tm'] - (p_tm_base + 8)), abs(x['gc'] - 50)))[0]
            break

    if found_candidate:
        return found_candidate

    print(f"   ⚠️ Diagnostic summary for failed design:")
    print(f"      - Total candidate windows analyzed: {diag['total_windows']}")
    print(f"      - Rejected due to mismatch on targets: {diag['failed_target_match']}")
    print(f"      - Rejected due to match on background (specificity issue): {diag['failed_bg_match']}")
    print(f"      - Rejected due to Tm rules: {diag['failed_tm']}")
    print(f"      - Rejected due to GC content rules: {diag['failed_gc']}")
    print(f"      - Rejected due to 5' G restriction: {diag['failed_5g']}")
    print(f"      - Rejected due to homopolymer restriction: {diag['failed_homopolymer']}")
    print(f"      - Rejected due to hairpin or dimer/cross-dimer: {diag['failed_hairpin_dimer']}")
            
    return None
