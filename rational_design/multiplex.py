import os
import itertools
import math
import pandas as pd
import json
import csv

MISSING_SEQUENCE_VALUES = {"", "nan", "none", "n/a", "na", "null", "0"}


def clean_sequence(value) -> str:
    """Normalize optional sequence fields from CSV/JSON inputs."""
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except Exception:
        pass
    text = str(value).strip()
    if text.lower() in MISSING_SEQUENCE_VALUES:
        return ""
    return text


def has_valid_probe(assay: dict) -> bool:
    return bool(clean_sequence(assay.get("Probe_Seq")))


# --- 1. BIOPHYSICAL ANALYSIS ALGORITHMS ---

def check_dimer(seq1: str, seq2: str) -> tuple[int, int]:
    """
    Search for the longest consecutive complementary region between seq1 and seq2 (Cross-Dimer/Self-Dimer).
    Particularly checks for complementary pairing at the 3' end as this is where Taq Polymerase acts.
    
    Args:
        seq1 (str): First DNA sequence (e.g., Forward primer).
        seq2 (str): Second DNA sequence (e.g., Reverse primer).
        
    Returns:
        tuple[int, int]: (max_any_match_length, 3_prime_match_length)
    """
    s1 = seq1.upper().strip()
    s2 = seq2.upper().strip()
    if not s1 or not s2:
        return 0, 0
        
    comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    
    n1, n2 = len(s1), len(s2)
    max_any = 0
    
    # 1. Scan for consecutive matches anywhere (Sliding window)
    for shift in range(-n2 + 1, n1):
        consec = 0
        for i in range(n1):
            j = i - shift
            if 0 <= j < n2:
                c1 = s1[i]
                c2 = s2[j]
                is_comp = False
                if c1 in 'ATGC' and c2 in 'ATGC':
                    if (c1 == 'A' and c2 == 'T') or (c1 == 'T' and c2 == 'A') or \
                       (c1 == 'G' and c2 == 'C') or (c1 == 'C' and c2 == 'G'):
                        is_comp = True
                
                if is_comp:
                    consec += 1
                    max_any = max(max_any, consec)
                else:
                    consec = 0

    # 2. Scan for complementary pairing directly at the 3' end (3'-end complementarity)
    # Taq Polymerase extends the primer from the 3' end (last base of the primer)
    consec_3prime = 0
    for k in range(1, min(n1, n2) + 1):
        c1 = s1[n1 - k]
        c2 = s2[n2 - k]
        is_comp = False
        if c1 in 'ATGC' and c2 in 'ATGC':
            is_comp = (c1 == 'A' and c2 == 'T') or (c1 == 'T' and c2 == 'A') or \
                      (c1 == 'G' and c2 == 'C') or (c1 == 'C' and c2 == 'G')
        if is_comp:
            consec_3prime += 1
        else:
            break
            
    return max_any, consec_3prime

def check_hairpin(seq: str) -> int:
    """
    Scan for self-hairpin loop structures.
    Finds complementary stems of size >= 3 bp separated by a loop of 3 to 10 bp.
    
    Args:
        seq (str): DNA sequence to check.
        
    Returns:
        int: Maximum length of the complementary stem forming a hairpin. Returns 0 if no hairpin found.
    """
    s = seq.upper().strip()
    n = len(s)
    if n < 8:
        return 0
        
    max_hairpin = 0
    comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    
    # Loop size from 3 to 10
    for loop in range(3, 11):
        # Maximum possible stem size
        for stem in range(3, (n - loop) // 2 + 1):
            for i in range(n - 2 * stem - loop + 1):
                st1 = s[i : i + stem]
                st2 = s[i + stem + loop : i + 2 * stem + loop]
                
                st2_rc = "".join([comp.get(c, c) for c in reversed(st2)])
                
                if st1 == st2_rc:
                    max_hairpin = max(max_hairpin, stem)
                    
    return max_hairpin

# --- 2. CROSS-TARGET COMPATIBILITY PRE-CHECK ---

def check_cross_target_compatibility(
    new_assay: dict,
    existing_assays: list[dict],
    assay_type: str = "qPCR",
    tm_span_threshold: float = 4.0
) -> dict:
    """
    Check whether a newly accepted assay (from a new target) is compatible
    with all previously accepted assays from earlier targets in the multiplex run.

    Called by the pipeline BEFORE writing ACCEPT for a target, allowing
    REJECT_AND_CONTINUE if the selected primers are incompatible with existing ones.

    Args:
        new_assay (dict): Assay dict with keys: Fwd_Primer, Rev_Primer, Probe_Seq, Target_Name.
        existing_assays (list[dict]): All previously accepted assay dicts from earlier targets.
        assay_type (str): "qPCR" or "Conventional".
        tm_span_threshold (float): Maximum allowed Tm span (°C) across ALL primers cross-target.

    Returns:
        dict: {
            "compatible": bool,
            "tm_span": float,
            "cross_dimer_issues": [str],
            "tm_issues": [str],
            "all_issues": [str]
        }
    """
    if not existing_assays:
        return {"compatible": True, "tm_span": 0.0, "cross_dimer_issues": [], "tm_issues": [], "all_issues": []}

    from Bio.SeqUtils import MeltingTemp as mt
    from Bio.Seq import Seq

    def get_tm(seq):
        clean = "".join([c for c in str(seq).upper() if c in 'ATGC'])
        if len(clean) < 8:
            return 60.0
        try:
            return mt.Tm_NN(
                clean, nn_table=mt.DNA_NN3,
                Na=50, Mg=3.0, dNTPs=0.8, dnac1=250, dnac2=0
            )
        except Exception:
            return 60.0

    cross_dimer_issues = []
    tm_issues = []

    # Collect sequences for new assay
    new_seqs = [("Fwd", new_assay.get("Fwd_Primer", "")), ("Rev", new_assay.get("Rev_Primer", ""))]
    if assay_type == "qPCR" and has_valid_probe(new_assay):
        new_seqs.append(("Probe", new_assay.get("Probe_Seq", "")))

    # Collect ALL primer sequences (new + existing) for Tm span check
    all_primer_tms = [get_tm(new_assay.get("Fwd_Primer", "")), get_tm(new_assay.get("Rev_Primer", ""))]

    for existing in existing_assays:
        ex_seqs = [("Fwd", existing.get("Fwd_Primer", "")), ("Rev", existing.get("Rev_Primer", ""))]
        if assay_type == "qPCR" and has_valid_probe(existing):
            ex_seqs.append(("Probe", existing.get("Probe_Seq", "")))

        # Add existing primer Tms to combined pool
        all_primer_tms.append(get_tm(existing.get("Fwd_Primer", "")))
        all_primer_tms.append(get_tm(existing.get("Rev_Primer", "")))

        # Check cross-dimer between new assay and each existing assay
        for new_role, new_seq in new_seqs:
            if not new_seq:
                continue
            for ex_role, ex_seq in ex_seqs:
                if not ex_seq:
                    continue
                any_m, three_m = check_dimer(new_seq, ex_seq)
                if three_m >= 4:
                    cross_dimer_issues.append(
                        f"🚨 Cross-Dimer 3': [{new_assay.get('Target_Name','New')}_{new_role}] ↔ "
                        f"[{existing.get('Target_Name','Existing')}_{ex_role}] — {three_m}bp tại đầu 3'"
                    )
                elif any_m >= 7:
                    cross_dimer_issues.append(
                        f"⚠️ Cross-Hybridization: [{new_assay.get('Target_Name','New')}_{new_role}] ↔ "
                        f"[{existing.get('Target_Name','Existing')}_{ex_role}] — {any_m}bp bổ sung nội tâm"
                    )

    # Tm span check across ALL primers combined
    tm_span = max(all_primer_tms) - min(all_primer_tms) if len(all_primer_tms) > 1 else 0.0
    if tm_span > tm_span_threshold:
        tm_issues.append(
            f"⚠️ Tm Span tổng thể vượt ngưỡng multiplex: {tm_span:.1f}°C "
            f"(ngưỡng IDT: ≤{tm_span_threshold}°C). "
            f"Khoảng Tm: {min(all_primer_tms):.1f}°C – {max(all_primer_tms):.1f}°C"
        )

    # Fatal issues: 3' cross-dimers = must reject
    fatal = [i for i in cross_dimer_issues if "🚨" in i]
    all_issues = cross_dimer_issues + tm_issues
    compatible = (len(fatal) == 0 and len(tm_issues) == 0)

    return {
        "compatible": compatible,
        "tm_span": round(tm_span, 1),
        "cross_dimer_issues": cross_dimer_issues,
        "tm_issues": tm_issues,
        "all_issues": all_issues
    }


# --- 3. MULTIPLEX COMBINATORIAL EVALUATION ---

def evaluate_multiplex_kit(assays: list[dict], assay_type: str = "qPCR") -> dict:
    """
    Evaluate the thermodynamic and physical compatibility of a Multiplex Assay combination.
    
    Args:
        assays (list[dict]): List of candidate dictionaries for each target species.
            Each dictionary must contain: Set_ID, Fwd_Primer, Rev_Primer, Probe_Seq, Probe_Tm, Target_Name, Amplicon_Size.
        assay_type (str): "qPCR" (TaqMan Probe) or "Conventional" (Gel Electrophoresis).
        
    Returns:
        dict: A dictionary containing compatibility score, verdict, Tm stats, and warnings.
    """
    target_names = [a.get("Target_Name", "Unknown") for a in assays]
    
    # 1. Evaluate Tm mismatch (Tm Balance)
    # Extract Tms for primers and probes
    fwd_tms = []
    rev_tms = []
    probe_tms = []
    
    from Bio.SeqUtils import MeltingTemp as mt
    def get_tm(seq):
        clean = "".join([c for c in str(seq).upper() if c in 'ATGC'])
        if len(clean) < 8: return 60.0
        try:
            # IDT OligoAnalyzer / ThermoFisher-compatible:
            # Na=50mM, Mg=3.0mM, dNTPs=0.8mM, [Oligo]=250nM, SantaLucia 1998
            return mt.Tm_NN(
                clean,
                nn_table=mt.DNA_NN3,
                Na=50,
                Mg=3.0,
                dNTPs=0.8,
                dnac1=250,
                dnac2=0
            )
        except Exception: return 60.0
        
    for a in assays:
        fwd_tms.append(get_tm(a["Fwd_Primer"]))
        rev_tms.append(get_tm(a["Rev_Primer"]))
        if assay_type == "qPCR" and has_valid_probe(a):
            probe_tms.append(get_tm(a["Probe_Seq"]))
            
    all_primer_tms = fwd_tms + rev_tms
    mean_primer_tm = sum(all_primer_tms) / len(all_primer_tms)
    tm_span = max(all_primer_tms) - min(all_primer_tms)
    
    mean_probe_tm = (sum(probe_tms) / len(probe_tms)) if probe_tms else 0.0
    probe_span = (max(probe_tms) - min(probe_tms)) if probe_tms else 0.0
    
    # Thermodynamic compatibility stats
    tm_penalty = tm_span * 2.0  # Penalty for wide primer Tm span
    
    probe_missing_penalty = 0.0
    probe_warnings = []
    size_warnings = []
    size_penalty = 0.0
    
    if assay_type == "qPCR":
        # qPCR/TaqMan mode cannot be used without a valid probe for every target.
        for a in assays:
            if not has_valid_probe(a):
                probe_missing_penalty += 100.0
                probe_warnings.append(
                    f"qPCR requires a valid probe for {a.get('Target_Name', 'Unknown')} "
                    f"({a.get('Set_ID', 'N/A')}); this kit is not usable as a probe-based assay."
                )
                
        if mean_probe_tm > 0:
            # Ideal Probe Tm should be 5-12°C higher than Primer Tm
            probe_diff = mean_probe_tm - mean_primer_tm
            if not (5.0 <= probe_diff <= 12.0):
                tm_penalty += abs(probe_diff - 8.0) * 1.5
                
        # Enforce small qPCR amplicons
        for a in assays:
            size_val = a.get("Amplicon_Size", 150)
            try:
                size = int(size_val)
            except Exception:
                size = 150
            if size > 250:
                size_warnings.append(
                    f"⚠️ Kích thước Sp PCR của {a['Target_Name']} hơi lớn ({size}bp). qPCR lý tưởng nên ngắn (70-200bp) để đạt hiệu suất nhân bản tối đa."
                )
    else:
        # Gel-based conventional mode: ignore probes (no penalties for missing probes).
        # Enforce product size separation of at least 20-50 bp!
        sizes = []
        for a in assays:
            size_val = a.get("Amplicon_Size", 150)
            if size_val == "N/A" or size_val is None:
                size_val = len(a.get("Amplicon_Sequence", ""))
                if size_val == 0: size_val = 150
            try:
                size_val = int(size_val)
            except Exception:
                size_val = 150
            sizes.append((a["Target_Name"], size_val))
            
        sorted_sizes = sorted(sizes, key=lambda x: x[1])
        for idx in range(len(sorted_sizes) - 1):
            t1, s1 = sorted_sizes[idx]
            t2, s2 = sorted_sizes[idx+1]
            gap = s2 - s1
            if gap < 25:
                size_penalty += (35 - gap) * 6.0 # Severe penalty!
                size_warnings.append(
                    f"🚨 Trùng/Sát kích thước Sp PCR: {t1} ({s1}bp) và {t2} ({s2}bp) cách nhau {gap}bp. Agarose gel điện di yêu cầu cách biệt tối thiểu 20-50bp để phân tách rõ ràng."
                )
            elif gap < 50:
                size_warnings.append(
                    f"⚠️ Cách biệt kích thước Sp PCR hơi hẹp: {t1} ({s1}bp) và {t2} ({s2}bp) cách nhau {gap}bp (có thể cần gel agarose độ phân giải cao)."
                )

    # 2. Evaluate Self-Dimer and Hairpin risks for each primer/probe
    self_dimer_warnings = []
    hairpin_warnings = []
    
    for a in assays:
        roles_to_check = [("Fwd", "Fwd_Primer"), ("Rev", "Rev_Primer")]
        if assay_type == "qPCR" and has_valid_probe(a):
            roles_to_check.append(("Probe", "Probe_Seq"))
            
        for role, seq_key in roles_to_check:
            seq = a.get(seq_key, "")
            if not seq: continue
            
            # Hairpin scan
            hp_stem = check_hairpin(seq)
            if hp_stem >= 4:
                hairpin_warnings.append(f"{a['Target_Name']} ({role}): Trình tự dễ gập Hairpin (stem={hp_stem}bp)")
                
            # Self-dimer scan
            any_m, three_m = check_dimer(seq, seq)
            if three_m >= 4 or any_m >= 7:
                self_dimer_warnings.append(f"{a['Target_Name']} ({role}): Tự bắt cặp Dimer (3'={three_m}bp, any={any_m}bp)")

    # 3. Evaluate Cross-Dimerization between different target species
    cross_dimer_warnings = []
    total_dimer_score = 0
    
    # Scan all cross-interacting pairs between species
    n_assays = len(assays)
    for i in range(n_assays):
        for j in range(i + 1, n_assays):
            a1 = assays[i]
            a2 = assays[j]
            
            # Scan all primer combinations of the two species
            sequences_a1 = [("Fwd", a1["Fwd_Primer"]), ("Rev", a1["Rev_Primer"])]
            if assay_type == "qPCR" and has_valid_probe(a1):
                sequences_a1.append(("Probe", a1["Probe_Seq"]))
            
            sequences_a2 = [("Fwd", a2["Fwd_Primer"]), ("Rev", a2["Rev_Primer"])]
            if assay_type == "qPCR" and has_valid_probe(a2):
                sequences_a2.append(("Probe", a2["Probe_Seq"]))
            
            for role1, seq1 in sequences_a1:
                for role2, seq2 in sequences_a2:
                    any_m, three_m = check_dimer(seq1, seq2)
                    
                    # Calculate cross-dimer penalty score
                    total_dimer_score += (three_m * 3.0 + any_m * 1.0)
                    
                    if three_m >= 4:
                        cross_dimer_warnings.append(
                            f"🚨 Nguy cơ Primer-Dimer 3': [{a1['Target_Name']}_{role1}] bắt cặp [{a2['Target_Name']}_{role2}] liên tiếp {three_m}bp ở đầu 3'."
                        )
                    elif any_m >= 7:
                        cross_dimer_warnings.append(
                            f"⚠️ Bám chéo nội tâm (Internal Cross-hybridization): [{a1['Target_Name']}_{role1}] và [{a2['Target_Name']}_{role2}] bổ sung {any_m}bp."
                        )

    # 4. Calculate final Compatibility Score
    score = 100.0 - tm_penalty - total_dimer_score - (len(hairpin_warnings) * 4.0) - (len(self_dimer_warnings) * 5.0) - probe_missing_penalty - size_penalty
    score = max(0.0, round(score, 1))
    
    blocking_reasons = []
    if assay_type == "qPCR" and probe_warnings:
        blocking_reasons.extend(probe_warnings)
        score = 0.0

    fatal_cross_dimers = [w for w in cross_dimer_warnings if "🚨" in w]
    fatal_size_warnings = [w for w in size_warnings if "🚨" in w]
    blocking_reasons.extend(fatal_cross_dimers)
    blocking_reasons.extend(fatal_size_warnings)

    # 5. Final Verdict
    if blocking_reasons:
        verdict = "POOR/UNUSABLE"
    elif score >= 85.0 and len(cross_dimer_warnings) == 0:
        verdict = "EXCELLENT"
    elif score >= 70.0:
        verdict = "GOOD"
    elif score >= 50.0:
        verdict = "MARGINAL"
    else:
        verdict = "POOR/UNUSABLE"
        
    return {
        "score": score,
        "verdict": verdict,
        "assay_type": assay_type,
        "mean_primer_tm": round(mean_primer_tm, 1),
        "mean_probe_tm": round(mean_probe_tm, 1),
        "primer_tm_span": round(tm_span, 1),
        "probe_tm_span": round(probe_span, 1),
        "hairpin_warnings": hairpin_warnings,
        "self_dimer_warnings": self_dimer_warnings,
        "cross_dimer_warnings": cross_dimer_warnings,
        "size_warnings": size_warnings,
        "probe_warnings": probe_warnings,
        "blocking_reasons": blocking_reasons,
        "has_required_probes": not (assay_type == "qPCR" and probe_warnings),
        "assays": assays
    }

# --- 3. COMBINATORIAL ASSEMBLY PIPELINE ---

def assemble_optimal_multiplex_kits(target_folders: list[str], max_kits: int = 3, assay_type: str = "qPCR") -> dict:
    """
    Scan individual target output folders, extract candidates, and perform optimal binary 
    combinations to generate Multiplex Kits.
    
    Args:
        target_folders (list[str]): List of absolute paths to the individual target run folders.
        max_kits (int): Maximum number of top multiplex kits to return.
        assay_type (str): "qPCR" (TaqMan Probe) or "Conventional" (Gel Electrophoresis).
        
    Returns:
        dict: Dictionary containing success status, top kits, and other metadata.
    """
    target_assays = {}
    
    for folder in target_folders:
        if not os.path.exists(folder):
            continue
        final_csv = os.path.join(folder, "FINAL_ASSAY.csv")
        if not os.path.exists(final_csv):
            continue
            
        # Robust target name extraction
        target_name = None
        
        # 1. Try to read from metadata.json
        metadata_path = os.path.join(folder, "metadata.json")
        if os.path.exists(metadata_path):
            try:
                with open(metadata_path, "r", encoding="utf-8") as f:
                    import json
                    meta_data = json.load(f)
                    if "target_name" in meta_data and meta_data["target_name"]:
                        target_name = meta_data["target_name"]
            except Exception:
                pass

        # 2. Try to parse from pipeline_log.txt
        if not target_name:
            log_path = os.path.join(folder, "pipeline_log.txt")
            if os.path.exists(log_path):
                try:
                    with open(log_path, "r", encoding="utf-8") as f:
                        for line in f:
                            if "Auto-detecting genome size for:" in line:
                                parts = line.split("Auto-detecting genome size for:")
                                query = parts[1].strip().rstrip(".")
                                if " AND " in query:
                                    query = query.split(" AND ")[0]
                                if "[Organism]" in query:
                                    query = query.replace("[Organism]", "")
                                query = query.strip('"\'')
                                if query:
                                    target_name = query
                                    break
                except Exception:
                    pass

        # 3. Fallback: Parse from folder structure
        if not target_name:
            basename = os.path.basename(folder)
            parent_name = os.path.basename(os.path.dirname(folder))
            grandparent_name = os.path.basename(os.path.dirname(os.path.dirname(folder)))
            
            # If parent starts with run_ or is generic timestamp, look at grandparent
            if parent_name.startswith("run_") or (len(parent_name) == 15 and parent_name.replace("_", "").isdigit()):
                if grandparent_name and grandparent_name != "runs":
                    if grandparent_name.startswith("multiplex_") or grandparent_name.startswith("local_multiplex_"):
                        target_name = basename.upper() # e.g. T1, T2
                    else:
                        target_name = grandparent_name
                else:
                    target_name = basename.upper()
            else:
                target_name = basename
                    
        target_name = target_name.replace("_", " ").title()
        
        try:
            df = pd.read_csv(final_csv)
            if df.empty: continue
            
            # Extract top 5 candidates for this target to avoid combinatorial explosion
            df_top = df.head(5)
            candidates = []
            seen_assay_keys = set()
            duplicate_count = 0
            for _, row in df_top.iterrows():
                fwd = row.get("Fwd_Primer", "")
                fwd = clean_sequence(fwd)
                rev = row.get("Rev_Primer", "")
                rev = clean_sequence(rev)
                probe = row.get("Probe_Seq", "")
                probe = clean_sequence(probe)
                assay_key = (fwd.upper(), rev.upper(), probe.upper())
                if assay_key in seen_assay_keys:
                    duplicate_count += 1
                    continue
                seen_assay_keys.add(assay_key)
                
                probe_tm = row.get("Probe_Tm", 0.0)
                probe_tm = 0.0 if pd.isna(probe_tm) else float(probe_tm)
                
                amp_size = row.get("Amplicon_Size", 150)
                amp_size = 150 if pd.isna(amp_size) else int(amp_size)
                
                candidates.append({
                    "Set_ID": row.get("Set_ID", "N/A"),
                    "Sensitivity": row.get("Sensitivity", 100.0),
                    "Spec": row.get("Spec", 100.0),
                    "Fwd_Primer": fwd,
                    "Rev_Primer": rev,
                    "Probe_Seq": probe,
                    "Probe_Tm": probe_tm,
                    "Target_Gene": row.get("Target_Gene", "N/A"),
                    "Amplicon_Size": amp_size,
                    "Target_Name": target_name
                })
            if duplicate_count:
                print(f"   ℹ️ Removed {duplicate_count} duplicate assay(s) for {target_name}.")
            if not candidates:
                continue
            target_assays[target_name] = candidates
        except Exception:
            pass
            
    if len(target_assays) < 2:
        return {"error": "Cần tối thiểu 2 target loài đạt kết quả chạy đơn để tổ hợp Multiplex Kit."}
        
    # Generate Cartesian product combinations of assays from each target
    target_names = list(target_assays.keys())
    lists_to_combine = [target_assays[name] for name in target_names]
    
    all_combos = list(itertools.product(*lists_to_combine))
    
    evaluated_kits = []
    for combo in all_combos:
        kit_details = evaluate_multiplex_kit(list(combo), assay_type=assay_type)
        evaluated_kits.append(kit_details)
        
    # Sort by compatibility score descending
    evaluated_kits.sort(key=lambda x: x["score"], reverse=True)
    top_kits = evaluated_kits[:max_kits]
    usable_kits = [
        kit for kit in evaluated_kits
        if kit.get("verdict") in ("EXCELLENT", "GOOD", "MARGINAL") and not kit.get("blocking_reasons")
    ]
    if any(kit.get("verdict") in ("EXCELLENT", "GOOD") for kit in usable_kits):
        overall_verdict = "ACCEPT"
    elif usable_kits:
        overall_verdict = "MARGINAL"
    else:
        overall_verdict = "REJECT"

    blocking_reasons = []
    for kit in top_kits:
        for reason in kit.get("blocking_reasons", []):
            if reason not in blocking_reasons:
                blocking_reasons.append(reason)
    
    return {
        "success": True,
        "total_combinations": len(all_combos),
        "usable_kit_count": len(usable_kits),
        "has_usable_kit": len(usable_kits) > 0,
        "overall_verdict": overall_verdict,
        "blocking_reasons": blocking_reasons,
        "target_species": target_names,
        "assay_type": assay_type,
        "top_kits": top_kits
    }
