import argparse
import sys
import json
import csv
import multiprocessing
import datetime
import time
import random
import gc
import ctypes
from pathlib import Path

from .fetcher import SequenceFetcher
from .constructor import LibraryConstructor
from .designer import PrimerDesigner
from .validator import InSilicoValidator, positive_product_mask
from .prober import ProbeSelector
import math
import shutil
import pandas as pd
from .utils import nuclear_ram_flush, DualLogger, generate_batch_analytical_summary


def load_json(path: str) -> dict:
    """Load JSON data from a file path."""
    with open(path) as f: return json.load(f)


def clean_optional_sequence(value) -> str:
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except Exception:
        pass
    text = str(value).strip()
    if text.lower() in {"", "nan", "none", "n/a", "na", "null", "0"}:
        return ""
    return text


def clean_optional_number(value, default=0.0):
    try:
        if value is None or pd.isna(value):
            return default
        number = float(value)
        if math.isnan(number) or math.isinf(number):
            return default
        return number
    except Exception:
        return default


def load_final_assays_for_audit(path_final_assay: Path) -> list[dict]:
    if not path_final_assay.exists():
        return []
    try:
        df_final = pd.read_csv(path_final_assay)
    except Exception:
        return []

    assays = []
    for _, row in df_final.iterrows():
        assays.append({
            "Set_ID": str(row.get("Set_ID", "N/A")),
            "Fwd_Primer": clean_optional_sequence(row.get("Fwd_Primer", "")),
            "Rev_Primer": clean_optional_sequence(row.get("Rev_Primer", "")),
            "Probe_Seq": clean_optional_sequence(row.get("Probe_Seq", "")) or "N/A",
            "Probe_Tm": clean_optional_number(row.get("Probe_Tm", 0.0)),
            "Sensitivity": str(row.get("Sensitivity", "")),
            "Specificity": str(row.get("Spec", "")),
            "Amplicon_Size": clean_optional_number(row.get("Amplicon_Size", 0), default=None),
            "Target_Gene": str(row.get("Target_Gene", "N/A")),
        })
    return assays


def write_audit_trail(
    base_dir: Path,
    base_params: dict,
    success: bool,
    path_final_assay: Path,
    validation_report: dict | None = None,
) -> None:
    audit_data = {
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "parameters": base_params,
        "success": success,
        "final_assay_file": str(path_final_assay) if success else None,
        "final_assays": load_final_assays_for_audit(path_final_assay) if success else [],
        "ai_report_included": False,
        "ai_decision": None,
        "validation_report": validation_report,
    }
    report_path = base_dir / "ai_expert_report.json"
    if report_path.exists():
        with open(report_path, "r", encoding="utf-8") as rf:
            audit_data["ai_report"] = json.load(rf)
            audit_data["ai_report_included"] = True
            audit_data["ai_decision"] = {
                "overall_verdict": audit_data["ai_report"].get("overall_verdict"),
                "next_action": audit_data["ai_report"].get("next_action"),
                "best_assay_name": audit_data["ai_report"].get("best_assay_name"),
            }

    audit_path = base_dir / "audit_trail.json"
    with open(audit_path, "w", encoding="utf-8") as f_audit:
        json.dump(audit_data, f_audit, ensure_ascii=False, indent=4)
    print(f"🔒 Immutable Audit Trail generated: {audit_path.name}")

def load_config_for_fetcher(json_path: str) -> dict:
    """
    Parse a JSON configuration file for the SequenceFetcher.

    Args:
        json_path (str): Path to the configuration JSON file.

    Returns:
        dict: Parsed dictionary mapping labels to tuple configurations.
    """
    data = load_json(json_path)
    parsed = {}
    for k, v in data.items():
        query = v[0]
        size = v[1]
        count = v[2] if len(v) > 2 else 0
        target_type = v[3] if len(v) > 3 else "genome"
        parsed[k] = (query, size, count, target_type)
    return parsed

def format_duration(seconds: float) -> str:
    """Format duration in seconds into a human-readable string."""
    m, s = divmod(int(seconds), 60)
    h, m = divmod(m, 60)
    if h > 0: return f"{h}h {m}m {s}s"
    return f"{m}m {s}s"

def run_full_pipeline(args: argparse.Namespace) -> None:
    """
    Execute the entire rational primer design pipeline end-to-end.

    Args:
        args (argparse.Namespace): Command line arguments.
    """
    # --- REPRODUCIBILITY ---
    random.seed(getattr(args, 'seed', 42))

    pipeline_start = time.time()

    base_dir = Path(args.out)
    base_dir.mkdir(parents=True, exist_ok=True)

    log_path = base_dir / "pipeline_log.txt"
    sys.stdout = DualLogger(str(log_path))
    sys.stderr = sys.stdout

    print(f"========================================================")
    print(f"   🧬 RATIONAL DESIGN LOG START")
    print(f"   🕒 Time: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"   📂 Project: {base_dir.resolve()}")
    print(f"========================================================")

    # --- DEFAULTS ---
    defaults = {
        "design_target_sampling_size": 0, "design_background_sampling_size": 100,
        "validation_target_sampling_size": 0, "validation_background_sampling_size": 200,
        "design_max_candidates": 50, "cpu_cores": 0,
        "primer_length_min": 18, "primer_length_max": 22,
        "primer_tm_min": 55.0, "primer_tm_max": 68.0,
        "product_size_min": 100, "product_size_max": 350,
        "design_min_conservation": 0.75, "min_sensitivity": 95.0, "max_mismatch": 3,
        "enable_blast": True,
        "auto_relax_constraints": True,
        "degenerate_primers": True,
        "max_iupac_per_primer": 2,
        "background_filter_mode": "auto",
        "strict_min_clean_kmers": 500,
        "strict_min_clean_ratio": 0.05,
        "rescue_top_kmers_for_pairing": 5000,
        "rescue_max_single_primer_background_freq": 0.25,
        "rescue_background_penalty_weight": 300.0,
        "final_max_background_amplicon_hits": 0

    }

    user_params = {}
    config_file = Path(args.params) if args.params else Path("config/parameters.json")
    if config_file.exists():
        print(f"\n[CONFIG] Attempting to load: {config_file}")
        try:
            user_params = load_json(config_file)
            print("   ✅ Configuration loaded.")
        except Exception as e:
            print(f"   ⚠️ Error reading config ({e}). Using INTERNAL DEFAULTS.")

    base_params = {**defaults, **user_params}

    # --- CROSS-TARGET MULTIPLEX CONTEXT (for sequential multiplex runs) ---
    shared_context_path = getattr(args, 'shared_context', None)
    cross_target_context = []
    if shared_context_path and Path(shared_context_path).exists():
        try:
            with open(shared_context_path, "r", encoding="utf-8") as scf:
                cross_target_context = json.load(scf)
            print(f"\n   🔗 Cross-target context loaded: {len(cross_target_context)} earlier target(s) already accepted.")
            for ctx in cross_target_context:
                print(f"      - {ctx.get('target_name', 'Unknown')}: Fwd={ctx.get('Fwd_Primer','?')[:20]}... Tm_Fwd={ctx.get('Fwd_Tm','?'):.1f}°C")
        except Exception as e:
            print(f"   ⚠️ Could not load shared context: {e}")
            cross_target_context = []
    if base_params["cpu_cores"] <= 0:
        base_params["cpu_cores"] = max(1, multiprocessing.cpu_count() - 2)

    # --- PATHS ---
    ws_dir = base_dir / "1_workspace"
    path_design_target = ws_dir / "design" / "target.fasta"
    path_design_bg = ws_dir / "design" / "background.fasta"
    path_val_target = ws_dir / "validate" / "target.fasta"
    path_val_bg = ws_dir / "validate" / "background.fasta"

    cand_dir = base_dir / "2_candidates"
    path_candidates_csv = cand_dir / "candidates.csv"
    path_candidates_fasta = cand_dir / "candidates.fasta"

    val_dir = base_dir / "3_validation"
    path_val_target_csv = val_dir / "results_target.csv"
    path_val_stats_csv = val_dir / "pcr_results_summary.csv"
    path_final_assay = base_dir / "FINAL_ASSAY.csv"

    # ==========================================
    # STAGE 0 & 1: DATA PREP
    # ==========================================
    if getattr(args, 'email', None) and getattr(args, 'target_config', None) and getattr(args, 'bg_config', None):
        print("\n--- [STAGE 0] FETCHING NCBI DATA ---")
        fetcher = SequenceFetcher(args.email)
        raw_dir = base_dir / "0_raw_data"

        t_conf = load_config_for_fetcher(args.target_config)
        b_conf = load_config_for_fetcher(args.bg_config)

        path_raw_target = raw_dir / "target"
        path_raw_bg = raw_dir / "background"

        path_raw_target.mkdir(parents=True, exist_ok=True)
        path_raw_bg.mkdir(parents=True, exist_ok=True)

        fetcher.fetch_and_save_all(t_conf, str(path_raw_target))
        fetcher.fetch_and_save_all(b_conf, str(path_raw_bg))

    elif args.local_target and args.local_bg:
        path_raw_target, path_raw_bg = Path(args.local_target), Path(args.local_bg)
    else: return print("❌ Error: Missing data inputs (need either local files or NCBI config).")

    print("\n--- [STAGE 1] BUILDING DATASETS (Strain-based) ---")
    cons = LibraryConstructor()
    cons.construct(str(path_raw_target), str(path_raw_bg), {
        "design_target": str(path_design_target), "design_background": str(path_design_bg),
        "validation_target": str(path_val_target), "validation_background": str(path_val_bg)
    }, {
        "design_target": base_params['design_target_sampling_size'],
        "design_background": base_params['design_background_sampling_size'],
        "validation_target": base_params['validation_target_sampling_size'],
        "validation_background": base_params['validation_background_sampling_size']
    })
    del cons # Delete constructor immediately
    nuclear_ram_flush()

    # ==========================================
    # AUTO-OPTIMIZATION LOOP
    # ==========================================
    user_max_cand = int(base_params.get("design_max_candidates", 50))
    strategies = [
        {"candidates": user_max_cand, "relax_cons": 0.0, "relax_sens": 0.0}
    ]
    if base_params.get("auto_relax_constraints", False):
        strategies.extend([
            {"candidates": user_max_cand, "relax_cons": 0.05, "relax_sens": 5.0},
            {"candidates": user_max_cand, "relax_cons": 0.10, "relax_sens": 10.0},
        ])

    success = False
    for i, strat in enumerate(strategies):
        current_params = base_params.copy()
        current_params.update({
            'design_max_candidates': strat['candidates'],
            'design_min_conservation': max(0.5, current_params['design_min_conservation'] - strat['relax_cons']),
            'min_sensitivity': current_params['min_sensitivity'] - strat['relax_sens']
        })

        print(f"\n" + "━"*60 + f"\n   🔄 CYCLE {i+1}/{len(strategies)}\n" + "━"*60)

        # --- [STAGE 2] DESIGNING PRIMERS (SCOPED) ---
        print("\n--- [STAGE 2] DESIGNING PRIMERS ---")
        cand_dir.mkdir(parents=True, exist_ok=True)
        if path_candidates_csv.exists(): path_candidates_csv.unlink()

        designer = PrimerDesigner(params=current_params)
        designer.design(str(path_design_target), str(path_design_bg), str(path_candidates_csv))

        # DESTROY DESIGNER IMMEDIATELY
        del designer
        nuclear_ram_flush()

        if not path_candidates_csv.exists(): continue

        # --- [STAGE 3] IN-SILICO VALIDATION (ITERATIVE BATCH) ---
        print("\n--- [STAGE 3] IN-SILICO VALIDATION (ITERATIVE BATCH) ---")
        val_dir.mkdir(parents=True, exist_ok=True)

        batch_size = int(current_params.get("validation_batch_size", 50))
        max_batches = int(current_params.get("validation_max_batches", 3))

        try:
            with open(path_candidates_csv, 'r') as f_in:
                all_candidates = list(csv.DictReader(f_in))
        except Exception as e:
            print(f"   ❌ Error reading candidates: {e}"); continue

        total_batches = min(max_batches, math.ceil(len(all_candidates) / batch_size))
        perfect_assay_found = False

        path_master_target = val_dir / "master_results_target.csv"
        path_master_bg = val_dir / "master_results_background.csv"
        path_master_stats = val_dir / "master_pcr_results_summary.csv"

        for p in [path_master_target, path_master_bg, path_master_stats]:
            if p.exists(): p.unlink()

        def append_csv(source_path, dest_path, is_first_batch):
            if not source_path.exists(): return
            with open(source_path, 'r') as f_src, open(dest_path, 'a') as f_dst:
                if not is_first_batch:
                    next(f_src, None)
                shutil.copyfileobj(f_src, f_dst)

        for batch_idx in range(total_batches):
            batch_candidates = all_candidates[batch_idx * batch_size : (batch_idx + 1) * batch_size]

            with open(path_candidates_fasta, 'w') as f_out:
                for row in batch_candidates:
                    f_out.write(f">{row['Set_ID']}_Fwd\n{row['Forward Primer']}\n")
                    f_out.write(f">{row['Set_ID']}_Rev\n{row['Reverse Primer']}\n")

            print(f"   ▶️ [Batch {batch_idx+1}/{total_batches}] Validating {len(batch_candidates)} candidates...")

            validator = InSilicoValidator()
            validator.validate(str(path_val_target), str(path_val_bg), str(path_candidates_fasta), str(val_dir), current_params)
            del validator
            nuclear_ram_flush()

            append_csv(val_dir / "results_target.csv", path_master_target, batch_idx == 0)
            append_csv(val_dir / "results_background.csv", path_master_bg, batch_idx == 0)
            append_csv(val_dir / "pcr_results_summary.csv", path_master_stats, batch_idx == 0)

            if path_master_stats.exists():
                df_stats = pd.read_csv(path_master_stats)
                perfect_cands = df_stats[df_stats['Sensitivity_Percent'] >= current_params['min_sensitivity']]

                # --- [AI COGNITIVE EVALUATION] AUTONOMOUS LOOP ---
                if getattr(args, 'ai_base_url', None):
                    print("   🤖 Consult AI Core (Autonomous Mode)...")
                    try:
                        from .ai_core import LocalBackend, AssayEvaluator
                        backend = LocalBackend(base_url=args.ai_base_url, model_name=getattr(args, 'ai_model', 'llama3'))
                        evaluator = AssayEvaluator(backend)
                        # Generate comparative analysis report with fully calculated parameters
                        analytical_report = generate_batch_analytical_summary(str(path_master_stats), current_params, language=getattr(args, 'language', 'en'))

                        # Pass cross-target context for Gate 4 evaluation in multiplex mode
                        report = evaluator.evaluate_candidates(
                            analytical_report,
                            language=getattr(args, 'language', 'en'),
                            cross_target_context=cross_target_context if cross_target_context else None
                        )
                        if "error" not in report:
                            print(f"      - Verdict: {report.get('overall_verdict', 'N/A')}")
                            print(f"      - Next Action: {report.get('next_action', 'N/A')}")
                            gate_check = report.get('biophysical_gate_check', {})
                            if gate_check:
                                print(f"      - Gate Check: {gate_check}")
                            print(f"      - Reason: {report.get('clinical_recommendation', 'N/A')}")

                            report_path = base_dir / "ai_expert_report.json"
                            with open(report_path, "w", encoding="utf-8") as rf:
                                json.dump(report, rf, ensure_ascii=False, indent=4)

                            # --- SHARED SUB-BATCH RETRY LOGIC (used by both REJECT and GATE 4 failure) ---
                            def _try_sub_batches() -> bool:
                                """Try remaining validated candidates in sub-batches of 10, up to 3 times. Returns True if accepted."""
                                df_sorted = df_stats.sort_values(by=['Sensitivity_Percent', 'Tm_Delta'], ascending=[False, True])
                                sub_batch_size = 10
                                max_sub_batches = 3

                                for sub_idx in range(max_sub_batches):
                                    start = 15 + sub_idx * sub_batch_size
                                    end = start + sub_batch_size
                                    sub_candidates = df_sorted.iloc[start:end]

                                    if sub_candidates.empty:
                                        print(f"      No more candidates to check after position {start}.")
                                        break

                                    sub_stats_path = val_dir / f"sub_batch_{batch_idx}_{sub_idx}.csv"
                                    sub_candidates.to_csv(sub_stats_path, index=False)
                                    sub_report = generate_batch_analytical_summary(str(sub_stats_path), current_params, language=getattr(args, 'language', 'en'))
                                    if sub_stats_path.exists(): sub_stats_path.unlink()

                                    sub_result = evaluator.evaluate_candidates(
                                        sub_report,
                                        language=getattr(args, 'language', 'en'),
                                        cross_target_context=cross_target_context if cross_target_context else None
                                    )

                                    if "error" not in sub_result and sub_result.get('next_action') == "ACCEPT_AND_STOP":
                                        print(f"      ✅ Sub-batch {sub_idx+1}/{max_sub_batches} ACCEPTED by AI!")
                                        gate4_ok = not cross_target_context
                                        if cross_target_context:
                                            try:
                                                from .multiplex import check_cross_target_compatibility
                                                top_3 = sub_result.get('top_3_assays', [])
                                                df_check = pd.read_csv(path_master_stats)
                                                id_col = 'Set_ID' if 'Set_ID' in df_check.columns else 'Primer Pair'
                                                best_candidate = None
                                                if top_3:
                                                    clean_id = str(top_3[0]).split()[0].strip()
                                                    match = df_check[df_check[id_col] == clean_id]
                                                    if not match.empty:
                                                        row = match.iloc[0]
                                                        best_candidate = {
                                                            "Fwd_Primer": row.get('Forward Primer', ''),
                                                            "Rev_Primer": row.get('Reverse Primer', ''),
                                                            "Probe_Seq": row.get('Probe_Seq', 'N/A'),
                                                            "Target_Name": row.get('Target_Name', 'Current')
                                                        }
                                                if not best_candidate and not df_check.empty:
                                                    row = df_check.iloc[0]
                                                    best_candidate = {
                                                        "Fwd_Primer": row.get('Forward Primer', ''),
                                                        "Rev_Primer": row.get('Reverse Primer', ''),
                                                        "Probe_Seq": row.get('Probe_Seq', 'N/A'),
                                                        "Target_Name": 'Current'
                                                    }
                                                if best_candidate:
                                                    existing_assays = [{
                                                        "Fwd_Primer": c.get("Fwd_Primer", ""),
                                                        "Rev_Primer": c.get("Rev_Primer", ""),
                                                        "Probe_Seq": c.get("Probe_Seq", "N/A"),
                                                        "Target_Name": c.get("target_name", "")
                                                    } for c in (cross_target_context or [])]
                                                    compat = check_cross_target_compatibility(
                                                        best_candidate, existing_assays,
                                                        assay_type="qPCR", tm_span_threshold=4.0
                                                    )
                                                    gate4_ok = compat["compatible"]
                                                    if not gate4_ok:
                                                        print(f"      🚨 Sub-batch GATE 4 FAILED — continuing...")
                                                        for issue in compat["all_issues"]:
                                                            print(f"         {issue}")
                                                    else:
                                                        print(f"      ✅ Sub-batch GATE 4 PASSED (Tm span={compat['tm_span']}°C)")
                                            except Exception as gate4_err:
                                                print(f"      ⚠️ Sub-batch Gate 4 check error (non-fatal): {gate4_err}")

                                        if gate4_ok:
                                            with open(report_path, "w", encoding="utf-8") as rf:
                                                json.dump(sub_result, rf, ensure_ascii=False, indent=4)
                                            return True
                                    else:
                                        r = sub_result.get('clinical_recommendation', sub_result.get('error', 'Unknown')) if "error" not in sub_result else sub_result.get('error')
                                        print(f"      ⚠️ Sub-batch {sub_idx+1}/{max_sub_batches} rejected: {str(r)[:120]}")

                                return False

                            if report.get('next_action') == "ACCEPT_AND_STOP":
                                # --- GATE 4 HARD CHECK (Python-level, bypasses AI hallucination) ---
                                gate4_failed = False
                                if cross_target_context:
                                    try:
                                        from .multiplex import check_cross_target_compatibility
                                        top_3 = report.get('top_3_assays', [])
                                        df_check = pd.read_csv(path_master_stats)
                                        id_col = 'Set_ID' if 'Set_ID' in df_check.columns else 'Primer Pair'
                                        best_candidate = None
                                        if top_3:
                                            clean_id = str(top_3[0]).split()[0].strip()
                                            match = df_check[df_check[id_col] == clean_id]
                                            if not match.empty:
                                                row = match.iloc[0]
                                                best_candidate = {
                                                    "Fwd_Primer": row.get('Forward Primer', ''),
                                                    "Rev_Primer": row.get('Reverse Primer', ''),
                                                    "Probe_Seq": row.get('Probe_Seq', 'N/A'),
                                                    "Target_Name": row.get('Target_Name', 'Current')
                                                }
                                        if not best_candidate and not df_check.empty:
                                            row = df_check.iloc[0]
                                            best_candidate = {
                                                "Fwd_Primer": row.get('Forward Primer', ''),
                                                "Rev_Primer": row.get('Reverse Primer', ''),
                                                "Probe_Seq": row.get('Probe_Seq', 'N/A'),
                                                "Target_Name": 'Current'
                                            }
                                        if best_candidate:
                                            existing_assays = [{
                                                "Fwd_Primer": c.get("Fwd_Primer", ""),
                                                "Rev_Primer": c.get("Rev_Primer", ""),
                                                "Probe_Seq": c.get("Probe_Seq", "N/A"),
                                                "Target_Name": c.get("target_name", "")
                                            } for c in cross_target_context]
                                            compat = check_cross_target_compatibility(
                                                best_candidate, existing_assays,
                                                assay_type="qPCR", tm_span_threshold=4.0
                                            )
                                            if not compat["compatible"]:
                                                print(f"   🚨 GATE 4 FAILED (Python hard check) — Cross-target incompatibility!")
                                                for issue in compat["all_issues"]:
                                                    print(f"      {issue}")
                                                gate4_failed = True
                                            else:
                                                print(f"   ✅ GATE 4 PASSED: Cross-target Tm span={compat['tm_span']}°C, no fatal dimers.")
                                    except Exception as gate4_err:
                                        print(f"   ⚠️ Gate 4 check error (non-fatal): {gate4_err}")

                                if not gate4_failed:
                                    print(f"   ✨ AI ACCEPTED primers. Triggering Early Stopping!")
                                    perfect_assay_found = True
                                    break
                                else:
                                    print(f"   ➡️ GATE 4 forced reject — checking subsequent candidates...")

                                # Fall through to sub-batch retry when GATE 4 fails
                            else:
                                print(f"   ⚠️ AI REJECTED primers. Checking subsequent candidates in sub-batches of 10...")

                            if _try_sub_batches():
                                print(f"   ✨ AI ACCEPTED via sub-batch retry. Triggering Early Stopping!")
                                perfect_assay_found = True
                                break
                        else:
                            print(f"   ❌ AI Error: {report.get('error')}")
                    except Exception as e:
                        print(f"   ❌ AI Integration Error: {e}")

                # If AI is disabled or fails, use Fallback Hard Rules
                if not perfect_assay_found and not perfect_cands.empty:
                    print(f"   ✨ Batch {batch_idx+1} meets technical standards (Hard Rules). Triggering Early Stopping!")
                    perfect_assay_found = True
                    break

        path_val_target_csv = path_master_target
        path_val_stats_csv = path_master_stats

        if not path_val_target_csv.exists(): continue
        if not perfect_assay_found:
            print("   ⚠️ No assays accepted by AI or hard rules across all batches. Skipping probe design.")
            continue

        # --- [STAGE 4 & 5] PROBE DESIGN & ANNOTATION (SCOPED) ---
        print("\n--- [STAGE 4] PROBE DESIGN ---")
        prober = ProbeSelector()
        if path_final_assay.exists(): path_final_assay.unlink()

        prober_params = current_params.copy()
        prober_params["background_results_csv"] = str(path_master_bg)
        found_probes = prober.design(str(path_val_target_csv), str(path_final_assay), str(path_val_stats_csv), prober_params)

        if found_probes:
            print("\n   ✨ SUCCESS! Valid assays generated.")

            # --- FILTER TOP 3 ASSAYS BEFORE BLAST ---
            try:
                df_final = pd.read_csv(str(path_final_assay))

                report_path = base_dir / "ai_expert_report.json"
                top_3_names = []
                if report_path.exists():
                    with open(report_path, "r", encoding="utf-8") as rf:
                        try:
                            ai_rep = json.load(rf)
                            top_3_names = ai_rep.get("top_3_assays", [])
                        except Exception:
                            pass

                if top_3_names and isinstance(top_3_names, list):
                    # AI sometimes appends sequences like 'Set_8 (ACTG/ACTG)'. We need to extract just 'Set_8'.
                    clean_names = [str(n).split()[0].strip() for n in top_3_names]
                    print(f"   🤖 Filtering to Top AI-selected assays: {clean_names}")
                    df_final = df_final[df_final["Set_ID"].isin(clean_names)]
                else:
                    print("   📉 AI Disabled: Exporting ALL valid assays for manual evaluation.")
                    df_final['Sens_num'] = df_final['Sensitivity'].str.rstrip('%').astype('float')
                    df_final = df_final.sort_values(by=['Sens_num', 'Max_Copy'], ascending=[False, False])
                    df_final = df_final.drop(columns=['Sens_num'])

                seq_cols = [c for c in ["Fwd_Primer", "Rev_Primer", "Probe_Seq"] if c in df_final.columns]
                if seq_cols:
                    before_dedup = len(df_final)
                    for col in seq_cols:
                        df_final[col] = df_final[col].apply(clean_optional_sequence)
                    df_final = df_final.drop_duplicates(subset=seq_cols, keep="first")
                    removed = before_dedup - len(df_final)
                    if removed:
                        print(f"   🧹 Removed {removed} duplicate final assay(s) with identical primer/probe sequences.")

                df_final.to_csv(str(path_final_assay), index=False)
            except Exception as e:
                print(f"   ⚠️ Error while filtering Top 3 assays: {e}")

            if current_params.get("enable_blast", False):
                print("\n--- [STAGE 5] BLAST ANNOTATION ---")
                prober.run_blast_annotation(str(path_final_assay))

            # --- WRITE ACCEPTED ASSAYS TO SHARED CONTEXT (Multiplex sequential mode) ---
            if shared_context_path:
                try:
                    from Bio.SeqUtils import MeltingTemp as mt
                    from Bio.Seq import Seq
                    def _calc_tm_shared(seq):
                        clean = "".join([c for c in str(seq).upper() if c in 'ATGC'])
                        if len(clean) < 8: return 0.0
                        try:
                            return round(mt.Tm_NN(
                                Seq(clean), nn_table=mt.DNA_NN3,
                                Na=50, Mg=3.0, dNTPs=0.8, dnac1=250, dnac2=0
                            ), 1)
                        except Exception: return 0.0

                    df_accepted = pd.read_csv(str(path_final_assay))
                    # Extract target name from metadata.json if available
                    meta_path = base_dir / "metadata.json"
                    target_name = base_dir.name
                    if meta_path.exists():
                        try:
                            with open(meta_path, "r") as mf:
                                target_name = json.load(mf).get("target_name", base_dir.name)
                        except Exception:
                            pass

                    # Load existing shared context
                    existing_ctx = []
                    if Path(shared_context_path).exists():
                        try:
                            with open(shared_context_path, "r") as scf:
                                existing_ctx = json.load(scf)
                        except Exception:
                            existing_ctx = []

                    # Append new accepted assays to context
                    for _, acc_row in df_accepted.iterrows():
                        fwd = str(acc_row.get('Fwd_Primer', ''))
                        rev = str(acc_row.get('Rev_Primer', ''))
                        probe = clean_optional_sequence(acc_row.get('Probe_Seq', '')) or "N/A"
                        existing_ctx.append({
                            "target_name": target_name,
                            "Set_ID": str(acc_row.get('Set_ID', 'N/A')),
                            "Fwd_Primer": fwd,
                            "Rev_Primer": rev,
                            "Probe_Seq": probe,
                            "Fwd_Tm": _calc_tm_shared(fwd),
                            "Rev_Tm": _calc_tm_shared(rev),
                            "Probe_Tm": _calc_tm_shared(probe) if probe not in ('N/A', '') else 0.0,
                        })

                    with open(shared_context_path, "w", encoding="utf-8") as scf:
                        json.dump(existing_ctx, scf, ensure_ascii=False, indent=4)
                    print(f"   🔗 Shared context updated: {len(existing_ctx)} assay(s) now in {shared_context_path}")
                except Exception as ctx_err:
                    print(f"   ⚠️ Could not update shared context: {ctx_err}")

            success = True
            break

        del prober
        nuclear_ram_flush()

    if success:
        print(f"\n✅✅✅ PIPELINE COMPLETE! Result: {path_final_assay.resolve()}")
    else:
        print("\n❌ FAILED after all cycles. No valid assays found.")

    # --- [STAGE 6] SAVE VALIDATION REPORT + CROSS-CONTAMINATION TRACEBACK ---
    validation_audit_report = None
    try:
        print("\n--- [STAGE 6] SAVING VALIDATION REPORT ---")
        xc_report = build_validation_report(
            base_dir=base_dir,
            val_dir=val_dir,
            path_final_assay=path_final_assay,
            success=success
        )
        validation_audit_report = xc_report
        if xc_report:
            # Append cross-contamination summary to ai_expert_report.json for AI chatbot context
            try:
                ai_rep_path = base_dir / "ai_expert_report.json"
                ai_rep = {}
                if ai_rep_path.exists():
                    with open(ai_rep_path, "r", encoding="utf-8") as rf:
                        ai_rep = json.load(rf)
                ai_rep["cross_contamination_traceback"] = xc_report
                with open(ai_rep_path, "w", encoding="utf-8") as rf:
                    json.dump(ai_rep, rf, ensure_ascii=False, indent=4)
                print(f"   📊 Cross-contamination traceback appended to ai_expert_report.json")
            except Exception as ae:
                print(f"   ⚠️ Could not append cross-contamination to AI report: {ae}")
    except Exception as stage6_err:
        print(f"   ⚠️ Stage 6 error (non-fatal): {stage6_err}")

    # --- GENERATE AUDIT TRAIL AFTER FINAL REPORTS ARE ATTACHED ---
    try:
        write_audit_trail(base_dir, base_params, success, path_final_assay, validation_audit_report)
    except Exception as e:
        print(f"   ⚠️ Could not generate audit trail: {e}")

    print(f"\n🚀 TOTAL TIME: {format_duration(time.time() - pipeline_start)}")


def build_validation_report(
    base_dir: Path,
    val_dir: Path,
    path_final_assay: Path,
    success: bool
) -> dict | None:
    """
    Save structured validation results into 4_validation_report/ folder and
    run cross-contamination traceback algorithm on background results.

    Folder structure created:
        {base_dir}/4_validation_report/
            validation_summary.json       — high-level run summary
            target_hits.csv               — all per-strain target PCR hits
            background_hits.csv           — all per-strain background PCR hits
            cross_contamination_report.json — ranked list of cross-reactive strains

    Returns:
        dict: cross_contamination_report dict (for AI context injection), or None on error.
    """
    report_dir = base_dir / "4_validation_report"
    report_dir.mkdir(parents=True, exist_ok=True)

    try:
        # --- 1. Copy master CSVs into report folder ---
        master_target = val_dir / "master_results_target.csv"
        master_bg     = val_dir / "master_results_background.csv"
        master_stats  = val_dir / "master_pcr_results_summary.csv"

        for src, dst_name in [
            (master_target, "target_hits.csv"),
            (master_bg,     "background_hits.csv"),
            (master_stats,  "pcr_results_summary.csv"),
        ]:
            if src.exists():
                shutil.copy2(str(src), str(report_dir / dst_name))

        if path_final_assay.exists():
            shutil.copy2(str(path_final_assay), str(report_dir / "FINAL_ASSAY.csv"))

        # --- 2. Run Cross-Contamination Traceback ---
        xc_report = None
        if master_bg.exists():
            print("   🔎 Running cross-contamination traceback...")
            xc_report = traceback_cross_contamination(
                bg_hits_csv=str(master_bg),
                stats_csv=str(master_stats) if master_stats.exists() else None,
                final_assay_csv=str(path_final_assay) if path_final_assay.exists() else None
            )
            xc_path = report_dir / "cross_contamination_report.json"
            with open(xc_path, "w", encoding="utf-8") as f_xc:
                json.dump(xc_report, f_xc, ensure_ascii=False, indent=4)
            print(f"   📊 Cross-contamination report: {xc_path.name}")
            n_xc = xc_report.get("accepted_total_cross_reactive_strains", xc_report.get("total_cross_reactive_strains", 0))
            pool_xc = xc_report.get("total_cross_reactive_strains", 0)
            if n_xc > 0:
                print(f"   ⚠️ {n_xc} cross-reactive background strains affect final accepted assays.")
                if pool_xc != n_xc:
                    print(f"      Candidate pool total before final filtering: {pool_xc} strain(s).")
                top = xc_report.get("accepted_cross_reactive_strains", [])[:3]
                for s in top:
                    primers = s.get("accepted_primer_pairs", s.get("primer_pairs", []))
                    print(f"      • {s['strain']} — hit by {len(primers)} accepted primer pair(s), primers: {', '.join(primers[:3])}")
            else:
                if pool_xc > 0:
                    print(f"   ✅ Final accepted assays show no background cross-reactivity. Candidate pool had {pool_xc} cross-reactive strain(s).")
                else:
                    print("   ✅ No cross-reactive background strains detected.")

        # --- 3. Write validation_summary.json ---
        summary = {
            "generated_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "pipeline_success": success,
            "report_folder": str(report_dir.resolve()),
            "files": {
                "target_hits": str(report_dir / "target_hits.csv"),
                "background_hits": str(report_dir / "background_hits.csv"),
                "pcr_results_summary": str(report_dir / "pcr_results_summary.csv"),
                "final_assay": str(report_dir / "FINAL_ASSAY.csv") if success else None,
                "cross_contamination_report": str(report_dir / "cross_contamination_report.json") if xc_report else None
            },
            "cross_contamination_summary": {
                "total_cross_reactive_strains": xc_report.get("total_cross_reactive_strains", 0) if xc_report else 0,
                "severity": xc_report.get("severity", "NONE") if xc_report else "NONE",
                "accepted_total_cross_reactive_strains": xc_report.get("accepted_total_cross_reactive_strains", 0) if xc_report else 0,
                "accepted_severity": xc_report.get("accepted_severity", "NONE") if xc_report else "NONE",
                "top_3_offenders": xc_report.get("top_cross_reactive_strains", [])[:3] if xc_report else []
            }
        }
        with open(report_dir / "validation_summary.json", "w", encoding="utf-8") as f_sum:
            json.dump(summary, f_sum, ensure_ascii=False, indent=4)
        print(f"   ✅ Validation report saved: {report_dir.resolve()}")
        return xc_report

    except Exception as e:
        print(f"   ❌ build_validation_report error: {e}")
        return None


def traceback_cross_contamination(
    bg_hits_csv: str,
    stats_csv: str | None = None,
    final_assay_csv: str | None = None
) -> dict:
    """
    Traceback algorithm: analyze background PCR results to identify which
    background strains are cross-reactive (amplified by target primers).

    Algorithm:
        1. Load results_background.csv — rows are per-strain PCR results.
        2. Group by Strain: count how many distinct primer pairs amplified each strain.
        3. Rank strains by hit_count descending (most cross-reactive first).
        4. For each cross-reactive strain, list which primer pairs hit it and
           what amplicon sizes were produced.
        5. Compute a severity score based on hit counts vs total primer pairs.
        6. Cross-reference with final accepted assays to flag which accepted
           primers are responsible for the cross-reactivity.

    Returns:
        dict with full traceback report suitable for AI chatbot context.
    """
    result = {
        "algorithm": "CrossContamination Traceback v1.0",
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "total_cross_reactive_strains": 0,
        "severity": "NONE",        # NONE / LOW / MODERATE / HIGH / CRITICAL
        "top_cross_reactive_strains": [],
        "per_primer_cross_reactivity": [],
        "accepted_primer_offenders": [],
        "accepted_cross_reactive_strains": [],
        "accepted_total_cross_reactive_strains": 0,
        "accepted_severity": "NONE",
        "ai_summary": ""
    }

    try:
        df_bg = pd.read_csv(bg_hits_csv)

        # Normalize column names (handle whitespace)
        df_bg.columns = df_bg.columns.str.strip()

        # Identify key columns. Do not treat "PCR Product Sequence" as a status column.
        normalized = {c: c.strip().lower().replace("_", " ") for c in df_bg.columns}
        strain_col = next((c for c in df_bg.columns if "strain" in normalized[c]), None)
        primer_col = next((c for c in df_bg.columns if normalized[c] in ("primer pair", "primer", "set id")), None)
        result_col = next((c for c in df_bg.columns if normalized[c] in (
            "result", "status", "pcr result", "pcr status", "amplification status"
        )), None)
        size_col = next((c for c in df_bg.columns if "amplicon" in normalized[c] and "size" in normalized[c]), None)
        product_col = next((c for c in df_bg.columns if normalized[c] in (
            "pcr product sequence", "product sequence", "amplicon sequence"
        )), None)

        if not strain_col or not primer_col:
            result["ai_summary"] = "Background hit CSV missing Strain or Primer Pair columns — traceback could not run."
            return result

        # Filter rows where PCR result is positive (amplification detected)
        if result_col:
            text = df_bg[result_col].astype("string").str.strip().str.lower()
            df_hits = df_bg[
                df_bg[result_col].notna()
                & ~text.isin(['absent', 'negative', 'nan', '', '0', 'none', 'n/a'])
            ].copy()
        else:
            # If no result column, assume presence of amplicon sequence means hit
            if product_col:
                product_df = (
                    df_bg
                    if product_col == "PCR Product Sequence"
                    else df_bg.rename(columns={product_col: "PCR Product Sequence"})
                )
                df_hits = df_bg[positive_product_mask(product_df)].copy()
            else:
                df_hits = df_bg.copy()  # conservative: treat all rows as potential hits

        if df_hits.empty:
            result["ai_summary"] = "No background cross-reactivity detected. All primers are specific."
            return result

        # --- PER-STRAIN TRACEBACK ---
        strain_groups = df_hits.groupby(strain_col)
        strain_records = []

        for strain_name, grp in strain_groups:
            primer_pairs = grp[primer_col].dropna().astype(str).unique().tolist()
            hit_count = len(primer_pairs)

            amplicon_sizes = []
            if size_col and size_col in grp.columns:
                sizes = grp[size_col].dropna().astype(str).tolist()
                amplicon_sizes = [s for s in sizes if s not in ('nan', '', '0')]
            elif product_col and product_col in grp.columns:
                amplicon_sizes = [str(len(s)) + "bp" for s in grp[product_col].dropna().astype(str) if len(s) > 5]

            strain_records.append({
                "strain": str(strain_name),
                "primer_hit_count": hit_count,
                "primer_pairs": primer_pairs,
                "amplicon_sizes": amplicon_sizes[:10],  # cap at 10
                "cross_reactivity_score": round(hit_count / max(1, df_hits[primer_col].nunique()) * 100, 1)
            })

        # Sort descending by hit count
        strain_records.sort(key=lambda x: x["primer_hit_count"], reverse=True)
        result["top_cross_reactive_strains"] = strain_records[:20]  # top 20
        result["total_cross_reactive_strains"] = len(strain_records)

        # --- PER-PRIMER CROSS-REACTIVITY ---
        primer_groups = df_hits.groupby(primer_col)
        primer_records = []
        for primer_name, grp in primer_groups:
            strains_hit = grp[strain_col].dropna().astype(str).unique().tolist()
            primer_records.append({
                "primer_pair": str(primer_name),
                "cross_reactive_strain_count": len(strains_hit),
                "affected_strains": strains_hit[:10]
            })
        primer_records.sort(key=lambda x: x["cross_reactive_strain_count"], reverse=True)
        result["per_primer_cross_reactivity"] = primer_records[:20]

        # --- SEVERITY SCORING ---
        total_bg_strains = df_bg[strain_col].nunique()
        def score_severity(records: list[dict]) -> tuple[str, float]:
            xc_ratio_local = len(records) / max(1, total_bg_strains)
            max_hits_local = records[0]["primer_hit_count"] if records else 0
            if xc_ratio_local == 0:
                return "NONE", xc_ratio_local
            if xc_ratio_local < 0.05 and max_hits_local <= 1:
                return "LOW", xc_ratio_local
            if xc_ratio_local < 0.15 or max_hits_local <= 2:
                return "MODERATE", xc_ratio_local
            if xc_ratio_local < 0.30 or max_hits_local <= 4:
                return "HIGH", xc_ratio_local
            return "CRITICAL", xc_ratio_local

        severity, xc_ratio = score_severity(strain_records)
        result["severity"] = severity

        # --- CROSS-REFERENCE WITH ACCEPTED FINAL ASSAY ---
        accepted_ids = set()
        if final_assay_csv and Path(final_assay_csv).exists():
            try:
                df_final = pd.read_csv(final_assay_csv)
                df_final.columns = df_final.columns.str.strip()
                accepted_ids = set(df_final.get("Set_ID", pd.Series(dtype=str)).astype(str).tolist())
                offenders = [p for p in primer_records if p["primer_pair"] in accepted_ids]
                result["accepted_primer_offenders"] = offenders

                accepted_strain_records = []
                for record in strain_records:
                    accepted_pairs = [p for p in record["primer_pairs"] if p in accepted_ids]
                    if accepted_pairs:
                        accepted_record = dict(record)
                        accepted_record["accepted_primer_pairs"] = accepted_pairs
                        accepted_record["primer_hit_count"] = len(accepted_pairs)
                        accepted_record["accepted_primer_hit_count"] = len(accepted_pairs)
                        accepted_record["cross_reactivity_score"] = round(
                            len(accepted_pairs) / max(1, len(accepted_ids)) * 100, 1
                        )
                        accepted_strain_records.append(accepted_record)
                accepted_strain_records.sort(key=lambda x: x["accepted_primer_hit_count"], reverse=True)
                accepted_severity, accepted_ratio = score_severity(accepted_strain_records)
                result["accepted_cross_reactive_strains"] = accepted_strain_records[:20]
                result["accepted_total_cross_reactive_strains"] = len(accepted_strain_records)
                result["accepted_severity"] = accepted_severity
            except Exception:
                pass

        # --- BUILD AI-READABLE SUMMARY ---
        top_strains_str = "; ".join(
            f"{s['strain']} (hit by {s['primer_hit_count']} primer(s))"
            for s in strain_records[:5]
        )
        top_primers_str = "; ".join(
            f"{p['primer_pair']} → {p['cross_reactive_strain_count']} strain(s)"
            for p in primer_records[:5]
        )
        offender_str = ""
        if accepted_ids and result["accepted_primer_offenders"]:
            offender_str = " ACCEPTED PRIMERS WITH CROSS-REACTIVITY: " + "; ".join(
                p["primer_pair"] for p in result["accepted_primer_offenders"]
            ) + "."
        elif accepted_ids:
            offender_str = " Accepted final assays show no background cross-reactivity; detected issues belong to rejected/non-final candidates."

        result["ai_summary"] = (
            f"Cross-contamination traceback complete. "
            f"Candidate-pool severity: {severity}. "
            f"Final-assay severity: {result['accepted_severity']}. "
            f"{len(strain_records)} background strain(s) showed cross-reactivity "
            f"out of {total_bg_strains} total background strains "
            f"({round(xc_ratio*100,1)}% cross-reactivity rate). "
            f"Most affected strains: {top_strains_str or 'None'}. "
            f"Most promiscuous primers: {top_primers_str or 'None'}."
            f"{offender_str}"
        )

    except Exception as e:
        result["ai_summary"] = f"Traceback failed: {e}"

    return result


def run_multiplex_analysis(args: argparse.Namespace) -> None:
    """
    Execute the multiplex PCR compatibility analysis and kit assembly.

    Args:
        args (argparse.Namespace): Command line arguments containing folder targets and assay settings.
    """
    print("\n--- [STAGE 7] MULTIPLEX PCR COMPATIBILITY ANALYSIS ---")
    from .multiplex import assemble_optimal_multiplex_kits

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    folders = args.folders
    print(f"📂 Scanning {len(folders)} candidate target folders...")
    for f in folders:
        print(f"   - {f}")

    assay_type = getattr(args, 'assay_type', 'qPCR')
    print(f"🔬 Multiplex reaction type: {assay_type}")

    result = assemble_optimal_multiplex_kits(folders, assay_type=assay_type)
    if "error" in result:
        print(f"❌ Multiplex combination error: {result['error']}")
        return

    print(
        f"✨ Evaluated {result['total_combinations']} Multiplex Kit combination(s); "
        f"usable: {result.get('usable_kit_count', 0)}."
    )
    if not result.get("has_usable_kit", False):
        print(f"🚫 No usable multiplex kit passed hard gates. Overall verdict: {result.get('overall_verdict', 'REJECT')}.")
        for reason in result.get("blocking_reasons", [])[:5]:
            print(f"   - {reason}")

    # Write CSV report
    csv_path = out_dir / "MULTIPLEX_KITS.csv"
    with open(csv_path, "w", newline="", encoding="utf-8") as f_csv:
        writer = csv.writer(f_csv)
        writer.writerow(["Kit_Index", "Assay_Type", "Compatibility_Score", "Verdict", "Mean_Primer_Tm", "Mean_Probe_Tm", "Primer_Tm_Span", "Probe_Tm_Span", "Composition"])

        for idx, kit in enumerate(result["top_kits"]):
            comp_details = []
            for a in kit["assays"]:
                size_str = f",Size:{a.get('Amplicon_Size', 150)}bp"
                comp_details.append(f"{a['Target_Name']}(Set:{a['Set_ID']},Fwd:{a['Fwd_Primer']},Rev:{a['Rev_Primer']}{size_str})")
            comp_str = " | ".join(comp_details)
            writer.writerow([
                idx+1, kit.get("assay_type", assay_type), kit["score"], kit["verdict"], kit["mean_primer_tm"], kit["mean_probe_tm"],
                kit["primer_tm_span"], kit["probe_tm_span"], comp_str
            ])

    print(f"💾 Best combinations table written to: {csv_path.resolve()}")

    # Write detailed JSON report for GUI
    details_path = out_dir / "multiplex_details.json"
    try:
        with open(details_path, "w", encoding="utf-8") as f_det:
            json.dump(result, f_det, ensure_ascii=False, indent=4)
        print(f"💾 Biological combination details written to: {details_path.resolve()}")
    except Exception as e:
        print(f"❌ Error writing multiplex_details.json: {e}")

    # Integrate AI Expert consultation
    if not result.get("has_usable_kit", False):
        report_path = out_dir / "ai_expert_report.json"
        hard_report = {
            "best_assay_name": None,
            "overall_verdict": "REJECT",
            "next_action": "REJECT_AND_STOP",
            "tm_analysis": "No usable multiplex kit passed mandatory hard gates.",
            "specificity_sensitivity_balance": "The current multiplex assembly is not suitable for wet-lab use without redesign.",
            "structural_risks": result.get("blocking_reasons", []),
            "clinical_recommendation": (
                "Do not proceed to wet-lab validation with this multiplex set. "
                "For qPCR, redesign or recover valid probes for every target before repeating multiplex evaluation."
            ),
        }
        with open(report_path, "w", encoding="utf-8") as rf:
            json.dump(hard_report, rf, ensure_ascii=False, indent=4)
        print(f"🤖 Deterministic rejection report written to: {report_path.resolve()}")
        return

    if getattr(args, 'ai_base_url', None):
        print("🤖 Consult AI Expert to evaluate Multiplex Kit...")
        try:
            from .ai_core import LocalBackend, AssayEvaluator
            backend = LocalBackend(base_url=args.ai_base_url, model_name=getattr(args, 'ai_model', 'llama3'))
            evaluator = AssayEvaluator(backend)

            # Generate text summary of Top Multiplex Kits
            summary = f"=== MULTIPLEX COMPATIBILITY REPORT ({assay_type.upper()} MODE) ===\n"
            summary += f"Assay Type: {assay_type}\n"
            summary += f"Target Species: {', '.join(result['target_species'])}\n"
            summary += f"Total Cartesians evaluated: {result['total_combinations']}\n\n"
            summary += f"Usable kit count after hard gates: {result.get('usable_kit_count', 0)}\n"
            summary += f"Overall hard-gate verdict: {result.get('overall_verdict', 'UNKNOWN')}\n\n"

            for idx, kit in enumerate(result["top_kits"]):
                summary += f"Multiplex Kit Combo #{idx+1} (Score: {kit['score']}/100, Verdict: {kit['verdict']}):\n"
                summary += f"   - Mean Primer Tm: {kit['mean_primer_tm']}°C | Tm Span: {kit['primer_tm_span']}°C\n"
                summary += f"   - Mean Probe Tm: {kit['mean_probe_tm']}°C | Probe Tm Span: {kit['probe_tm_span']}°C\n"
                summary += f"   - Composition:\n"
                for a in kit["assays"]:
                    probe_str = f" / {a['Probe_Seq']} (Probe)" if assay_type == "qPCR" else ""
                    summary += f"     * {a['Target_Name']}: {a['Fwd_Primer']} (Fwd) / {a['Rev_Primer']} (Rev){probe_str} (Size: {a.get('Amplicon_Size', 150)}bp, Gene: {a['Target_Gene']})\n"

                if kit.get("size_warnings"):
                    summary += f"   - Size/Amplicon Warnings: 📏 {'; '.join(kit['size_warnings'])}\n"
                if kit["hairpin_warnings"]:
                    summary += f"   - Hairpin Warnings: ⚠️ {'; '.join(kit['hairpin_warnings'])}\n"
                if kit["cross_dimer_warnings"]:
                    summary += f"   - Cross-Dimer Warnings: 🚨 {'; '.join(kit['cross_dimer_warnings'])}\n"
                if kit.get("probe_warnings"):
                    summary += f"   - Probe Warnings: 🚫 {'; '.join(kit['probe_warnings'])}\n"
                if kit.get("blocking_reasons"):
                    summary += f"   - Hard Blocking Reasons: 🚫 {'; '.join(kit['blocking_reasons'])}\n"
                summary += "\n"

            prompt = (
                f"Below is the multiplex PCR primer/probe compatibility report (Reaction type: {assay_type}) on multiple target agents.\n"
                "Perform a comprehensive cross-comparison analysis between these Multiplex Kits, evaluate physicochemical stability (Tm balance, Primer-Dimers, Hairpins, and product size separation on gel if Conventional) to approve the most optimal Multiplex Kit for multi-agent molecular diagnostics.\n\n"
                f"{summary}"
            )

            lang = getattr(args, 'language', 'en')
            if lang == "en":
                instruction = (
                    "You are a Molecular Biology Diagnostics Consultant. Evaluate the provided multiplex qPCR combinations report.\n"
                    f"The assay type is: {assay_type}.\n"
                    + ("For qPCR: TaqMan probes must be separated by different fluorophores (e.g. FAM/HEX), product sizes are not required to be different but should ideally be small (70-200bp) for high efficiency.\n" if assay_type == "qPCR" else
                       "For Conventional PCR: product sizes MUST be clearly separated (by at least 20-50 bp gap between targets) so they can be resolved as individual bands on an agarose gel, and sizes can go up to 1000 bp.\n") +
                    "Recommend the single optimal Multiplex Kit combination by biophysical alignment.\n"
                    "Explain your reasoning in detail and suggest a qPCR/PCR Wet-lab multiplex cycling protocol (dye assignments, primer concentrations, salt mix) to optimize reaction efficiency.\n"
                    "Return strict JSON with all fields written in English:\n"
                    "{\n"
                    "  \"best_assay_name\": \"Name or Index of optimal Multiplex Combo\",\n"
                    "  \"overall_verdict\": \"ACCEPT | MARGINAL | REJECT\",\n"
                    "  \"next_action\": \"ACCEPT_AND_STOP\",\n"
                    "  \"tm_analysis\": \"Review of Tm balances and size separation alignments of the recommended kit\",\n"
                    "  \"specificity_sensitivity_balance\": \"Cross-hybridization and dimer risk overview among the targets\",\n"
                    "  \"structural_risks\": [],\n"
                    "  \"clinical_recommendation\": \"Detailed qPCR/PCR protocol suggestions and clinical evaluation\"\n"
                    "}"
                )
            else:
                instruction = (
                    "Bạn là Chuyên gia Tư vấn Chẩn đoán Sinh học Phân tử. Hãy đánh giá báo cáo tổ hợp multiplex PCR được cung cấp.\n"
                    f"Kiểu phản ứng Multiplex: {assay_type}.\n"
                    + ("Đối với qPCR: bắt buộc dùng probe (TaqMan), các kênh màu fluorophore khác nhau (ví dụ FAM/HEX), kích thước sản phẩm không bắt buộc khác nhau nhưng nên nhỏ (70-200bp) để tối ưu hiệu suất.\n" if assay_type == "qPCR" else
                       "Đối với Conventional PCR: kích thước sản phẩm (PCR product sizes) BẮT BUỘC phải khác nhau và cách biệt từ 20-50 bp để có thể phân tách thành các vạch riêng biệt trên agarose gel điện di, kích thước có thể dài lên tới 1000 bp.\n") +
                    "Đề xuất bộ Multiplex Kit tối ưu nhất dựa trên độ cân bằng Tm, kích thước sản phẩm và giảm thiểu rủi ro tạo dimer.\n"
                    "Giải thích lập luận chi tiết và đề xuất quy trình qPCR/PCR Wet-lab multiplex (nhãn huỳnh quang nếu có, nồng độ mồi, muối Mg2+) để tối ưu hiệu suất phản ứng.\n"
                    "Trả về JSON nghiêm ngặt với các trường viết bằng Tiếng Việt:\n"
                    "{\n"
                    "  \"best_assay_name\": \"Tên hoặc số của Multiplex Combo tối ưu nhất\",\n"
                    "  \"overall_verdict\": \"ACCEPT | MARGINAL | REJECT\",\n"
                    "  \"next_action\": \"ACCEPT_AND_STOP\",\n"
                    "  \"tm_analysis\": \"Phân tích độ lệch Tm và khoảng cách kích thước amplicon của bộ mồi được chọn\",\n"
                    "  \"specificity_sensitivity_balance\": \"Đánh giá nguy cơ dimer chéo và bắt cặp sai giữa các loài đích\",\n"
                    "  \"structural_risks\": [],\n"
                    "  \"clinical_recommendation\": \"Hướng dẫn chi tiết wet-lab multiplex qPCR/PCR và lập luận khoa học\"\n"
                    "}"
                )

            report = evaluator.backend.generate_json(prompt, instruction)

            # Clean up Markdown JSON
            if isinstance(report, str):
                cleaned = report.replace("```json", "").replace("```", "").strip()
                report = json.loads(cleaned)

            report_path = out_dir / "ai_expert_report.json"
            with open(report_path, "w", encoding="utf-8") as rf:
                json.dump(report, rf, ensure_ascii=False, indent=4)

            print(f"🤖 AI Expert consultation report written to: {report_path.resolve()}")
        except Exception as e:
            print(f"❌ AI Expert Multiplex Error: {e}")

def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command")
    cmd_pipe = subparsers.add_parser("pipeline")
    cmd_pipe.add_argument("--out", required=True)
    cmd_pipe.add_argument("--params")
    cmd_pipe.add_argument("--local_target")
    cmd_pipe.add_argument("--local_bg")
    cmd_pipe.add_argument("--email")
    cmd_pipe.add_argument("--target_config")
    cmd_pipe.add_argument("--bg_config")
    cmd_pipe.add_argument("--ai_base_url", help="Local AI Base URL (e.g. http://localhost:11434/v1)")
    cmd_pipe.add_argument("--ai_model", help="Model Name (e.g. llama3)")
    cmd_pipe.add_argument("--language", default="en", help="Output language for AI (vi or en)")
    cmd_pipe.add_argument("--shared_context", default=None,
                          help="Path to shared_primer_context.json for cross-target multiplex sequential gating")

    # Register validate_primers command
    cmd_val = subparsers.add_parser("validate_primers")
    cmd_val.add_argument("-c", "--config", required=True, help="Tệp primers.csv")
    cmd_val.add_argument("-t", "--target", help="Thư mục Target Genomes")
    cmd_val.add_argument("-b", "--background", help="Thư mục Background Genomes")
    cmd_val.add_argument("-o", "--output", default="PCR_Advanced_Report.csv")
    cmd_val.add_argument("-s", "--seq", action="store_true", help="Trích xuất trình tự Amplicon")
    cmd_val.add_argument("-e", "--error", type=int, default=4, help="Mismatch tối đa")
    cmd_val.add_argument("-w", "--workers", type=int, default=8, help="Số nhân CPU")
    cmd_val.add_argument("--max_len", type=int, default=1500, help="Độ dài Amplicon tối đa")

    # Register design_probes command
    cmd_probe = subparsers.add_parser("design_probes")
    cmd_probe.add_argument("-c", "--config", required=True, help="Tệp primers.csv chứa danh sách cặp mồi")
    cmd_probe.add_argument("-t", "--target", required=True, help="Thư mục Target Genomes")
    cmd_probe.add_argument("-b", "--background", help="Thư mục Background Genomes (tùy chọn)")
    cmd_probe.add_argument("-o", "--output", default="designed_probes.csv", help="Tệp CSV đầu ra")
    cmd_probe.add_argument("-e", "--error", type=int, default=4, help="Mismatch tối đa")
    cmd_probe.add_argument("--max_len", type=int, default=1500, help="Độ dài Amplicon tối đa")

    # Register multiplex_analyze command
    cmd_multi = subparsers.add_parser("multiplex_analyze")
    cmd_multi.add_argument("-f", "--folders", nargs="+", required=True, help="Danh sách các thư mục targets kết quả thiết kế")
    cmd_multi.add_argument("-o", "--out", required=True, help="Thư mục đầu ra ghi báo cáo multiplex")
    cmd_multi.add_argument("--ai_base_url")
    cmd_multi.add_argument("--ai_model")
    cmd_multi.add_argument("--language", default="en")
    cmd_multi.add_argument("--assay_type", default="qPCR", choices=["qPCR", "Conventional"], help="Kiểu phản ứng Multiplex: qPCR (Probe-based) hoặc Conventional (Gel-based)")

    # Register terminal interactive mode
    subparsers.add_parser("term", help="Interactive terminal wizard")

    args = parser.parse_args()
    if args.command == "pipeline":
        run_full_pipeline(args)
    elif args.command == "validate_primers":
        import sys
        # Cut out 'validate_primers' command to pass arguments to insilico_pcr_advanced
        sys.argv = [sys.argv[0]] + sys.argv[2:]
        from .insilico_pcr_advanced import main as run_advanced_pcr
        run_advanced_pcr()
    elif args.command == "design_probes":
        from .probe_designer import design_probes_for_primers
        design_probes_for_primers(
            primers_csv=args.config,
            target_dir=args.target,
            bg_dir=args.background,
            output_csv=args.output,
            max_error=args.error,
            max_len=args.max_len
        )
    elif args.command == "multiplex_analyze":
        run_multiplex_analysis(args)
    elif args.command == "term":
        from .term import main as run_terminal
        run_terminal()

if __name__ == "__main__":
    multiprocessing.freeze_support()
    main()
