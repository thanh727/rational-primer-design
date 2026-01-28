import argparse
import sys
import json
import csv
import multiprocessing
import datetime
import time
import random 
from pathlib import Path

from .fetcher import SequenceFetcher
from .constructor import LibraryConstructor
from .designer import PrimerDesigner
from .validator import InSilicoValidator
from .prober import ProbeSelector

class DualLogger(object):
    def __init__(self, filename):
        self.terminal = sys.stdout
        self.log = open(filename, "w", encoding='utf-8')
    def write(self, message):
        self.terminal.write(message); self.log.write(message); self.log.flush()
    def flush(self):
        self.terminal.flush(); self.log.flush()

def load_json(path):
    with open(path) as f: return json.load(f)

def load_config_for_fetcher(json_path):
    data = load_json(json_path)
    parsed = {}
    for k, v in data.items():
        query = v[0]
        size = v[1]
        count = v[2] if len(v) > 2 else 0
        parsed[k] = (query, size, count)
    return parsed

def format_duration(seconds):
    m, s = divmod(int(seconds), 60)
    h, m = divmod(m, 60)
    if h > 0: return f"{h}h {m}m {s}s"
    return f"{m}m {s}s"

def run_full_pipeline(args):
    # --- REPRODUCIBILITY ---
    if hasattr(args, 'seed'):
        random.seed(args.seed)
    else:
        random.seed(42) 
    
    pipeline_start = time.time()
    timing_log = []

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
        "design_target_sampling_size": 0,
        "design_background_sampling_size": 100,
        "validation_target_sampling_size": 0,
        "validation_background_sampling_size": 200,
        "design_max_candidates": 50, # Will be overwritten by config/slider
        "cpu_cores": 0,
        "primer_length": 20, "primer_opt_tm": 60.0, "primer_gc_min": 20.0, "primer_gc_max": 80.0,
        "probe_tm_plus": 10.0, "product_size_min": 100, "product_size_max": 350,
        "design_min_conservation": 0.90, "design_max_bg_prevalence": 0.05,
        "min_sensitivity": 95.0, "max_mismatch": 3, "validation_max_cross_reactivity": 5.0,
        "probe_max_mismatch": 3,
        "enable_blast": True
    }

    user_params = {}
    config_file = None
    if args.params: config_file = Path(args.params)
    elif Path("config/parameters.json").exists(): config_file = Path("config/parameters.json")

    if config_file:
        print(f"\n[CONFIG] Attempting to load: {config_file}")
        try:
            user_params = load_json(config_file)
            print("   ✅ Configuration loaded.")
        except Exception as e:
            print(f"   ⚠️ Error reading config ({e}). Using INTERNAL DEFAULTS.")

    base_params = {**defaults, **user_params}
    if base_params["cpu_cores"] <= 0:
        base_params["cpu_cores"] = max(1, multiprocessing.cpu_count() - 2)

    # --- PATHS ---
    path_design_target = base_dir / "1_workspace" / "design" / "target.fasta"
    path_design_bg = base_dir / "1_workspace" / "design" / "background.fasta"
    path_val_target = base_dir / "1_workspace" / "validate" / "target.fasta"
    path_val_bg = base_dir / "1_workspace" / "validate" / "background.fasta"
    path_candidates_csv = base_dir / "2_candidates" / "candidates.csv"
    path_candidates_fasta = base_dir / "2_candidates" / "candidates.fasta"
    path_val_results = base_dir / "3_validation"
    path_val_target_csv = path_val_results / "results_target.csv"
    path_val_stats_csv = path_val_results / "pcr_results_summary.csv"
    path_final_assay = base_dir / "FINAL_ASSAY.csv"

    # ==========================================
    # STAGE 0 & 1: DATA PREP
    # ==========================================
    path_raw_target = base_dir / "0_raw_data" / "target"
    path_raw_bg = base_dir / "0_raw_data" / "background"
    
    if args.local_target and args.local_bg:
        print("\n--- [STAGE 0] LOCAL DATA MODE ---")
        path_raw_target = Path(args.local_target)
        path_raw_bg = Path(args.local_bg)
    elif args.target_config and args.bg_config and args.email:
        print("\n--- [STAGE 0] DOWNLOADING DATA ---")
        t0 = time.time()
        fetcher = SequenceFetcher(email=args.email)
        fetcher.fetch_and_save_all(load_config_for_fetcher(args.target_config), str(path_raw_target))
        fetcher.fetch_and_save_all(load_config_for_fetcher(args.bg_config), str(path_raw_bg))
        timing_log.append(("Download", time.time() - t0))
    else: return print("❌ Error: Missing inputs.")

    print("\n--- [STAGE 1] BUILDING DATASETS ---")
    t0 = time.time()
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
    timing_log.append(("Dataset Construction", time.time() - t0))

    # ==========================================
    # AUTO-OPTIMIZATION LOOP (DYNAMIC)
    # ==========================================
    # 1. Read User's Slider Value
    user_max_cand = int(base_params.get("design_max_candidates", 2000))
    print(f"\n[INFO] User selected Max Candidates: {user_max_cand}")

    # 2. Build Strategies dynamically based on User Input
    # Cycle 1: Strict Settings, User's Candidate Count
    # Cycle 2+: Relax Conservation, Increase Candidates slightly (if possible)
    
    strategies = [
        {"candidates": user_max_cand, "relax_cons": 0.00, "relax_sens": 0.0},
        {"candidates": int(user_max_cand * 1.5), "relax_cons": 0.05, "relax_sens": 0.0},
        {"candidates": int(user_max_cand * 2.0), "relax_cons": 0.10, "relax_sens": 5.0},
        {"candidates": int(user_max_cand * 2.0), "relax_cons": 0.15, "relax_sens": 5.0},
        {"candidates": int(user_max_cand * 3.0), "relax_cons": 0.20, "relax_sens": 10.0}
    ]
    
    success = False
    
    for i, strat in enumerate(strategies):
        attempt_num = i + 1
        current_params = base_params.copy()
        
        current_params['design_max_candidates'] = strat['candidates']
        current_params['design_min_conservation'] -= strat['relax_cons']
        current_params['min_sensitivity'] -= strat['relax_sens']
        
        # Ensure we don't go below 0 or above reasonable limits
        if current_params['design_min_conservation'] < 0.5: current_params['design_min_conservation'] = 0.5

        print(f"\n" + "━"*60)
        print(f"   🔄 CYCLE {attempt_num}/{len(strategies)}: Testing Top {strat['candidates']} Candidates")
        print(f"   🎯 Specs: Cons >={current_params['design_min_conservation']*100:.0f}% | Sens >={current_params['min_sensitivity']:.0f}%")
        print("━"*60)

        with open(base_dir / f"run_parameters_cycle_{attempt_num}.json", 'w') as f:
            json.dump(current_params, f, indent=4)

        # STAGE 2
        print("\n--- [STAGE 2] DESIGNING PRIMERS ---")
        path_candidates_csv.parent.mkdir(parents=True, exist_ok=True)
        if path_candidates_csv.exists(): path_candidates_csv.unlink()
        
        designer = PrimerDesigner(params=current_params) 
        designer.design(str(path_design_target), str(path_design_bg), str(path_candidates_csv))
        
        if not path_candidates_csv.exists():
            print("   ⚠️ No candidates found. Next strategy...")
            continue

        try:
            with open(path_candidates_csv, 'r') as f_in, open(path_candidates_fasta, 'w') as f_out:
                reader = csv.DictReader(f_in)
                count = 0
                for row in reader:
                    f_out.write(f">{row['Set_ID']}_Fwd\n{row['Forward Primer']}\n")
                    f_out.write(f">{row['Set_ID']}_Rev\n{row['Reverse Primer']}\n")
                    count += 1
            print(f"   ✅ [Bridge] Testing {count} candidates.")
        except: 
            print("   ❌ Error in bridge.")
            continue

        # STAGE 3
        print("\n--- [STAGE 3] IN-SILICO VALIDATION ---")
        validator = InSilicoValidator()
        validator.validate(str(path_val_target), str(path_val_bg), str(path_candidates_fasta), str(path_val_results), current_params)
        
        if not path_val_target_csv.exists(): 
            print("   ⚠️ Validation failed. Expanding search...")
            continue

        # STAGE 4
        print("\n--- [STAGE 4] PROBE DESIGN ---")
        prober = ProbeSelector()
        if path_final_assay.exists(): path_final_assay.unlink()
        
        found_probes = prober.design(str(path_val_target_csv), str(path_final_assay), str(path_val_stats_csv), current_params)

        if found_probes:
            print("\n   ✨ SUCCESS! Valid assays generated.")
            success = True
            break
        else:
            print("   ⚠️ No valid probes found. Expanding search...")
            continue

    if success:
        if current_params.get("enable_blast", True):
            print("\n--- [STAGE 5] BLAST ANNOTATION ---")
            prober.run_blast_annotation(str(path_final_assay))
        print(f"\n✅✅✅ PIPELINE COMPLETE! Output: {path_final_assay}")
    else:
        print("\n❌ FAILED after all attempts.")

    print("\n" + "="*40 + "\n   ⏱️  PERFORMANCE SUMMARY\n" + "="*40)
    for stage, duration in timing_log: print(f"   {stage:<25}: {format_duration(duration)}")
    print("-" * 40 + f"\n   🚀 TOTAL TIME             : {format_duration(time.time() - pipeline_start)}\n" + "="*40)
    print(f"\n📝 Log saved to: {log_path}")

def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command")
    cmd_pipe = subparsers.add_parser("pipeline")
    cmd_pipe.add_argument("--out", required=True)
    cmd_pipe.add_argument("--params")
    cmd_pipe.add_argument("--seed", type=int, default=42)
    cmd_pipe.add_argument("--target_config")
    cmd_pipe.add_argument("--bg_config")
    cmd_pipe.add_argument("--email")
    cmd_pipe.add_argument("--local_target")
    cmd_pipe.add_argument("--local_bg")
    subparsers.add_parser("fetch") 
    args = parser.parse_args()
    if args.command == "pipeline": run_full_pipeline(args)

if __name__ == "__main__":
    multiprocessing.freeze_support()
    main()