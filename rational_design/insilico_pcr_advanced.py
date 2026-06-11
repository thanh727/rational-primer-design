import os
import csv
import gc
import argparse
import multiprocessing
import pandas as pd
from concurrent.futures import ProcessPoolExecutor
from Bio import SeqIO
from Bio.Seq import Seq
import time

# --- SYSTEM OPTIMIZATION FOR MACOS ---
try:
    multiprocessing.set_start_method('spawn', force=True)
except RuntimeError:
    pass

# Prevent math library conflicts when running multicore
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"

# ================= IUPAC CONFIG =================
IUPAC = {
    "A": {"A"}, "T": {"T"}, "G": {"G"}, "C": {"C"},
    "R": {"A", "G"}, "Y": {"C", "T"}, "S": {"G", "C"},
    "W": {"A", "T"}, "K": {"G", "T"}, "M": {"A", "C"},
    "B": {"C", "G", "T"}, "D": {"A", "G", "T"},
    "H": {"A", "C", "T"}, "V": {"A", "C", "G"},
    "N": {"A", "T", "G", "C"}
}

def is_match(q, t):
    """Check primer match with IUPAC characters."""
    return t in IUPAC.get(q, {q})

def _expand_degenerate_seed(seed_seq):
    seed_seq = str(seed_seq).upper()
    results = [""]
    for char in seed_seq:
        bases = IUPAC.get(char, {char})
        next_results = []
        for r in results:
            for b in bases:
                next_results.append(r + b)
        results = next_results
    return [s.encode() for s in results]

class IndustrialEngine:
    def __init__(self, name, fwd, rev, template_path='', max_error=4, max_p=1500, extract_seq=False, probe=''):
        self.name = name
        self.fwd, self.rev = fwd.strip().upper(), rev.strip().upper()
        self.probe = probe.strip().upper() if probe else ''
        self.f_rc = str(Seq(self.fwd).reverse_complement()).upper()
        self.r_rc = str(Seq(self.rev).reverse_complement()).upper()
        self.max_error, self.max_p, self.extract_seq = max_error, max_p, extract_seq
        self.template_path = template_path
        
        # 8bp seeds at 3' end, expanded if they contain degenerate IUPAC bases
        self.seeds = {
            "f": _expand_degenerate_seed(self.fwd[-8:]), 
            "f_rc": _expand_degenerate_seed(self.f_rc[:8]),
            "r": _expand_degenerate_seed(self.rev[-8:]), 
            "r_rc": _expand_degenerate_seed(self.r_rc[:8])
        }

    def check_probe_match(self, amplicon_seq):
        """Check if the probe matches inside the amplicon sequence.
        Supports IUPAC degenerate bases.
        Allows up to 2 mismatches by default.
        """
        if not self.probe:
            return True
        
        probe_len = len(self.probe)
        amp_len = len(amplicon_seq)
        if amp_len < probe_len:
            return False
            
        probe_rc = str(Seq(self.probe).reverse_complement()).upper()
        
        for i in range(amp_len - probe_len + 1):
            sub = amplicon_seq[i : i + probe_len]
            err_fwd = sum(1 for j in range(probe_len) if not is_match(self.probe[j], sub[j]))
            err_rev = sum(1 for j in range(probe_len) if not is_match(probe_rc[j], sub[j]))
            if err_fwd <= 2 or err_rev <= 2:
                return True
        return False

    def get_mm(self, seq_b, primer_str, seed_pos):
        """Calculate detailed mismatch for a binding site."""
        p_len = len(primer_str)
        start = max(0, seed_pos - (p_len - 8) - 2)
        target = seq_b[start : seed_pos + 10].decode()
        best_err = p_len
        for shift in range(min(5, len(target) - p_len + 1)):
            sub = target[shift : shift + p_len]
            err = sum(1 for j in range(p_len) if not is_match(primer_str[j], sub[j]))
            if err < best_err: best_err = err
        return best_err

    def process_genome(self, args):
        """Process PCR on a specific genome."""
        sid, path, is_target = args
        group_label = "Target" if is_target else "Background"
        all_hits = []
        
        try:
            # 1. Fast scan using Binary Read to find Seed
            with open(path, 'rb') as f:
                raw_content = f.read().upper()

            has_f = any(s in raw_content for s in self.seeds["f"]) or any(s in raw_content for s in self.seeds["f_rc"])
            has_r = any(s in raw_content for s in self.seeds["r"]) or any(s in raw_content for s in self.seeds["r_rc"])

            if not (has_f and has_r) and not is_target:
                return [sid, group_label, self.name, "Absent", "N/A", 0, "N/A", 0, "N/A", "N/A"]

            # 2. Detailed analysis of each Contig
            for rec in SeqIO.parse(path, "fasta"):
                # Iterate both original strand (+) and complementary strand (-)
                for strand_type, seq_obj in [("+", rec.seq), ("-", rec.seq.reverse_complement())]:
                    seq_b = str(seq_obj).upper().encode()
                    
                    # Check forward primer pair (F -> R_rc)
                    has_f_seed = any(s in seq_b for s in self.seeds["f"])
                    has_r_rc_seed = any(s in seq_b for s in self.seeds["r_rc"])
                    if has_f_seed and has_r_rc_seed:
                        f_pos = []
                        for s in self.seeds["f"]:
                            f_pos.extend([i for i in range(len(seq_b)) if seq_b.startswith(s, i)])
                        r_pos = []
                        for s in self.seeds["r_rc"]:
                            r_pos.extend([i for i in range(len(seq_b)) if seq_b.startswith(s, i)])
                        for fi in f_pos:
                            for ri in r_pos:
                                start_pos = fi - (len(self.fwd) - 8)
                                end_pos = ri + len(self.rev)
                                if 0 < (end_pos - start_pos) < self.max_p:
                                    fwd_bind = seq_b[start_pos : fi + 8].decode()
                                    rev_bind = seq_b[ri : end_pos].decode()
                                    mm = sum(1 for j in range(len(self.fwd)) if not is_match(self.fwd[j], fwd_bind[j])) + \
                                         sum(1 for j in range(len(self.rev)) if not is_match(self.r_rc[j], rev_bind[j]))
                                    if mm <= self.max_error:
                                        amp_seq = seq_b[start_pos:end_pos].decode()
                                        if self.check_probe_match(amp_seq):
                                            size = end_pos - start_pos
                                            seq = amp_seq if self.extract_seq else "N/A"
                                            all_hits.append([mm, size, "Standard", seq])

                    # Check reverse primer pair (R -> F_rc)
                    has_r_seed = any(s in seq_b for s in self.seeds["r"])
                    has_f_rc_seed = any(s in seq_b for s in self.seeds["f_rc"])
                    if has_r_seed and has_f_rc_seed:
                        r_pos = []
                        for s in self.seeds["r"]:
                            r_pos.extend([i for i in range(len(seq_b)) if seq_b.startswith(s, i)])
                        f_pos = []
                        for s in self.seeds["f_rc"]:
                            f_pos.extend([i for i in range(len(seq_b)) if seq_b.startswith(s, i)])
                        for ri in r_pos:
                            for fi in f_pos:
                                start_pos = ri - (len(self.rev) - 8)
                                end_pos = fi + len(self.fwd)
                                if 0 < (end_pos - start_pos) < self.max_p:
                                    rev_bind = seq_b[start_pos : ri + 8].decode()
                                    fwd_bind = seq_b[fi : end_pos].decode()
                                    mm = sum(1 for j in range(len(self.rev)) if not is_match(self.rev[j], rev_bind[j])) + \
                                         sum(1 for j in range(len(self.fwd)) if not is_match(self.f_rc[j], fwd_bind[j]))
                                    if mm <= self.max_error:
                                        amp_seq = seq_b[start_pos:end_pos].decode()
                                        if self.check_probe_match(amp_seq):
                                            size = end_pos - start_pos
                                            seq = amp_seq if self.extract_seq else "N/A"
                                            all_hits.append([mm, size, "Reverse-Pair", seq])

            if all_hits:
                b = sorted(all_hits, key=lambda x: x[0])[0]
                return [sid, group_label, self.name, "Positive", b[0], b[1], b[2], len(all_hits), "100%", b[3]]
            
            # Diagnose Divergent if PCR is negative but high identity (For Target)
            if is_target and self.template_path and os.path.exists(self.template_path):
                return [sid, group_label, self.name, "Divergent", "N/A", 0, "N/A", 0, "Check Identity", "N/A"]
                            
        except Exception as e:
            return [sid, group_label, self.name, "Error", "N/A", 0, "N/A", 0, "N/A", str(e)]
        
        return [sid, group_label, self.name, "Absent", "N/A", 0, "N/A", 0, "N/A", "N/A"]

# ================= ORCHESTRATOR (MAIN) =================
def main():
    parser = argparse.ArgumentParser(description="TMRC Industrial PCR Engine v4.5")
    parser.add_argument("-c", "--config", required=True, help="primers.csv file")
    parser.add_argument("-t", "--target", help="Target Genomes directory")
    parser.add_argument("-b", "--background", help="Background Genomes directory")
    parser.add_argument("-o", "--output", default="PCR_Advanced_Report.csv")
    parser.add_argument("-s", "--seq", action="store_true", help="Extract Amplicon sequence")
    parser.add_argument("-e", "--error", type=int, default=4, help="Total maximum mismatch")
    parser.add_argument("-w", "--workers", type=int, default=8, help="Number of CPU cores")
    parser.add_argument("--max_len", type=int, default=1500, help="Maximum Amplicon size")
    args = parser.parse_args()

    try:
        with open(args.config, 'r', encoding='utf-8') as f:
            first_line = f.readline()
        delim = ';' if ';' in first_line else ','
        primers = pd.read_csv(args.config, sep=delim).fillna('').to_dict('records')
    except Exception:
        primers = pd.read_csv(args.config).fillna('').to_dict('records')

    tasks = []
    for folder, is_target in [(args.target, True), (args.background, False)]:
        if folder and os.path.exists(folder):
            tasks.extend([(f, os.path.join(folder, f), is_target) 
                          for f in os.listdir(folder) if f.lower().endswith((".fasta", ".fna", ".fa"))])

    print(f"📊 Total genomes to analyze: {len(tasks)}")
    print(f"🧬 Total Markers: {len(primers)}")

    with open(args.output, "w", newline="") as f_out:
        writer = csv.writer(f_out)
        writer.writerow(["id", "group", "marker", "status", "mismatch", "amplicon_size", "pair_type", "copy_number", "identity", "amplicon_sequence"])
        
        for p in primers:
            engine = IndustrialEngine(
                p['name'], 
                p['fwd'], 
                p['rev'], 
                p.get('template', ''), 
                args.error, 
                args.max_len, 
                args.seq,
                probe=p.get('probe', '')
            )
            
            workers = args.workers if args.workers > 0 else max(1, multiprocessing.cpu_count() - 2)
            with ProcessPoolExecutor(max_workers=workers) as executor:
                for res in executor.map(engine.process_genome, tasks):
                    writer.writerow(res)
            gc.collect()

    print(f"✨ Done! Report at: {args.output}")

    # --- POST-PROCESSING: Save structured validation report with cross-contamination traceback ---
    try:
        save_validation_report(args.output, primers)
    except Exception as e:
        print(f"   ⚠️ Post-processing error (non-fatal): {e}")


def save_validation_report(pcr_report_csv: str, primers: list) -> None:
    """
    Create a structured validation report folder alongside the PCR report CSV.
    
    Creates:
        {parent_dir}/4_validation_report/
            PCR_Advanced_Report.csv          — copy of the full report
            target_results.csv               — filtered target-only results
            background_results.csv           — filtered background-only results  
            per_primer_summary.json          — per-primer sensitivity/specificity stats
            cross_contamination_report.json  — traceback of cross-reactive strains
            validation_summary.json          — overall summary for AI context
    """
    import json
    import shutil
    from pathlib import Path

    report_path = Path(pcr_report_csv)
    if not report_path.exists():
        return

    # Create report folder next to the CSV
    parent_dir = report_path.parent
    report_dir = parent_dir / "4_validation_report"
    report_dir.mkdir(parents=True, exist_ok=True)

    print(f"\n--- [POST] SAVING VALIDATION REPORT ---")

    # Copy full report
    shutil.copy2(str(report_path), str(report_dir / report_path.name))

    # Load and split by group
    df = pd.read_csv(pcr_report_csv)
    df.columns = df.columns.str.strip()

    group_col = next((c for c in df.columns if "group" in c.lower()), None)
    if not group_col:
        print("   ⚠️ No 'group' column found in report — skipping structured save.")
        return

    df_target = df[df[group_col].str.strip().str.lower() == "target"]
    df_bg = df[df[group_col].str.strip().str.lower() == "background"]

    df_target.to_csv(str(report_dir / "target_results.csv"), index=False)
    df_bg.to_csv(str(report_dir / "background_results.csv"), index=False)
    print(f"   📁 Target hits: {len(df_target)} rows, Background hits: {len(df_bg)} rows")

    # --- PER-PRIMER SUMMARY ---
    marker_col = next((c for c in df.columns if "marker" in c.lower() or "name" in c.lower()), None)
    status_col = next((c for c in df.columns if "status" in c.lower()), None)

    per_primer_stats = []
    if marker_col and status_col:
        for primer_name in df[marker_col].unique():
            df_p = df[df[marker_col] == primer_name]
            df_p_target = df_p[df_p[group_col].str.strip().str.lower() == "target"]
            df_p_bg = df_p[df_p[group_col].str.strip().str.lower() == "background"]

            total_target = len(df_p_target)
            hits_target = len(df_p_target[df_p_target[status_col].str.strip().str.lower().isin(["positive"])])
            total_bg = len(df_p_bg)
            hits_bg = len(df_p_bg[df_p_bg[status_col].str.strip().str.lower().isin(["positive"])])

            sensitivity = round(hits_target / max(1, total_target) * 100, 1)
            specificity = round((1 - hits_bg / max(1, total_bg)) * 100, 1) if total_bg > 0 else 100.0

            # Extract primer sequences from config if available
            primer_info = next((p for p in primers if str(p.get("name", "")) == str(primer_name)), {})

            per_primer_stats.append({
                "primer_name": str(primer_name),
                "fwd_sequence": primer_info.get("fwd", "N/A"),
                "rev_sequence": primer_info.get("rev", "N/A"),
                "total_target_strains": total_target,
                "target_hits": hits_target,
                "sensitivity_pct": sensitivity,
                "total_background_strains": total_bg,
                "background_hits": hits_bg,
                "specificity_pct": specificity,
                "verdict": "PASS" if sensitivity >= 90.0 and specificity >= 95.0 else "MARGINAL" if sensitivity >= 70.0 else "FAIL"
            })

    with open(report_dir / "per_primer_summary.json", "w", encoding="utf-8") as f:
        json.dump(per_primer_stats, f, ensure_ascii=False, indent=4)

    # Print summary table
    if per_primer_stats:
        print(f"\n   {'Primer':<20} {'Sens%':>8} {'Spec%':>8} {'Verdict':>10}")
        print(f"   {'─'*20} {'─'*8} {'─'*8} {'─'*10}")
        for ps in per_primer_stats:
            v_icon = "✅" if ps["verdict"] == "PASS" else "⚠️" if ps["verdict"] == "MARGINAL" else "❌"
            print(f"   {ps['primer_name']:<20} {ps['sensitivity_pct']:>7.1f}% {ps['specificity_pct']:>7.1f}% {v_icon} {ps['verdict']:>7}")

    # --- CROSS-CONTAMINATION TRACEBACK ---
    xc_report = _traceback_validation_cross_contamination(df_bg, group_col, marker_col, status_col)

    with open(report_dir / "cross_contamination_report.json", "w", encoding="utf-8") as f:
        json.dump(xc_report, f, ensure_ascii=False, indent=4)

    if xc_report["total_cross_reactive_strains"] > 0:
        print(f"\n   ⚠️ Cross-contamination: {xc_report['total_cross_reactive_strains']} background strain(s) amplified.")
        for s in xc_report.get("top_cross_reactive_strains", [])[:3]:
            print(f"      • {s['strain']} — hit by {s['primer_hit_count']} marker(s): {', '.join(s['primer_pairs'][:3])}")
    else:
        print(f"   ✅ No cross-contamination detected.")

    # --- VALIDATION SUMMARY (AI-readable) ---
    summary = {
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "validation_type": "user_primer_validation",
        "total_primers_tested": len(per_primer_stats),
        "primers_passed": len([p for p in per_primer_stats if p["verdict"] == "PASS"]),
        "primers_marginal": len([p for p in per_primer_stats if p["verdict"] == "MARGINAL"]),
        "primers_failed": len([p for p in per_primer_stats if p["verdict"] == "FAIL"]),
        "per_primer_summary": per_primer_stats,
        "cross_contamination": {
            "severity": xc_report.get("severity", "NONE"),
            "total_cross_reactive_strains": xc_report.get("total_cross_reactive_strains", 0),
            "top_offenders": xc_report.get("top_cross_reactive_strains", [])[:5],
            "ai_summary": xc_report.get("ai_summary", "")
        },
        "report_folder": str(report_dir.resolve()),
        "files": {
            "full_report": str(report_dir / report_path.name),
            "target_results": str(report_dir / "target_results.csv"),
            "background_results": str(report_dir / "background_results.csv"),
            "per_primer_summary": str(report_dir / "per_primer_summary.json"),
            "cross_contamination_report": str(report_dir / "cross_contamination_report.json")
        }
    }

    with open(report_dir / "validation_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=4)

    # Also write to parent dir as validation_report.json for easy AI pickup
    with open(parent_dir / "validation_report.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=4)

    print(f"\n   ✅ Validation report saved: {report_dir.resolve()}")


def _traceback_validation_cross_contamination(
    df_bg: pd.DataFrame,
    group_col: str,
    marker_col: str | None,
    status_col: str | None
) -> dict:
    """
    Traceback cross-contamination from validation background results.
    Groups by strain ID to find which background strains were amplified
    and by which markers.
    """
    result = {
        "algorithm": "Validation CrossContamination Traceback v1.0",
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "total_cross_reactive_strains": 0,
        "severity": "NONE",
        "top_cross_reactive_strains": [],
        "per_marker_cross_reactivity": [],
        "ai_summary": ""
    }

    if df_bg.empty or not marker_col or not status_col:
        result["ai_summary"] = "No background data available for cross-contamination analysis."
        return result

    # Identify strain column (usually 'id')
    strain_col = next((c for c in df_bg.columns if c.lower() in ("id", "strain", "genome")), df_bg.columns[0])

    # Filter positive hits only
    df_hits = df_bg[df_bg[status_col].str.strip().str.lower().isin(["positive"])].copy()

    if df_hits.empty:
        result["ai_summary"] = "No background cross-reactivity detected. All tested primers are specific to target genomes only."
        return result

    # --- PER-STRAIN TRACEBACK ---
    strain_records = []
    for strain_name, grp in df_hits.groupby(strain_col):
        markers = grp[marker_col].dropna().astype(str).unique().tolist()

        # Collect mismatch info if available
        mm_col = next((c for c in grp.columns if "mismatch" in c.lower()), None)
        mismatches = []
        if mm_col:
            mismatches = grp[mm_col].dropna().astype(str).tolist()[:5]

        size_col = next((c for c in grp.columns if "amplicon_size" in c.lower() or "size" in c.lower()), None)
        sizes = []
        if size_col:
            sizes = grp[size_col].dropna().astype(str).tolist()[:5]

        strain_records.append({
            "strain": str(strain_name),
            "primer_hit_count": len(markers),
            "primer_pairs": markers,   # Normalized: same key as cli.py traceback
            "mismatches": mismatches,
            "amplicon_sizes": sizes,
            "cross_reactivity_score": round(len(markers) / max(1, df_hits[marker_col].nunique()) * 100, 1)
        })

    strain_records.sort(key=lambda x: x["primer_hit_count"], reverse=True)
    result["top_cross_reactive_strains"] = strain_records[:20]
    result["total_cross_reactive_strains"] = len(strain_records)

    # --- PER-MARKER CROSS-REACTIVITY ---
    marker_records = []
    for marker_name, grp in df_hits.groupby(marker_col):
        strains_hit = grp[strain_col].dropna().astype(str).unique().tolist()
        marker_records.append({
            "primer_pair": str(marker_name),   # Normalized: same key as cli.py traceback
            "cross_reactive_strain_count": len(strains_hit),
            "affected_strains": strains_hit[:10]
        })
    marker_records.sort(key=lambda x: x["cross_reactive_strain_count"], reverse=True)
    result["per_primer_cross_reactivity"] = marker_records[:20]   # Normalized key

    # --- SEVERITY ---
    total_bg = df_bg[strain_col].nunique()
    xc_ratio = len(strain_records) / max(1, total_bg)
    max_hits = strain_records[0]["primer_hit_count"] if strain_records else 0

    if xc_ratio == 0:
        severity = "NONE"
    elif xc_ratio < 0.05 and max_hits <= 1:
        severity = "LOW"
    elif xc_ratio < 0.15 or max_hits <= 2:
        severity = "MODERATE"
    elif xc_ratio < 0.30 or max_hits <= 4:
        severity = "HIGH"
    else:
        severity = "CRITICAL"
    result["severity"] = severity

    # --- AI SUMMARY ---
    top_str = "; ".join(
        f"{s['strain']} (by {s['primer_hit_count']} marker(s))"
        for s in strain_records[:5]
    )
    marker_str = "; ".join(
        f"{m['primer_pair']} → {m['cross_reactive_strain_count']} strain(s)"
        for m in marker_records[:5]
    )
    result["ai_summary"] = (
        f"Validation cross-contamination traceback complete. "
        f"Severity: {severity}. "
        f"{len(strain_records)} background strain(s) showed cross-reactivity "
        f"out of {total_bg} total background strains "
        f"({round(xc_ratio * 100, 1)}% cross-reactivity rate). "
        f"Most affected strains: {top_str or 'None'}. "
        f"Most promiscuous markers: {marker_str or 'None'}."
    )

    return result


if __name__ == "__main__":
    main()
