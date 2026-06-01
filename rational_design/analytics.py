import os
import json
import pandas as pd
from datetime import datetime

class ResultAnalysisEngine:
    """
    Result Analysis Engine: Manages, analyzes, and indexes 
    all historical primer design runs in the runs/ directory.
    """
    def __init__(self, runs_dir="runs"):
        self.runs_dir = os.path.abspath(runs_dir)
        os.makedirs(self.runs_dir, exist_ok=True)

    def scan_historical_runs(self) -> list:
        """
        Scan the entire runs/ directory to find all completed runs.
        Returns a list of runs with metadata and summary results.
        """
        historical_runs = []
        if not os.path.exists(self.runs_dir):
            return historical_runs

        # Level 1 folder: Target species/gene name
        for target_name in os.listdir(self.runs_dir):
            target_path = os.path.join(self.runs_dir, target_name)
            if not os.path.isdir(target_path) or target_name.startswith("."):
                continue

            # Level 2 folder: Run sessions (run_YYYYMMDD_HHMMSS)
            for run_name in os.listdir(target_path):
                run_path = os.path.join(target_path, run_name)
                if not os.path.isdir(run_path) or run_name.startswith("."):
                    continue

                final_csv = os.path.join(run_path, "FINAL_ASSAY.csv")
                multiplex_csv = os.path.join(run_path, "MULTIPLEX_KITS.csv")
                ai_json = os.path.join(run_path, "ai_expert_report.json")

                if os.path.exists(multiplex_csv):
                    run_info = {
                        "target_name": target_name.replace("_", " ").title() + " (Multiplex)",
                        "run_folder_name": run_name,
                        "path": run_path,
                        "timestamp": self._parse_timestamp(run_name),
                        "best_assay": "Multiplex Kit",
                        "sensitivity": "N/A",
                        "specificity": "N/A",
                        "target_gene": "Multiple Targets",
                        "total_candidates": 0,
                        "is_multiplex": True
                    }
                    try:
                        df = pd.read_csv(multiplex_csv)
                        run_info["total_candidates"] = len(df)
                        if not df.empty:
                            top_row = df.iloc[0]
                            run_info["best_assay"] = f"Combo #{top_row.get('Kit_Index', 1)}"
                            run_info["sensitivity"] = f"Score: {top_row.get('Compatibility_Score', 'N/A')}"
                            run_info["specificity"] = top_row.get("Verdict", "N/A")
                            run_info["target_gene"] = top_row.get("Composition", "N/A")[:100] + "..."
                    except Exception:
                        pass

                    if os.path.exists(ai_json):
                        try:
                            with open(ai_json, "r", encoding="utf-8") as f:
                                report = json.load(f)
                                run_info["ai_verdict"] = report.get("overall_verdict", "N/A")
                                run_info["clinical_recommendation"] = report.get("clinical_recommendation", "N/A")
                        except Exception:
                            pass

                    historical_runs.append(run_info)

                elif os.path.exists(final_csv):
                    run_info = {
                        "target_name": target_name.replace("_", " ").title(),
                        "run_folder_name": run_name,
                        "path": run_path,
                        "timestamp": self._parse_timestamp(run_name),
                        "best_assay": "N/A",
                        "sensitivity": "N/A",
                        "specificity": "N/A",
                        "target_gene": "N/A",
                        "total_candidates": 0,
                        "is_multiplex": False
                    }

                    try:
                        df = pd.read_csv(final_csv)
                        run_info["total_candidates"] = len(df)
                        if not df.empty:
                            top_row = df.iloc[0]
                            run_info["best_assay"] = top_row.get("Set_ID", "N/A")
                            run_info["sensitivity"] = top_row.get("Sensitivity", "N/A")
                            run_info["specificity"] = top_row.get("Spec", "N/A")
                            run_info["target_gene"] = top_row.get("Target_Gene", "N/A")
                    except Exception:
                        pass

                    if os.path.exists(ai_json):
                        try:
                            with open(ai_json, "r", encoding="utf-8") as f:
                                report = json.load(f)
                                run_info["ai_verdict"] = report.get("overall_verdict", "N/A")
                                run_info["clinical_recommendation"] = report.get("clinical_recommendation", "N/A")
                        except Exception:
                            pass

                    historical_runs.append(run_info)

        # Sort by latest timestamp first
        historical_runs.sort(key=lambda x: x["timestamp"], reverse=True)
        return historical_runs

    def get_run_details(self, run_path: str) -> dict:
        """Read details of all designed assays for a specific run."""
        details = {"assays": [], "report": {}, "is_multiplex": False, "multiplex_details": {}}
        final_csv = os.path.join(run_path, "FINAL_ASSAY.csv")
        multiplex_csv = os.path.join(run_path, "MULTIPLEX_KITS.csv")
        multiplex_json = os.path.join(run_path, "multiplex_details.json")
        ai_json = os.path.join(run_path, "ai_expert_report.json")

        if os.path.exists(multiplex_csv):
            details["is_multiplex"] = True
            try:
                df = pd.read_csv(multiplex_csv)
                details["assays"] = df.to_dict(orient="records")
            except Exception:
                pass
            if os.path.exists(multiplex_json):
                try:
                    with open(multiplex_json, "r", encoding="utf-8") as f:
                        details["multiplex_details"] = json.load(f)
                except Exception:
                    pass
        elif os.path.exists(final_csv):
            try:
                df = pd.read_csv(final_csv)
                details["assays"] = df.to_dict(orient="records")
            except Exception:
                pass

        if os.path.exists(ai_json):
            try:
                with open(ai_json, "r", encoding="utf-8") as f:
                    details["report"] = json.load(f)
            except Exception:
                pass

        return details

    def generate_ai_summary(self, historical_runs: list) -> str:
        """
        Generate a text summary string of the historical result database 
        to pass into the system prompt of the AI Expert.
        """
        if not historical_runs:
            return "No historical primer design runs have been recorded in the system yet."

        summary = f"--- HISTORICAL DESIGN DATABASE ---\n"
        summary += f"The system currently stores a total of {len(historical_runs)} successful primer design runs. Below is the list of results:\n\n"

        for idx, run in enumerate(historical_runs):
            summary += f"{idx+1}. Target: {run['target_name']} (Folder: {run['run_folder_name']})\n"
            summary += f"   - Run Time: {run['timestamp'].strftime('%Y-%m-%d %H:%M:%S')}\n"
            summary += f"   - Best Assay: {run['best_assay']} | Target BLAST Gene: {run['target_gene']}\n"
            summary += f"   - Inclusivity (Sensitivity): {run['sensitivity']} | Exclusivity (Specificity): {run['specificity']}\n"
            if "ai_verdict" in run:
                summary += f"   - AI Expert Verdict: {run['ai_verdict']}\n"
                summary += f"   - Clinical Recommendation: {run['clinical_recommendation']}\n"
            summary += f"   - Storage Path: {run['path']}\n\n"

        summary += (
            "ANALYSIS RULES:\n"
            "1. If the user asks about previously run species, exactly list the species above.\n"
            "2. If the user requests to compare or contrast assays of different species, use the sensitivity, specificity, BLAST target gene, or primer sequences (if detailed) to perform scientific comparative analysis.\n"
            "3. When the user wants to see details or wet-lab instructions for a specific completed species, use the corresponding information above to answer."
        )
        return summary

    def _parse_timestamp(self, run_name: str) -> datetime:
        """Extract datetime from folder format run_YYYYMMDD_HHMMSS."""
        try:
            parts = run_name.split("_")
            if len(parts) >= 3:
                ts_str = f"{parts[1]}_{parts[2]}"
                return datetime.strptime(ts_str, "%Y%m%d_%H%M%S")
        except Exception:
            pass
        return datetime.now()
