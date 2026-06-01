"""Run probe selection from an existing validation output directory.

Example:
    python examples/run_prober_only.py --base-dir results_local --annotate
"""

import argparse
import os

from rational_design.prober import ProbeSelector


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-dir", default="results_local", help="Pipeline output directory.")
    parser.add_argument(
        "--detail-csv",
        help="Target validation detail CSV. Defaults to <base-dir>/3_validation/results_target.csv.",
    )
    parser.add_argument(
        "--summary-csv",
        help="Target validation summary CSV. Defaults to <base-dir>/3_validation/pcr_results_summary.csv.",
    )
    parser.add_argument(
        "--background-csv",
        help="Background validation CSV. Defaults to <base-dir>/3_validation/master_results_background.csv if present.",
    )
    parser.add_argument(
        "--output-csv",
        help="Output assay CSV. Defaults to <base-dir>/FINAL_ASSAY_ANNOTATED.csv.",
    )
    parser.add_argument(
        "--annotate",
        action="store_true",
        help="Run remote NCBI BLASTX annotation after probe selection.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    validation_dir = os.path.join(args.base_dir, "3_validation")
    detail_csv = args.detail_csv or os.path.join(validation_dir, "results_target.csv")
    summary_csv = args.summary_csv or os.path.join(validation_dir, "pcr_results_summary.csv")
    background_csv = args.background_csv or os.path.join(validation_dir, "master_results_background.csv")
    output_csv = args.output_csv or os.path.join(args.base_dir, "FINAL_ASSAY_ANNOTATED.csv")

    params = {}
    if os.path.exists(background_csv):
        params["background_results_csv"] = background_csv

    prober = ProbeSelector()
    success = prober.design(detail_csv, output_csv, summary_csv, params=params)
    if not success:
        raise SystemExit("No assay could be generated from the validation files.")

    if args.annotate:
        prober.run_blast_annotation(output_csv)

    print(f"Wrote {output_csv}")


if __name__ == "__main__":
    main()
