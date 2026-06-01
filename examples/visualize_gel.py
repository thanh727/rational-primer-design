"""Plot a virtual gel-size map from pipeline validation results.

Example:
    python examples/visualize_gel.py --base-dir results_local
"""

import argparse
import os

import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-dir", default="results_local", help="Pipeline output directory.")
    parser.add_argument(
        "--results-csv",
        help="Target validation detail CSV. Defaults to <base-dir>/3_validation/results_target.csv.",
    )
    parser.add_argument(
        "--final-csv",
        help="Final assay CSV. Defaults to <base-dir>/FINAL_ASSAY.csv.",
    )
    parser.add_argument("--top", type=int, default=10, help="Number of assays to plot.")
    parser.add_argument("--output", default="virtual_gel_top10.png", help="PNG output path.")
    parser.add_argument("--show", action="store_true", help="Display the plot interactively.")
    return parser.parse_args()


def _amplicon_size_series(df):
    for column in ("Band_Size_bp", "Amplicon_Size", "amplicon_size", "Product_Size"):
        if column in df.columns:
            return pd.to_numeric(df[column], errors="coerce")
    if "PCR Product Sequence" in df.columns:
        seq = df["PCR Product Sequence"].astype("string").str.strip()
        invalid = seq.str.lower().isin(["", "nan", "none", "n/a", "0", "absent"])
        return seq.mask(invalid).str.len()
    raise ValueError("Cannot infer product size from validation CSV.")


def main():
    args = parse_args()
    validation_dir = os.path.join(args.base_dir, "3_validation")
    results_csv = args.results_csv or os.path.join(validation_dir, "results_target.csv")
    final_csv = args.final_csv or os.path.join(args.base_dir, "FINAL_ASSAY.csv")

    df_results = pd.read_csv(results_csv)
    df_final = pd.read_csv(final_csv)
    id_column = "Set_ID" if "Set_ID" in df_final.columns else "Primer Pair"
    top_pairs = df_final.head(args.top)[id_column].astype(str).tolist()

    df_plot = df_results[df_results["Primer Pair"].astype(str).isin(top_pairs)].copy()
    df_plot["Amplicon_Size_bp"] = _amplicon_size_series(df_plot)
    df_plot = df_plot.dropna(subset=["Amplicon_Size_bp"])
    if df_plot.empty:
        raise SystemExit("No positive amplicons found for the selected assays.")

    import matplotlib.pyplot as plt
    import seaborn as sns

    plt.figure(figsize=(12, 8))
    sns.set_style("darkgrid")
    sns.stripplot(
        data=df_plot,
        x="Amplicon_Size_bp",
        y="Primer Pair",
        hue="Primer Pair",
        palette="viridis",
        alpha=0.45,
        jitter=0.2,
        size=4,
    )

    for i, pair in enumerate(top_pairs):
        median_size = df_plot[df_plot["Primer Pair"].astype(str) == pair]["Amplicon_Size_bp"].median()
        if pd.notna(median_size):
            plt.vlines(
                x=median_size,
                ymin=i - 0.2,
                ymax=i + 0.2,
                colors="red",
                linestyles="solid",
                lw=2,
                label="Median" if i == 0 else "",
            )

    plt.title(f"Virtual gel map for top {len(top_pairs)} assays", fontsize=15, fontweight="bold")
    plt.xlabel("PCR product size (bp)", fontsize=12)
    plt.ylabel("Primer pair", fontsize=12)
    plt.legend(title="Assay ID", bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.tight_layout()
    plt.savefig(args.output, dpi=300)
    if args.show:
        plt.show()
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
