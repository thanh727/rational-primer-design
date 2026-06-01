import csv
import random

import pandas as pd
from Bio import SeqIO

from rational_design.cli import traceback_cross_contamination
from rational_design.constructor import LibraryConstructor
from rational_design.designer import is_low_complexity_sequence
from rational_design.multiplex import assemble_optimal_multiplex_kits, evaluate_multiplex_kit
from rational_design.prober import ProbeSelector
from rational_design.validator import find_primer_hits, positive_product_mask


def write_fasta(path, record_id, seq):
    path.write_text(f">{record_id}\n{seq}\n")


def test_pigeonhole_hit_finds_mismatch_inside_old_seed_region():
    primer = "ACGTACGTACGT"
    sequence = "TTTACGTACGTACGATTT"

    hits = find_primer_hits(sequence, primer, max_mm=1)

    assert hits
    assert hits[0]["pos"] == 3
    assert hits[0]["mm"] == 1


def test_positive_product_mask_handles_nullable_nan_values():
    df = pd.DataFrame({"PCR Product Sequence": [pd.NA, None, "", "N/A", "ACTGACTG"]})

    assert positive_product_mask(df).tolist() == [False, False, False, False, True]


def test_traceback_does_not_treat_empty_pcr_product_column_as_positive(tmp_path):
    bg_csv = tmp_path / "background.csv"
    with bg_csv.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["Strain", "Primer Pair", "PCR Product Sequence"])
        writer.writerow(["bg1", "Set_1", ""])
        writer.writerow(["bg2", "Set_1", ""])

    report = traceback_cross_contamination(str(bg_csv))

    assert report["total_cross_reactive_strains"] == 0
    assert report["severity"] == "NONE"


def test_traceback_accepts_amplicon_sequence_alias(tmp_path):
    bg_csv = tmp_path / "background.csv"
    with bg_csv.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["Strain", "Primer Pair", "Amplicon Sequence"])
        writer.writerow(["bg1", "Set_1", ""])
        writer.writerow(["bg2", "Set_1", "ACTGACTGACTG"])

    report = traceback_cross_contamination(str(bg_csv))

    assert report["total_cross_reactive_strains"] == 1
    assert report["top_cross_reactive_strains"][0]["strain"] == "bg2"


def test_traceback_separates_candidate_pool_from_final_assays(tmp_path):
    bg_csv = tmp_path / "background.csv"
    final_csv = tmp_path / "FINAL_ASSAY.csv"
    with bg_csv.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["Strain", "Primer Pair", "PCR Product Sequence"])
        writer.writerow(["bg_pool_only", "Set_1", "ACTGACTGACTG"])
        writer.writerow(["bg_final", "Set_2", "TGCATGCATGCA"])
    pd.DataFrame([{"Set_ID": "Set_2"}]).to_csv(final_csv, index=False)

    report = traceback_cross_contamination(str(bg_csv), final_assay_csv=str(final_csv))

    assert report["total_cross_reactive_strains"] == 2
    assert report["accepted_total_cross_reactive_strains"] == 1
    assert report["accepted_cross_reactive_strains"][0]["strain"] == "bg_final"
    assert report["accepted_primer_offenders"][0]["primer_pair"] == "Set_2"


def test_multiplex_qpcr_missing_probe_is_hard_reject():
    kit = evaluate_multiplex_kit(
        [
            {
                "Set_ID": "Set_1",
                "Target_Name": "Target A",
                "Fwd_Primer": "ACGTACGTACGTACGTAC",
                "Rev_Primer": "TGCATGCATGCATGCATG",
                "Probe_Seq": "nan",
                "Amplicon_Size": 120,
            },
            {
                "Set_ID": "Set_2",
                "Target_Name": "Target B",
                "Fwd_Primer": "GCTAGCTAGCTAGCTAGC",
                "Rev_Primer": "CGATCGATCGATCGATCG",
                "Probe_Seq": "",
                "Amplicon_Size": 130,
            },
        ],
        assay_type="qPCR",
    )

    assert kit["score"] == 0.0
    assert kit["verdict"] == "POOR/UNUSABLE"
    assert kit["has_required_probes"] is False
    assert len(kit["probe_warnings"]) == 2


def test_assemble_multiplex_deduplicates_and_rejects_qpcr_without_probes(tmp_path):
    target_a = tmp_path / "target_a"
    target_b = tmp_path / "target_b"
    target_a.mkdir()
    target_b.mkdir()
    common_rows = [
        {
            "Set_ID": "Set_1",
            "Sensitivity": "100%",
            "Spec": "100%",
            "Fwd_Primer": "ACGTACGTACGTACGTAC",
            "Rev_Primer": "TGCATGCATGCATGCATG",
            "Probe_Seq": "N/A",
            "Probe_Tm": 0.0,
            "Target_Gene": "gene",
            "Amplicon_Size": 120,
        },
        {
            "Set_ID": "Set_1_dup",
            "Sensitivity": "100%",
            "Spec": "100%",
            "Fwd_Primer": "ACGTACGTACGTACGTAC",
            "Rev_Primer": "TGCATGCATGCATGCATG",
            "Probe_Seq": "N/A",
            "Probe_Tm": 0.0,
            "Target_Gene": "gene",
            "Amplicon_Size": 120,
        },
    ]
    pd.DataFrame(common_rows).to_csv(target_a / "FINAL_ASSAY.csv", index=False)
    pd.DataFrame([{**common_rows[0], "Set_ID": "Set_2", "Fwd_Primer": "GCTAGCTAGCTAGCTAGC"}]).to_csv(
        target_b / "FINAL_ASSAY.csv",
        index=False,
    )

    result = assemble_optimal_multiplex_kits([str(target_a), str(target_b)], assay_type="qPCR")

    assert result["total_combinations"] == 1
    assert result["usable_kit_count"] == 0
    assert result["overall_verdict"] == "REJECT"
    assert result["top_kits"][0]["verdict"] == "POOR/UNUSABLE"


def test_constructor_respects_sampling_counts(tmp_path):
    target = tmp_path / "target"
    background = tmp_path / "background"
    target.mkdir()
    background.mkdir()
    for i in range(5):
        write_fasta(target / f"t_{i}.fasta", f"target_{i}", "ACGT" * 20)
        write_fasta(background / f"b_{i}.fasta", f"background_{i}", "TGCA" * 20)

    out = {
        "design_target": str(tmp_path / "design_target.fasta"),
        "design_background": str(tmp_path / "design_background.fasta"),
        "validation_target": str(tmp_path / "validation_target.fasta"),
        "validation_background": str(tmp_path / "validation_background.fasta"),
    }
    random.seed(1)
    LibraryConstructor().construct(str(target), str(background), out, {
        "design_target": 2,
        "design_background": 3,
        "validation_target": 0,
        "validation_background": 0,
    })

    assert len(list(SeqIO.parse(out["design_target"], "fasta"))) == 2
    assert len(list(SeqIO.parse(out["design_background"], "fasta"))) == 3
    assert len(list(SeqIO.parse(out["validation_target"], "fasta"))) == 5
    assert len(list(SeqIO.parse(out["validation_background"], "fasta"))) == 5


def test_short_tandem_repeat_primers_are_low_complexity():
    assert is_low_complexity_sequence("ATGC" * 5)
    assert not is_low_complexity_sequence("ATGCGTACGTTAGCCATGTA")


def test_prober_specificity_ignores_empty_background_products(tmp_path, monkeypatch):
    detail_csv = tmp_path / "results_target.csv"
    summary_csv = tmp_path / "pcr_results_summary.csv"
    background_csv = tmp_path / "master_results_background.csv"
    output_csv = tmp_path / "FINAL_ASSAY.csv"

    amplicon = "ACGT" * 40
    pd.DataFrame(
        [
            {"Strain": "t1", "Primer Pair": "Set_1", "PCR Product Sequence": amplicon},
            {"Strain": "t2", "Primer Pair": "Set_1", "PCR Product Sequence": amplicon},
        ]
    ).to_csv(detail_csv, index=False)
    pd.DataFrame(
        [
            {
                "Primer Pair": "Set_1",
                "Sensitivity_Percent": 100.0,
                "Max_Copy_Number": 1,
                "Forward Primer": "ACGTACGTACGTACGTAC",
                "Reverse Primer": "TGCATGCATGCATGCATG",
            }
        ]
    ).to_csv(summary_csv, index=False)
    pd.DataFrame(
        [
            {"Strain": "bg1", "Primer Pair": "Set_1", "PCR Product Sequence": ""},
            {"Strain": "bg2", "Primer Pair": "Set_1", "PCR Product Sequence": "N/A"},
            {"Strain": "bg3", "Primer Pair": "Set_1", "PCR Product Sequence": amplicon},
        ]
    ).to_csv(background_csv, index=False)

    monkeypatch.setattr(
        ProbeSelector,
        "_find_best_probe",
        lambda self, *args, **kwargs: {"seq": "ACGTACGTACGTACGTACGT", "tm": 68.0, "gc": 50.0},
    )

    assert ProbeSelector().design(
        str(detail_csv),
        str(output_csv),
        str(summary_csv),
        {"background_results_csv": str(background_csv)},
    )

    df = pd.read_csv(output_csv)
    assert df.loc[0, "Spec"] == "66.67%"
