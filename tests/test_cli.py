import csv
import json
import sys
import pandas as pd
import pytest
from pathlib import Path
from rational_design.cli import main as cli_main

def test_cli_validate_primers(tmp_path, monkeypatch):
    # Create target directory and a mock FASTA
    target_dir = tmp_path / "target"
    target_dir.mkdir()
    target_fasta = target_dir / "target.fasta"
    target_fasta.write_text(">target_strain\nGCTGATGCAAAATTACTCGGTTTTGAGTACCCACTTCGGTTTGGTTTTTATTTACAGAGGATGAACAGG\n")

    # Create background directory and a mock FASTA
    bg_dir = tmp_path / "background"
    bg_dir.mkdir()
    bg_fasta = bg_dir / "bg.fasta"
    bg_fasta.write_text(">background_strain\nCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n")

    # Create primers CSV
    primers_csv = tmp_path / "primers.csv"
    with open(primers_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["name", "fwd", "rev"])
        writer.writerow(["Set_1", "GCTGATGCAAAATTACTCGG", "CCTGTTCATCCTCTGTAAAT"])

    output_csv = tmp_path / "PCR_Report.csv"

    # Mock sys.argv for the subcommand validate_primers
    args = [
        "rational-design",
        "validate_primers",
        "-c", str(primers_csv),
        "-t", str(target_dir),
        "-b", str(bg_dir),
        "-o", str(output_csv),
        "-e", "2",
        "-w", "1"
    ]
    monkeypatch.setattr(sys, "argv", args)

    # Call CLI main
    cli_main()

    # Assertions
    assert output_csv.exists()
    df_report = pd.read_csv(output_csv)
    assert len(df_report) > 0
    # Check that validation report directory was created
    report_dir = tmp_path / "4_validation_report"
    assert report_dir.exists()
    assert (report_dir / "validation_summary.json").exists()


def test_cli_design_probes(tmp_path, monkeypatch):
    # Create target directory and a mock FASTA
    target_dir = tmp_path / "target"
    target_dir.mkdir()
    target_fasta = target_dir / "target.fasta"
    target_fasta.write_text(">target_strain\nGCTGATGCAAAATTACTCGGTTTTGAGTACCCACTTCGGTTTGGTTTTTATTTACAGAGGATGAACAGG\n")

    # Create primers CSV
    primers_csv = tmp_path / "primers.csv"
    with open(primers_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["name", "fwd", "rev"])
        writer.writerow(["Set_1", "GCTGATGCAAAATTACTCGG", "CCTGTTCATCCTCTGTAAAT"])

    output_csv = tmp_path / "designed_probes.csv"

    # Mock sys.argv for the subcommand design_probes
    args = [
        "rational-design",
        "design_probes",
        "-c", str(primers_csv),
        "-t", str(target_dir),
        "-o", str(output_csv),
        "-e", "2"
    ]
    monkeypatch.setattr(sys, "argv", args)

    # Call CLI main
    cli_main()

    # Assertions
    assert output_csv.exists()
    df_report = pd.read_csv(output_csv)
    assert len(df_report) == 1
    assert df_report.loc[0, "probe"] == "TGAGTACCCACTTCGGTTTGGTT"


def test_cli_multiplex_analyze(tmp_path, monkeypatch):
    # Create mock target directory results
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
            "Probe_Seq": "ATCGATCGATCGATCG",
            "Probe_Tm": 65.0,
            "Target_Gene": "gene",
            "Amplicon_Size": 120,
        }
    ]
    pd.DataFrame(common_rows).to_csv(target_a / "FINAL_ASSAY.csv", index=False)
    pd.DataFrame([{**common_rows[0], "Set_ID": "Set_2", "Fwd_Primer": "GCTAGCTAGCTAGCTAGC", "Probe_Seq": "TTAGGCTTAGGCTTAG"}]).to_csv(
        target_b / "FINAL_ASSAY.csv",
        index=False,
    )

    out_dir = tmp_path / "multiplex_out"

    # Mock sys.argv for the subcommand multiplex_analyze
    args = [
        "rational-design",
        "multiplex_analyze",
        "-f", str(target_a), str(target_b),
        "-o", str(out_dir),
        "--assay_type", "qPCR"
    ]
    monkeypatch.setattr(sys, "argv", args)

    # Call CLI main
    cli_main()

    # Assertions
    assert (out_dir / "MULTIPLEX_KITS.csv").exists()
    assert (out_dir / "multiplex_details.json").exists()


def test_cli_pipeline(tmp_path, monkeypatch):
    import random
    random.seed(42)
    bases = ['A', 'C', 'G', 'T']
    random_seq = lambda n: "".join(random.choice(bases) for _ in range(n))

    fwd = "GCTGATGCAAAATTACTCGG"
    rev_rc = "ATTTACAGAGGATGAACAGG"
    probe = "TGAGTACCCACTTCGGTTTGGTT"
    
    # Create target directory and a mock FASTA
    target_dir = tmp_path / "target"
    target_dir.mkdir()
    target_fasta = target_dir / "target.fasta"
    # Construct a high-complexity target sequence with fwd, probe, and rev_rc
    target_sequence = fwd + "AA" + probe + "GG" + random_seq(150) + rev_rc
    target_fasta.write_text(f">target_strain\n{target_sequence}\n")

    # Create background directory and a mock FASTA
    bg_dir = tmp_path / "background"
    bg_dir.mkdir()
    bg_fasta = bg_dir / "bg.fasta"
    # Background sequence of same length without target primers/probes
    bg_sequence = random_seq(len(target_sequence))
    bg_fasta.write_text(f">background_strain\n{bg_sequence}\n")

    out_dir = tmp_path / "pipeline_out"

    # Create a small parameters.json file
    params_json = tmp_path / "parameters.json"
    with open(params_json, "w") as f:
        json.dump({
            "design_target_sampling_size": 1,
            "design_background_sampling_size": 1,
            "validation_target_sampling_size": 1,
            "validation_background_sampling_size": 1,
            "design_max_candidates": 5,
            "cpu_cores": 1,
            "enable_blast": False
        }, f)

    # Mock sys.argv for the subcommand pipeline
    args = [
        "rational-design",
        "pipeline",
        "--out", str(out_dir),
        "--local_target", str(target_dir),
        "--local_bg", str(bg_dir),
        "--params", str(params_json),
    ]
    monkeypatch.setattr(sys, "argv", args)

    # Call CLI main
    cli_main()

    # Assertions
    assert (out_dir / "FINAL_ASSAY.csv").exists()

