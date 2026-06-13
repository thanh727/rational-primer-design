from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path
from typing import Any
import atexit
import psutil


from rational_design.utils import generate_batch_analytical_summary

from rich.console import Console
from rich.panel import Panel
from rich.prompt import Prompt, IntPrompt, Confirm
from rich.table import Table
from rich.progress import Progress, SpinnerColumn, TextColumn, BarColumn, TimeElapsedColumn
from rich.text import Text
from rich import box

console = Console()

FASTA_EXTS = {".fasta", ".fa", ".fna", ".fas"}
RUNS_DIR = Path(os.getcwd()) / "runs"
CONFIG_DIR = Path(__file__).resolve().parents[1] / "config"

ACTIVE_PROCESSES: list[subprocess.Popen] = []

def cleanup_active_processes():
    for p in ACTIVE_PROCESSES:
        if p.poll() is None:
            try:
                parent = psutil.Process(p.pid)
                for child in parent.children(recursive=True):
                    child.kill()
                p.kill()
            except Exception:
                pass

atexit.register(cleanup_active_processes)



def main() -> None:
    try:
        _check_rich()
        while True:
            choice = _main_menu()
            if choice == "exit":
                _print_exit()
                break
            _route(choice)
    except KeyboardInterrupt:
        console.print("\n[yellow]⚠\uFE0F  Interrupted. Goodbye![/]")
        sys.exit(1)


def _check_rich() -> None:
    try:
        import rich
    except ImportError:
        console.print("[red]\u274C 'rich' package is required. Run: pip install rich[/]")
        sys.exit(1)


def _print_header() -> None:
    lang = _get_language()
    title = Text("\U0001F9EC  Rational Primer Design  \U0001F9EC", style="bold cyan")
    subtitle = Text(f"   Interactive Terminal Mode  [dim](\U0001F310 {lang.upper()})[/]", style="dim white")
    console.print(Panel(f"{title}\n{subtitle}", box=box.ROUNDED, border_style="cyan", padding=(1, 2)))
    console.print()


def _print_exit() -> None:
    console.print("\n[green]\U0001F44B Goodbye![/]")


# ─────────────────────────────────────────────
#  LANGUAGE
# ─────────────────────────────────────────────
def _get_language() -> str:
    path = CONFIG_DIR / "language.json"
    try:
        with open(path, "r", encoding="utf-8") as f:
            return json.load(f).get("language", "en")
    except Exception:
        return "en"


def _save_language(lang: str) -> None:
    CONFIG_DIR.mkdir(parents=True, exist_ok=True)
    with open(CONFIG_DIR / "language.json", "w", encoding="utf-8") as f:
        json.dump({"language": lang}, f)


def _configure_language() -> str:
    current = _get_language()
    choice = Prompt.ask("Language / Ng\u00F4n ng\u1EEF", choices=["vi", "en"], default=current)
    if choice != current:
        _save_language(choice)
    return choice


# ─────────────────────────────────────────────
#  MAIN MENU
# ─────────────────────────────────────────────
def _main_menu() -> str:
    _print_header()
    console.print("[bold]Select mode:[/]\n")
    console.print("  [1]  \U0001F9EC  [bold]NCBI Auto Mode[/]           \u2014 Download genomes & design primers")
    console.print("  [2]  \U0001F4C2  [bold]Local Mode[/]                \u2014 Use your own FASTA files")
    console.print("  [3]  \U0001F52C  [bold]Validate Primers[/]          \u2014 Test existing primers in-silico")
    console.print("  [4]  \U0001F9EA  [bold]Multiplex Analysis[/]        \u2014 Combine designs into multiplex kits")
    console.print("  [5]  \U0001F9EC  [bold]Multi-Target Design (Local)  \u2014 Design primers for N targets from local files[/]")
    console.print("  [6]  \U0001F9EC  [bold]Multi-Target Design (NCBI)   \u2014 Design primers for N targets from NCBI[/]")
    console.print("  [7]  \U0001F916  [bold]AI Expert[/]                 \u2014 Evaluate existing results with AI")
    console.print("  [8]  \U0001F4AC  [bold]AI Chat[/]                   \u2014 Interactive AI Expert with auto-run")
    console.print("  [9]  \U0001F4CA  [bold]Run History[/]               \u2014 View past runs & results")
    console.print("  [10] \u2753  [bold]Language[/]                   \u2014 Switch language (vi/en)")
    console.print("  [11] \u2699\uFE0F  [bold]AI Setup[/]                 \u2014 Detect & configure local AI models")
    console.print("  [12] \u274C  [bold]Exit[/]\n")

    choice = Prompt.ask("Enter your choice", choices=[str(i) for i in range(1, 13)], default="1")
    mapping = {
        "1": "ncbi", "2": "local", "3": "validate", "4": "multiplex",
        "5": "multiplex_local", "6": "multiplex_ncbi",
        "7": "ai_expert", "8": "ai_chat", "9": "history", "10": "language",
        "11": "ai_setup", "12": "exit",
    }
    return mapping[choice]


def _route(choice: str) -> None:
    if choice == "ncbi":
        _wizard_ncbi()
    elif choice == "local":
        _wizard_local()
    elif choice == "validate":
        _wizard_validate()
    elif choice == "multiplex":
        _wizard_multiplex()
    elif choice == "multiplex_local":
        _wizard_multiplex_local()
    elif choice == "multiplex_ncbi":
        _wizard_multiplex_ncbi()
    elif choice == "ai_expert":
        _wizard_ai_expert()
    elif choice == "ai_chat":
        _wizard_ai_chat()
    elif choice == "history":
        _show_history()
    elif choice == "language":
        _configure_language()
        console.print(f"[green]\u2705 Language set to {_get_language().upper()}[/]")
    elif choice == "ai_setup":
        _wizard_ai_setup()


# ─────────────────────────────────────────────
#  LANGUAGE HELPER FOR PIPELINE COMMANDS
# ─────────────────────────────────────────────
def _lang() -> str:
    return _get_language()


# ─────────────────────────────────────────────
#  NCBI AUTO MODE
# ─────────────────────────────────────────────
def _wizard_ncbi() -> None:
    _print_header()
    console.print(Panel("[bold]\U0001F9EC  NCBI Auto Mode[/]\nDownload genomes from NCBI and design primers.", border_style="cyan"))

    email = Prompt.ask("NCBI email", default=_load_saved_email())
    _save_email(email)

    targets = _collect_ncbi_queries("Target", "t")
    if not targets:
        console.print("[yellow]\u26A0\uFE0F  At least one target species required.[/]")
        _wait_and_return()
        return

    backgrounds = _collect_ncbi_queries("Background", "b", show_sample_size=True)

    params_path = _configure_params()
    project_name = Prompt.ask("Project name", default=f"ncbi_{datetime.now().strftime('%Y%m%d_%H%M%S')}")
    out_dir = RUNS_DIR / project_name
    out_dir.mkdir(parents=True, exist_ok=True)

    ai_base_url, ai_model = _configure_ai()

    t_conf_path = CONFIG_DIR / "_term_t_conf.json"
    b_conf_path = CONFIG_DIR / "_term_b_conf.json"
    _save_fetcher_config(t_conf_path, targets)
    _save_fetcher_config(b_conf_path, backgrounds)

    cmd = [
        sys.executable, "-m", "rational_design.cli", "pipeline",
        "--out", str(out_dir),
        "--email", email,
        "--target_config", str(t_conf_path),
        "--bg_config", str(b_conf_path),
        "--params", str(params_path),
        "--language", _lang(),
    ]
    if ai_base_url:
        cmd.extend(["--ai_base_url", ai_base_url])
    if ai_model:
        cmd.extend(["--ai_model", ai_model])

    _run_subprocess(cmd, "\U0001F9EC NCBI Auto Pipeline")
    _show_run_summary(out_dir)
    _wait_and_return()


def _collect_ncbi_queries(label: str, prefix: str, show_sample_size: bool = False) -> dict[str, list]:
    items: dict[str, list] = {}
    console.print(f"\n[bold]{label} species[/] (leave query empty to finish)")
    idx = 1
    while True:
        query = Prompt.ask(f"  {label} #{idx} query", default="")
        if not query:
            if not items and idx == 1:
                console.print(f"[yellow]\u26A0\uFE0F  Need at least one {label}.[/]")
                continue
            break
        size = 0.0
        count = 500
        if show_sample_size:
            size_str = Prompt.ask(f"  Min genome size (Mb)", default="0.0")
            try:
                size = float(size_str)
            except ValueError:
                size = 0.0
        count_str = Prompt.ask(f"  Max genomes to download", default="500")
        try:
            count = int(count_str)
        except ValueError:
            count = 500
        items[f"{prefix}{idx}"] = [query, size, count, "genome"]
        idx += 1
    return items


# ─────────────────────────────────────────────
#  LOCAL MODE
# ─────────────────────────────────────────────
def _wizard_local() -> None:
    _print_header()
    console.print(Panel("[bold]\U0001F4C2  Local Mode[/]\nUse your own FASTA files for primer design.", border_style="green"))

    target_path = _prompt_path("Target FASTA file or directory", must_exist=True)
    bg_path = _prompt_path("Background FASTA file or directory", must_exist=True, allow_empty=False)

    params_path = _configure_params()
    project_name = Prompt.ask("Project name", default=f"local_{datetime.now().strftime('%Y%m%d_%H%M%S')}")
    out_dir = RUNS_DIR / project_name
    out_dir.mkdir(parents=True, exist_ok=True)

    ai_base_url, ai_model = _configure_ai()

    cmd = [
        sys.executable, "-m", "rational_design.cli", "pipeline",
        "--out", str(out_dir),
        "--local_target", str(target_path),
        "--local_bg", str(bg_path),
        "--params", str(params_path),
        "--language", _lang(),
    ]
    if ai_base_url:
        cmd.extend(["--ai_base_url", ai_base_url])
    if ai_model:
        cmd.extend(["--ai_model", ai_model])

    _run_subprocess(cmd, "\U0001F4C2 Local Pipeline")
    _show_run_summary(out_dir)
    _wait_and_return()


# ─────────────────────────────────────────────
#  VALIDATE PRIMERS
# ─────────────────────────────────────────────
def _wizard_validate() -> None:
    _print_header()
    console.print(Panel("[bold]\U0001F52C  Validate Primers[/]\nTest existing primers against target & background genomes.", border_style="yellow"))

    csv_path = _prompt_path("Primers CSV file", must_exist=True, extensions={".csv"})
    target_path = _prompt_path("Target genomes directory", must_exist=True)
    bg_path = _prompt_path("Background genomes directory", must_exist=True, allow_empty=False)
    output_path = _prompt_path("Output CSV", default="PCR_Advanced_Report.csv")
    max_mm_str = Prompt.ask("Max mismatches", default="4")
    workers_str = Prompt.ask("CPU workers", default="8")

    cmd = [
        sys.executable, "-m", "rational_design.cli", "validate_primers",
        "-c", str(csv_path),
        "-t", str(target_path),
        "-b", str(bg_path),
        "-o", str(output_path),
        "-e", max_mm_str,
        "-w", workers_str,
        "--max_len", "1500",
    ]

    _run_subprocess(cmd, "\U0001F52C Primer Validation")
    _wait_and_return()


# ─────────────────────────────────────────────
#  MULTIPLEX ANALYSIS (re-analyze existing folders)
# ─────────────────────────────────────────────
def _wizard_multiplex() -> None:
    _print_header()
    console.print(Panel("[bold]\U0001F9EA  Multiplex Analysis[/]\nCombine multiple design results into multiplex kits.", border_style="magenta"))

    folders: list[str] = []
    console.print("\n[bold]Design result folders[/] (each folder should contain FINAL_ASSAY.csv)")
    while True:
        p = _prompt_path(f"  Folder #{len(folders) + 1}", must_exist=True, allow_empty=True)
        if not p:
            if not folders:
                console.print("[yellow]\u26A0\uFE0F  At least one folder required.[/]")
                continue
            break
        folders.append(str(p))

    project_name = Prompt.ask("Project name", default=f"multiplex_{datetime.now().strftime('%Y%m%d_%H%M%S')}")
    out_dir = RUNS_DIR / project_name
    out_dir.mkdir(parents=True, exist_ok=True)

    assay_type = Prompt.ask("Assay type", choices=["qPCR", "Conventional"], default="qPCR")
    ai_base_url, ai_model = _configure_ai()

    cmd = [
        sys.executable, "-m", "rational_design.cli", "multiplex_analyze",
        "--folders", *folders,
        "--out", str(out_dir),
        "--assay_type", assay_type,
        "--language", _lang(),
    ]
    if ai_base_url:
        cmd.extend(["--ai_base_url", ai_base_url])
    if ai_model:
        cmd.extend(["--ai_model", ai_model])

    _run_subprocess(cmd, "\U0001F9EA Multiplex Analysis")
    _show_run_summary(out_dir)
    _wait_and_return()


# ─────────────────────────────────────────────
#  MULTI-TARGET DESIGN (LOCAL)
# ─────────────────────────────────────────────
def _wizard_multiplex_local() -> None:
    _print_header()
    console.print(Panel("[bold]\U0001F9EC  Multi-Target Design (Local)[/]\nDesign primers for multiple targets sequentially using local FASTA files.", border_style="cyan"))

    targets: list[str] = []
    console.print("\n[bold]Target FASTA files/directories[/] (at least 2, leave empty to finish)")
    while True:
        p = _prompt_path(f"  Target #{len(targets) + 1}", must_exist=True, allow_empty=True)
        if not p:
            if len(targets) < 2:
                console.print("[yellow]\u26A0\uFE0F  At least 2 targets required for multiplex.[/]")
                continue
            break
        targets.append(p)

    bg_path = _prompt_path("Shared Background FASTA file or directory", must_exist=True, allow_empty=False)
    params_path = _configure_params()
    project_name = Prompt.ask("Project name", default=f"multiplex_local_{datetime.now().strftime('%Y%m%d_%H%M%S')}")
    out_dir = RUNS_DIR / project_name
    out_dir.mkdir(parents=True, exist_ok=True)

    shared_context_path = out_dir / "shared_primer_context.json"

    _run_sequential_multiplex(targets, bg_path, params_path, out_dir, shared_context_path, is_ncbi=False)


def _run_sequential_multiplex(
    targets: list[str],
    bg_path: str,
    params_path: Path,
    out_dir: Path,
    shared_context_path: Path,
    is_ncbi: bool,
    email: str = "",
) -> None:
    lang = _lang()
    ai_base_url, ai_model = _configure_ai()
    target_folders: list[Path] = []

    for idx, target_source in enumerate(targets, 1):
        target_folder = out_dir / f"t{idx}"
        target_folder.mkdir(parents=True, exist_ok=True)
        target_folders.append(target_folder)

        _write_json(target_folder / "metadata.json", {"target_name": Path(target_source).stem})

        virtual_bg = out_dir / "0_raw_data" / f"bg_for_t{idx}"
        virtual_bg.mkdir(parents=True, exist_ok=True)

        if Path(bg_path).exists():
            _copy_sequence_files(Path(bg_path), virtual_bg)

        for other_idx, other_source in enumerate(targets, 1):
            if other_idx == idx:
                continue
            _copy_sequence_files(Path(other_source), virtual_bg, prefix=f"target_bg_t{other_idx}_")

        console.print(f"\n[bold cyan]\u23F3  Target #{idx}: {Path(target_source).stem}[/]")

        if is_ncbi and email:
            fetcher_cmd = [
                sys.executable, "-m", "rational_design.cli", "pipeline",
                "--out", str(target_folder),
                "--email", email,
                "--target_config", str(target_source),
                "--bg_config", str(_write_temp_fetcher_config(virtual_bg)),
                "--params", str(params_path),
                "--language", lang,
            ]
            if ai_base_url:
                fetcher_cmd.extend(["--ai_base_url", ai_base_url])
            if ai_model:
                fetcher_cmd.extend(["--ai_model", ai_model])
            _run_subprocess(fetcher_cmd, f"\U0001F9EC NCBI Target #{idx}")
        else:
            cmd = [
                sys.executable, "-m", "rational_design.cli", "pipeline",
                "--out", str(target_folder),
                "--local_target", str(target_source),
                "--local_bg", str(virtual_bg),
                "--params", str(params_path),
                "--language", lang,
                "--shared_context", str(shared_context_path),
            ]
            if ai_base_url:
                cmd.extend(["--ai_base_url", ai_base_url])
            if ai_model:
                cmd.extend(["--ai_model", ai_model])
            _run_subprocess(cmd, f"\U0001F4C2 Local Target #{idx}")

        _accumulate_shared_context(target_folder, shared_context_path)

    console.print("\n[bold]\U0001F9EA Running multiplex analysis across all targets...[/]")
    assay_type = Prompt.ask("Assay type", choices=["qPCR", "Conventional"], default="qPCR")
    mux_cmd = [
        sys.executable, "-m", "rational_design.cli", "multiplex_analyze",
        "--folders", *[str(f) for f in target_folders],
        "--out", str(out_dir),
        "--assay_type", assay_type,
        "--language", lang,
    ]
    if ai_base_url:
        mux_cmd.extend(["--ai_base_url", ai_base_url])
    if ai_model:
        mux_cmd.extend(["--ai_model", ai_model])
    _run_subprocess(mux_cmd, "\U0001F9EA Multiplex Analysis")

    _show_run_summary(out_dir)
    _wait_and_return()


def _accumulate_shared_context(target_folder: Path, shared_context_path: Path) -> None:
    final_csv = target_folder / "FINAL_ASSAY.csv"
    ai_report = target_folder / "ai_expert_report.json"
    if not final_csv.exists():
        return
    try:
        import pandas as pd
        df = pd.read_csv(final_csv)
        if df.empty:
            return
        row = df.iloc[0]
        target_name = "Unknown"
        meta_path = target_folder / "metadata.json"
        if meta_path.exists():
            try:
                with open(meta_path) as f:
                    meta = json.load(f)
                    target_name = meta.get("target_name", target_folder.name)
            except Exception:
                pass

        entry = {
            "target_name": target_name,
            "Fwd_Primer": str(row.get("Fwd_Primer", "")),
            "Rev_Primer": str(row.get("Rev_Primer", "")),
            "Probe_Seq": str(row.get("Probe_Seq", "")),
            "Fwd_Tm": float(row.get("Fwd_Tm", 0) or 0),
            "Rev_Tm": float(row.get("Rev_Tm", 0) or 0),
        }
        context = []
        if shared_context_path.exists():
            try:
                with open(shared_context_path, "r", encoding="utf-8") as f:
                    context = json.load(f)
            except Exception:
                context = []
        context.append(entry)
        with open(shared_context_path, "w", encoding="utf-8") as f:
            json.dump(context, f, ensure_ascii=False, indent=2)
        console.print(f"   [dim]\u2705 Shared context updated: {len(context)} target(s)[/]")
    except Exception:
        pass


def _copy_sequence_files(source: Path, destination: Path, prefix: str = "") -> None:
    if source.is_file() and source.suffix.lower() in FASTA_EXTS:
        shutil.copy2(source, destination / f"{prefix}{source.name}")
        return
    if not source.is_dir():
        return
    for path in source.iterdir():
        if path.is_file() and path.suffix.lower() in FASTA_EXTS:
            shutil.copy2(path, destination / f"{prefix}{path.name}")


def _write_json(path: Path, data: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, ensure_ascii=False, indent=2)


def _write_temp_fetcher_config(source_dir: Path) -> str:
    config = {}
    for idx, f in enumerate(sorted(source_dir.iterdir()), 1):
        if f.suffix.lower() in FASTA_EXTS:
            config[f"b{idx}"] = [f.stem, 0.0, 0, "genome"]
    path = source_dir / "_fetcher_config.json"
    _save_fetcher_config(path, config)
    return str(path)


# ─────────────────────────────────────────────
#  MULTI-TARGET DESIGN (NCBI)
# ─────────────────────────────────────────────
def _wizard_multiplex_ncbi() -> None:
    _print_header()
    console.print(Panel("[bold]\U0001F9EC  Multi-Target Design (NCBI)[/]\nDesign primers for multiple targets from NCBI, each using others as background.", border_style="cyan"))

    email = Prompt.ask("NCBI email", default=_load_saved_email())
    _save_email(email)

    targets = _collect_ncbi_queries("Target", "t")
    if not targets or len(targets) < 2:
        console.print("[yellow]\u26A0\uFE0F  At least 2 targets required.[/]")
        _wait_and_return()
        return

    backgrounds = _collect_ncbi_queries("Background", "b", show_sample_size=True)

    params_path = _configure_params()
    project_name = Prompt.ask("Project name", default=f"multiplex_ncbi_{datetime.now().strftime('%Y%m%d_%H%M%S')}")
    out_dir = RUNS_DIR / project_name
    out_dir.mkdir(parents=True, exist_ok=True)

    shared_context_path = out_dir / "shared_primer_context.json"

    _run_sequential_multiplex_ncbi(targets, backgrounds, email, params_path, out_dir, shared_context_path)


def _run_sequential_multiplex_ncbi(
    targets: dict[str, list],
    backgrounds: dict[str, list],
    email: str,
    params_path: Path,
    out_dir: Path,
    shared_context_path: Path,
) -> None:
    lang = _lang()
    ai_base_url, ai_model = _configure_ai()

    raw_dir = out_dir / "0_raw_data"
    shared_bg = raw_dir / "shared_bg"
    shared_bg.mkdir(parents=True, exist_ok=True)

    if backgrounds:
        console.print(f"\n[cyan]Downloading shared background from NCBI...[/]")
        b_conf_path = CONFIG_DIR / "_term_mux_b_conf.json"
        _save_fetcher_config(b_conf_path, backgrounds)
        fetcher = _run_subprocess([
            sys.executable, "-m", "rational_design.cli", "pipeline",
            "--out", str(out_dir / "_fetch_bg_only"),
            "--email", email,
            "--target_config", str(b_conf_path),
            "--bg_config", str(b_conf_path),
            "--params", str(params_path),
        ], "\U0001F9EC Downloading shared background")

    target_folders: list[Path] = []
    target_dirs: list[Path] = []

    for idx, (key, val) in enumerate(targets.items(), 1):
        t_conf_path = CONFIG_DIR / f"_term_mux_t{idx}_conf.json"
        _save_fetcher_config(t_conf_path, {key: val})

        t_dir = raw_dir / f"t{idx}"
        t_dir.mkdir(parents=True, exist_ok=True)

        fetch_cmd = [
            sys.executable, "-m", "rational_design.cli", "pipeline",
            "--out", str(out_dir / f"_fetch_t{idx}"),
            "--email", email,
            "--target_config", str(t_conf_path),
            "--bg_config", str(t_conf_path),
            "--params", str(params_path),
        ]
        _run_subprocess(fetch_cmd, f"\U0001F9EC Downloading Target #{idx}")

        if t_dir.exists():
            target_dirs.append(t_dir)

    for idx in range(1, len(targets) + 1):
        target_folder = out_dir / f"t{idx}"
        target_folder.mkdir(parents=True, exist_ok=True)
        target_folders.append(target_folder)

        target_name = list(targets.keys())[idx - 1]
        query = targets[target_name][0]
        _write_json(target_folder / "metadata.json", {"target_name": query})

        virtual_bg = raw_dir / f"bg_for_t{idx}"
        virtual_bg.mkdir(parents=True, exist_ok=True)

        if shared_bg.exists():
            _copy_sequence_files(shared_bg, virtual_bg)
        for jdx, td in enumerate(target_dirs):
            if jdx + 1 == idx:
                continue
            _copy_sequence_files(td, virtual_bg, prefix=f"target_bg_t{jdx+1}_")

        t_src = target_dirs[idx - 1] if idx - 1 < len(target_dirs) else raw_dir / f"t{idx}"
        console.print(f"\n[bold cyan]\u23F3  Target #{idx}: {query}[/]")
        cmd = [
            sys.executable, "-m", "rational_design.cli", "pipeline",
            "--out", str(target_folder),
            "--local_target", str(t_src),
            "--local_bg", str(virtual_bg),
            "--params", str(params_path),
            "--language", lang,
            "--shared_context", str(shared_context_path),
        ]
        if ai_base_url:
            cmd.extend(["--ai_base_url", ai_base_url])
        if ai_model:
            cmd.extend(["--ai_model", ai_model])
        _run_subprocess(cmd, f"\U0001F9EC NCBI Target #{idx}")

        _accumulate_shared_context(target_folder, shared_context_path)

    console.print("\n[bold]\U0001F9EA Running multiplex analysis across all targets...[/]")
    assay_type = Prompt.ask("Assay type", choices=["qPCR", "Conventional"], default="qPCR")
    mux_cmd = [
        sys.executable, "-m", "rational_design.cli", "multiplex_analyze",
        "--folders", *[str(f) for f in target_folders],
        "--out", str(out_dir),
        "--assay_type", assay_type,
        "--language", lang,
    ]
    if ai_base_url:
        mux_cmd.extend(["--ai_base_url", ai_base_url])
    if ai_model:
        mux_cmd.extend(["--ai_model", ai_model])
    _run_subprocess(mux_cmd, "\U0001F9EA Multiplex Analysis")

    _show_run_summary(out_dir)
    _wait_and_return()


# ─────────────────────────────────────────────
#  RUN HISTORY
# ─────────────────────────────────────────────
def _show_history() -> None:
    _print_header()
    console.print(Panel("[bold]\U0001F4CA  Run History[/]", border_style="blue"))

    if not RUNS_DIR.exists():
        console.print("[yellow]No runs found yet.[/]")
        _wait_and_return()
        return

    runs = sorted(RUNS_DIR.iterdir(), key=lambda p: p.stat().st_mtime, reverse=True)
    if not runs:
        console.print("[yellow]No runs found in 'runs/' directory.[/]")
        _wait_and_return()
        return

    table = Table(box=box.ROUNDED, border_style="blue")
    table.add_column("#", style="dim")
    table.add_column("Project", style="cyan")
    table.add_column("Latest Run", style="green")
    table.add_column("Has Results", style="yellow")

    for idx, proj_dir in enumerate(runs[:20], 1):
        if not proj_dir.is_dir():
            continue
        run_folders = sorted(proj_dir.iterdir(), key=lambda p: p.stat().st_mtime, reverse=True)
        latest = ""
        has_results = "\u274C"
        if run_folders:
            latest_run = run_folders[0]
            latest = datetime.fromtimestamp(latest_run.stat().st_mtime).strftime("%Y-%m-%d %H:%M")
            final_csv = latest_run / "FINAL_ASSAY.csv"
            audit = latest_run / "audit_trail.json"
            if final_csv.exists() or audit.exists():
                has_results = "\u2705"
        table.add_row(str(idx), proj_dir.name, latest, has_results)

    console.print(table)

    detail = Confirm.ask("View details of a run?", default=False)
    if detail:
        idx_str = Prompt.ask("Enter run number", default="1")
        try:
            idx = int(idx_str) - 1
            if 0 <= idx < len(runs):
                _show_run_detail(runs[idx])
        except ValueError:
            pass

    _wait_and_return()


def _show_run_detail(proj_dir: Path) -> None:
    run_folders = sorted(proj_dir.iterdir(), key=lambda p: p.stat().st_mtime, reverse=True)
    if not run_folders:
        return

    latest = run_folders[0]
    console.print(f"\n[bold cyan]\U0001F4C2 {proj_dir.name} / {latest.name}[/]")
    console.print(f"   Path: [dim]{latest}[/]")

    _show_run_summary(latest)

    report_path = latest / "ai_expert_report.json"
    if report_path.exists():
        try:
            with open(report_path, "r", encoding="utf-8") as f:
                report = json.load(f)
            console.print(f"\n[bold]\U0001F916 AI Expert Report[/]")
            console.print(f"   Best assay: [green]{report.get('best_assay_name', 'N/A')}[/]")
            console.print(f"   Verdict: [bold]{report.get('overall_verdict', 'N/A')}[/]")
            rec = report.get('clinical_recommendation', '')
            if rec:
                console.print(f"   Recommendation: {rec[:200]}")
        except Exception:
            pass


def _show_run_summary(out_dir: Path) -> None:
    final_csv = out_dir / "FINAL_ASSAY.csv"
    audit_path = out_dir / "audit_trail.json"

    if final_csv.exists():
        try:
            import pandas as pd
            df = pd.read_csv(final_csv)
            console.print(f"\n[bold green]\u2705 FINAL_ASSAY.csv generated[/] \u2014 {len(df)} assay(s)")
            preview = df.head(5)
            table = Table(box=box.SIMPLE, border_style="green")
            for col in preview.columns[:6]:
                table.add_column(str(col), style="cyan")
            for _, row in preview.iterrows():
                vals = [str(row[col])[:30] for col in preview.columns[:6]]
                table.add_row(*vals)
            console.print(table)
            if len(df) > 5:
                console.print(f"   ... and {len(df) - 5} more rows")
        except Exception:
            pass

    if audit_path.exists():
        console.print(f"   \U0001F4C4 audit_trail.json available")


# ─────────────────────────────────────────────
#  SHARED HELPERS
# ─────────────────────────────────────────────
def _prompt_path(
    label: str,
    must_exist: bool = False,
    allow_empty: bool = False,
    default: str = "",
    extensions: set[str] | None = None,
) -> str:
    while True:
        val = Prompt.ask(f"{label}", default=default)
        if not val and allow_empty:
            return ""
        if not val:
            continue
        p = Path(val).expanduser().resolve()
        if must_exist and not p.exists():
            console.print(f"[red]\u274C Path does not exist: {p}[/]")
            continue
        if extensions and p.is_file() and p.suffix.lower() not in extensions:
            console.print(f"[red]\u274C Expected extension: {', '.join(extensions)}[/]")
            continue
        return str(p)


def _save_email(email: str) -> None:
    CONFIG_DIR.mkdir(parents=True, exist_ok=True)
    path = CONFIG_DIR / "user_settings.json"
    try:
        with open(path, "w") as f:
            json.dump({"saved_email": email}, f)
    except Exception:
        pass


def _load_saved_email() -> str:
    path = CONFIG_DIR / "user_settings.json"
    try:
        with open(path, "r") as f:
            data = json.load(f)
            return data.get("saved_email", "")
    except Exception:
        return ""


def _save_fetcher_config(path: Path, data: dict[str, list]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, ensure_ascii=False, indent=2)
    return path


def _configure_params() -> Path:
    path = CONFIG_DIR / "parameters.json"
    if path.exists():
        if not Confirm.ask("Use existing parameters.json?", default=True):
            return _create_params_interactive()
        return path
    return _create_params_interactive()


def _create_params_interactive() -> Path:
    CONFIG_DIR.mkdir(parents=True, exist_ok=True)
    params: dict[str, Any] = {}

    console.print("\n[bold]Pipeline parameters[/] (press Enter to accept default)\n")
    params["design_target_sampling_size"] = int(Prompt.ask("  Design target sample size", default="0"))
    params["design_background_sampling_size"] = int(Prompt.ask("  Design background sample size", default="50"))
    params["validation_target_sampling_size"] = int(Prompt.ask("  Validation target sample size", default="0"))
    params["validation_background_sampling_size"] = int(Prompt.ask("  Validation background sample size", default="0"))
    params["design_max_candidates"] = int(Prompt.ask("  Max candidates", default="50"))
    params["cpu_cores"] = int(Prompt.ask("  CPU cores (0=auto)", default="0"))
    params["design_min_conservation"] = float(Prompt.ask("  Min conservation", default="0.75"))
    params["auto_relax_constraints"] = Confirm.ask("  Auto-relax constraints?", default=True)
    params["min_sensitivity"] = float(Prompt.ask("  Min sensitivity (%)", default="90.0"))
    params["product_size_min"] = int(Prompt.ask("  Min product size (bp)", default="120"))
    params["product_size_max"] = int(Prompt.ask("  Max product size (bp)", default="400"))
    params["max_mismatch"] = int(Prompt.ask("  Max mismatches", default="2"))
    params["validation_max_cross_reactivity"] = float(Prompt.ask("  Max cross-reactivity (%)", default="5.0"))
    params["enable_blast"] = Confirm.ask("  Enable BLAST annotation?", default=True)
    params["degenerate_primers"] = Confirm.ask("  Enable degenerate primers?", default=True)
    if params["degenerate_primers"]:
        params["max_iupac_per_primer"] = int(Prompt.ask("  Max IUPAC positions per primer", default="2"))

    path = CONFIG_DIR / "parameters.json"
    with open(path, "w", encoding="utf-8") as f:
        json.dump(params, f, indent=2, ensure_ascii=False)
    console.print(f"[dim]   Saved to {path}[/]")
    return path


def _load_ai_config() -> dict[str, str]:
    path = CONFIG_DIR / "ai_settings.json"
    try:
        with open(path, "r", encoding="utf-8") as f:
            return json.load(f)
    except Exception:
        return {}


def _save_ai_config(base_url: str, model: str) -> None:
    CONFIG_DIR.mkdir(parents=True, exist_ok=True)
    path = CONFIG_DIR / "ai_settings.json"
    with open(path, "w", encoding="utf-8") as f:
        json.dump({"base_url": base_url, "model": model}, f, indent=2)


def _configure_ai() -> tuple[str | None, str | None]:
    saved = _load_ai_config()
    saved_url = saved.get("base_url", "http://localhost:11434/v1")
    saved_model = saved.get("model", "llama3")

    if saved and Confirm.ask(f"Use saved AI config ({saved_url}, {saved_model})?", default=True):
        return saved_url, saved_model

    if not Confirm.ask("Enable AI Expert evaluation?", default=False):
        return None, None

    base_url = Prompt.ask("  AI API base URL", default=saved_url)
    model = Prompt.ask("  AI model name", default=saved_model)
    _save_ai_config(base_url, model)
    return base_url, model


# ─────────────────────────────────────────────
#  AI SETUP WIZARD
# ─────────────────────────────────────────────
def _wizard_ai_setup() -> None:
    """Detect available local AI servers and let user select one."""
    import urllib.request
    import urllib.error

    _print_header()
    console.print(Panel("[bold]\u2699\uFE0F  AI Setup[/]\nDetect local AI servers and configure model.", border_style="cyan"))

    known_endpoints = [
        ("http://localhost:11434/v1", "Ollama"),
        ("http://localhost:1234/v1",  "LM Studio"),
        ("http://localhost:8000/v1",  "LocalAI"),
        ("http://localhost:8080/v1",  "GPT4All"),
    ]

    with Progress(
        SpinnerColumn(), TextColumn("[progress.description]{task.description}"), transient=True,
    ) as progress:
        progress.add_task("[cyan]Scanning for local AI servers...", total=None)
        servers = []
        for base_url, provider in known_endpoints:
            try:
                from openai import OpenAI
                client = OpenAI(api_key="sk-local", base_url=base_url, timeout=3.0)
                models = sorted({m.id for m in client.models.list().data})
                if models:
                    servers.append((base_url, provider, models))
                    console.print(f"   [green]\u2713[/] {provider} at {base_url} [dim]({len(models)} models)[/]")
            except Exception:
                try:
                    root = base_url[:-3] if base_url.endswith("/v1") else base_url
                    req = urllib.request.Request(f"{root}/api/tags", headers={"Accept": "application/json"})
                    with urllib.request.urlopen(req, timeout=3.0) as resp:
                        payload = json.loads(resp.read().decode("utf-8"))
                    models = sorted({str(item.get("name") or item.get("model")) for item in payload.get("models", []) if item.get("name") or item.get("model")})
                    if models:
                        servers.append((base_url, f"{provider} (raw API)", models))
                        console.print(f"   [green]\u2713[/] {provider} (raw) at {base_url} [dim]({len(models)} models)[/]")
                except Exception:
                    console.print(f"   [dim]\u2014[/] {provider} at {base_url} [dim]no response[/]")
                    continue

    if not servers:
        console.print("\n[red]\u2717 No local AI servers detected.[/]")
        console.print("   [yellow]Start Ollama (https://ollama.com) or LM Studio first.[/]")
        console.print("   [yellow]Run this setup again after starting your AI server.[/]")
        return

    if len(servers) == 1:
        base_url, provider, models = servers[0]
        console.print(f"\n[green]\u2713 Found 1 server:[/] {provider} [dim]({base_url})[/]")
    else:
        console.print(f"\n[green]\u2713 Found {len(servers)} servers:[/]")
        choices = [f"{prov} \u2014 {url}" for url, prov, _ in servers]
        sel = Prompt.ask("Select server", choices=[str(i + 1) for i in range(len(servers))], default="1")
        base_url, provider, models = servers[int(sel) - 1]

    console.print(f"\nModels available on [bold]{provider}[/]:")
    for i, m in enumerate(models, 1):
        console.print(f"   [{i}] {m}")

    model_choices = [str(i + 1) for i in range(len(models))]
    sel = Prompt.ask("Select model number", choices=model_choices, default="1")
    selected_model = models[int(sel) - 1]

    _save_ai_config(base_url, selected_model)
    console.print(f"\n[green]\u2705 AI settings saved![/]")
    console.print(f"   Base URL: [cyan]{base_url}[/]")
    console.print(f"   Model:    [cyan]{selected_model}[/]")
    console.print(f"\n[dim]   File: {CONFIG_DIR / 'ai_settings.json'}[/]")


# ─────────────────────────────────────────────
#  AI EXPERT MODE (evaluate existing results)
# ─────────────────────────────────────────────
def _wizard_ai_expert() -> None:
    _print_header()
    console.print(Panel("[bold]\U0001F916  AI Expert[/]\nLoad existing results & consult the AI Expert.", border_style="magenta"))

    ai_base_url, ai_model = _configure_ai()
    if not ai_base_url:
        _wait_and_return()
        return

    run_dir = _select_run_or_path()
    if not run_dir:
        _wait_and_return()
        return

    _display_run_results(run_dir)
    _run_ai_evaluation(run_dir, ai_base_url, ai_model)
    _wait_and_return()


# ─────────────────────────────────────────────
#  AI CHAT MODE (interactive with auto-run)
# ─────────────────────────────────────────────
def _wizard_ai_chat() -> None:
    _print_header()
    console.print(Panel("[bold]\U0001F4AC  AI Chat[/]\nInteractive chat with AI Expert. AI can propose and auto-run pipeline jobs.", border_style="magenta"))

    ai_base_url, ai_model = _configure_ai()
    if not ai_base_url:
        _wait_and_return()
        return

    lang = _lang()
    auto_run = Confirm.ask("Allow AI to auto-run pipeline proposals?", default=True)
    email = _load_saved_email()
    if not email and auto_run:
        email = Prompt.ask("NCBI email (required for NCBI proposals)", default="")

    system_prompt = (
        "B\u1EA1n l\u00E0 chuy\u00EAn gia AI h\u1ED7 tr\u1EE3 thi\u1EBFt k\u1EBF m\u1ED3i PCR. "
        "H\u00E3y tr\u1EA3 l\u1EDDi b\u1EB1ng ti\u1EBFng Vi\u1EC7t." if lang == "vi"
        else "You are an AI expert assistant for PCR primer design. Respond in English."
    )

    messages = [{"role": "system", "content": system_prompt}]
    console.print(f"\n[dim]AI: {ai_base_url} / {ai_model} | Auto-run: {auto_run} | Language: {lang.upper()}[/]")
    console.print("[dim]Type 'exit' to return to menu, 'clear' to reset conversation.[/]\n")

    try:
        from .ai_core import LocalBackend, AIExpertAgent
        backend = LocalBackend(base_url=ai_base_url, model_name=ai_model)
        agent = AIExpertAgent(backend)
    except Exception as e:
        console.print(f"[red]\u274C AI module unavailable: {e}[/]")
        _wait_and_return()
        return

    with Progress(SpinnerColumn(), TextColumn("[progress.description]{task.description}"), console=console) as progress:
        task = progress.add_task("[cyan]Connecting to AI...", total=None)
        try:
            greeting = agent.chat(messages, expert_report=None, language=lang)
        except Exception:
            greeting = "Hello! How can I help you with primer design today?"
        progress.update(task, visible=False)

    console.print(f"[bold cyan]AI:[/] {greeting}\n")

    while True:
        user_input = Prompt.ask("[bold green]You[/]")
        if user_input.lower() in ("exit", "quit"):
            break
        if user_input.lower() == "clear":
            messages = [{"role": "system", "content": system_prompt}]
            console.print("[dim]Conversation cleared.[/]")
            continue

        messages.append({"role": "user", "content": user_input})

        with Progress(SpinnerColumn(), TextColumn("[progress.description]{task.description}"), console=console) as progress:
            task = progress.add_task("[cyan]AI thinking...", total=None)
            try:
                raw_reply = agent.chat(messages, expert_report=None, language=lang)
            except Exception as e:
                console.print(f"[red]\u274C AI Error: {e}[/]")
                continue
            progress.update(task, visible=False)

        proposal = _extract_proposal_from_reply(raw_reply)
        display_reply = _strip_proposal_from_reply(raw_reply) if proposal else raw_reply
        messages.append({"role": "assistant", "content": display_reply})

        console.print(f"[bold cyan]AI:[/] {display_reply}\n")

        if proposal and auto_run and proposal.get("run_immediately", False):
            console.print(f"\n[bold]\U0001F916 AI proposed an action: [cyan]{proposal.get('action')}[/][/]")
            
            if proposal.get("action") in {"propose_design", "propose_multiplex"}:
                t_list = _items_from_proposal(proposal.get("targets") or proposal.get("target") or [], default_count=500)
                b_list = _items_from_proposal(proposal.get("backgrounds") or proposal.get("background") or [], default_count=50)
                console.print("\n[bold cyan]📋 Proposed Targets:[/]")
                for idx, t in enumerate(t_list, 1):
                    console.print(f"  [{idx}] {t['query']} (Max download: {t['count']}, Min size: {t['size']} Mb)")
                if b_list:
                    console.print("[bold yellow]📋 Proposed Backgrounds (Exclusion):[/]")
                    for idx, b in enumerate(b_list, 1):
                        console.print(f"  [{idx}] {b['query']} (Max download: {b['count']}, Min size: {b['size']} Mb)")
                else:
                    console.print("[yellow]📋 No background exclusion species proposed.[/]")
            
            choice = Prompt.ask("\n  Choose action: [r]un proposal, [e]dit configuration, [c]ancel", choices=["r", "e", "c"], default="r")
            if choice == "r":
                _execute_ai_proposal(proposal, email, lang, ai_base_url, ai_model)
            elif choice == "e":
                proposal = _edit_ai_proposal(proposal)
                _execute_ai_proposal(proposal, email, lang, ai_base_url, ai_model)
            else:
                console.print("[yellow]⚠️ AI proposal cancelled.[/]")


def _extract_proposal_from_reply(text: str) -> dict[str, Any] | None:
    candidates = [
        block.strip()
        for block in re.findall(r"```json\s*([\s\S]*?)\s*```", text)
        if '"action"' in block
    ]
    decoder = json.JSONDecoder()
    for match in re.finditer(r"\{", text):
        try:
            parsed, _ = decoder.raw_decode(text[match.start():])
        except json.JSONDecodeError:
            continue
        if isinstance(parsed, dict) and "action" in parsed:
            candidates.append(json.dumps(parsed))
    for candidate in reversed(candidates):
        try:
            parsed = json.loads(candidate)
        except json.JSONDecodeError:
            continue
        if parsed.get("action") in {"propose_design", "propose_local_design", "propose_validation", "propose_multiplex"}:
            return parsed
    return None


def _strip_proposal_from_reply(text: str) -> str:
    def replace_block(match: re.Match) -> str:
        return "" if '"action"' in match.group(1) else match.group(0)
    return re.sub(r"```json\s*([\s\S]*?)\s*```", replace_block, text).strip()


def _edit_ai_proposal(proposal: dict[str, Any]) -> dict[str, Any]:
    targets = _items_from_proposal(proposal.get("targets") or proposal.get("target") or [], default_count=500)
    backgrounds = _items_from_proposal(proposal.get("backgrounds") or proposal.get("background") or [], default_count=50)

    # Edit targets
    if Confirm.ask("\nModify proposed Targets?", default=False):
        new_targets = []
        for idx, t in enumerate(targets, 1):
            console.print(f"\nTarget #{idx}: [bold]{t['query']}[/]")
            action = Prompt.ask("Action: [k]eep, [m]odify, [d]elete", choices=["k", "m", "d"], default="k")
            if action == "k":
                new_targets.append(t)
            elif action == "m":
                query = Prompt.ask("  Query", default=t["query"])
                count = IntPrompt.ask("  Max genomes to download", default=t["count"])
                new_targets.append({"query": query, "size": t["size"], "count": count, "type": t["type"]})
        
        while True:
            add_more = Confirm.ask("Add another Target species?", default=False)
            if not add_more:
                break
            query = Prompt.ask("  Query", default="")
            if query:
                count = IntPrompt.ask("  Max genomes to download", default=500)
                new_targets.append({"query": query, "size": 0.0, "count": count, "type": "genome"})
        targets = new_targets

    # Edit backgrounds
    if Confirm.ask("\nModify proposed Backgrounds (Exclusion)?", default=False):
        new_backgrounds = []
        for idx, b in enumerate(backgrounds, 1):
            console.print(f"\nBackground #{idx}: [bold]{b['query']}[/]")
            action = Prompt.ask("Action: [k]eep, [m]odify, [d]elete", choices=["k", "m", "d"], default="k")
            if action == "k":
                new_backgrounds.append(b)
            elif action == "m":
                query = Prompt.ask("  Query", default=b["query"])
                count = IntPrompt.ask("  Max genomes to download", default=b["count"])
                new_backgrounds.append({"query": query, "size": b["size"], "count": count, "type": b["type"]})
        
        while True:
            add_more = Confirm.ask("Add another Background species?", default=False)
            if not add_more:
                break
            query = Prompt.ask("  Query", default="")
            if query:
                count = IntPrompt.ask("  Max genomes to download", default=50)
                new_backgrounds.append({"query": query, "size": 0.0, "count": count, "type": "genome"})
        backgrounds = new_backgrounds

    proposal["targets"] = targets
    proposal["backgrounds"] = backgrounds
    return proposal


def _execute_ai_proposal(proposal: dict[str, Any], email: str, lang: str, ai_base_url: str | None, ai_model: str | None) -> None:
    action = proposal.get("action")
    params_path = _configure_params()

    if action == "propose_local_design":
        target = proposal.get("local_target") or proposal.get("target", "")
        background = proposal.get("local_background") or proposal.get("background", "")
        if not target or not background:
            console.print("[red]\u274C Proposal missing local_target or local_background[/]")
            return
        project_name = Prompt.ask("Project name", default=f"ai_local_{datetime.now().strftime('%Y%m%d_%H%M%S')}")
        out_dir = RUNS_DIR / project_name
        out_dir.mkdir(parents=True, exist_ok=True)
        cmd = [
            sys.executable, "-m", "rational_design.cli", "pipeline",
            "--out", str(out_dir),
            "--local_target", str(Path(target).expanduser().resolve()),
            "--local_bg", str(Path(background).expanduser().resolve()),
            "--params", str(params_path),
            "--language", lang,
        ]
        if ai_base_url:
            cmd.extend(["--ai_base_url", ai_base_url])
        if ai_model:
            cmd.extend(["--ai_model", ai_model])
        _run_subprocess(cmd, "\U0001F4C2 AI-Proposed Local Pipeline")
        _show_run_summary(out_dir)

    elif action == "propose_design":
        if not email:
            console.print("[red]\u274C NCBI email required for this proposal.[/]")
            return
        targets = _items_from_proposal(proposal.get("targets") or proposal.get("target") or [], default_count=500)
        backgrounds = _items_from_proposal(proposal.get("backgrounds") or proposal.get("background") or [], default_count=50)
        if not targets:
            console.print("[red]\u274C Proposal missing targets.[/]")
            return
        t_conf = {f"t{i+1}": [t["query"], t["size"], t["count"], t["type"]] for i, t in enumerate(targets)}
        b_conf = {f"b{i+1}": [b["query"], b["size"], b["count"], b["type"]] for i, b in enumerate(backgrounds)} if backgrounds else {}
        t_path = _save_fetcher_config(CONFIG_DIR / "_ai_t_conf.json", t_conf)
        b_path = _save_fetcher_config(CONFIG_DIR / "_ai_b_conf.json", b_conf)
        project_name = Prompt.ask("Project name", default=f"ai_ncbi_{datetime.now().strftime('%Y%m%d_%H%M%S')}")
        out_dir = RUNS_DIR / project_name
        out_dir.mkdir(parents=True, exist_ok=True)
        cmd = [
            sys.executable, "-m", "rational_design.cli", "pipeline",
            "--out", str(out_dir),
            "--email", email,
            "--target_config", str(t_path),
            "--bg_config", str(b_path),
            "--params", str(params_path),
            "--language", lang,
        ]
        if ai_base_url:
            cmd.extend(["--ai_base_url", ai_base_url])
        if ai_model:
            cmd.extend(["--ai_model", ai_model])
        _run_subprocess(cmd, "\U0001F9EC AI-Proposed NCBI Pipeline")
        _show_run_summary(out_dir)

    elif action == "propose_validation":
        primers = proposal.get("primers", [])
        if not primers:
            console.print("[red]\u274C Proposal missing primers.[/]")
            return
        csv_path = CONFIG_DIR / "_ai_primers.csv"
        with open(csv_path, "w", newline="") as f:
            import csv as csv_writer
            w = csv_writer.writer(f)
            w.writerow(["name", "forward", "reverse"])
            for p in primers:
                w.writerow([p.get("name", "Primer"), p.get("fwd", ""), p.get("rev", "")])
        console.print(f"[dim]   Saved {len(primers)} primers to {csv_path}[/]")
        output_path = Prompt.ask("Output CSV", default="PCR_Advanced_Report.csv")
        cmd = [
            sys.executable, "-m", "rational_design.cli", "validate_primers",
            "-c", str(csv_path),
            "-t", str(CONFIG_DIR),
            "-b", str(CONFIG_DIR),
            "-o", output_path,
        ]
        _run_subprocess(cmd, "\U0001F52C AI-Proposed Validation")

    elif action == "propose_multiplex":
        if not email:
            console.print("[red]\u274C NCBI email required for multiplex proposals.[/]")
            return
        targets = _items_from_proposal(proposal.get("targets") or proposal.get("target") or [], default_count=500)
        backgrounds = _items_from_proposal(proposal.get("backgrounds") or proposal.get("background") or [], default_count=50)
        if len(targets) < 2:
            console.print("[red]\u274C At least 2 targets required for multiplex.[/]")
            return
        t_conf = {f"t{i+1}": [t["query"], t["size"], t["count"], t["type"]] for i, t in enumerate(targets)}
        b_conf = {f"b{i+1}": [b["query"], b["size"], b["count"], b["type"]] for i, b in enumerate(backgrounds)} if backgrounds else {}
        project_name = Prompt.ask("Project name", default=f"ai_multiplex_{datetime.now().strftime('%Y%m%d_%H%M%S')}")
        out_dir = RUNS_DIR / project_name
        out_dir.mkdir(parents=True, exist_ok=True)
        shared_context_path = out_dir / "shared_primer_context.json"
        t_paths: list[str] = []
        for idx, (key, val) in enumerate(t_conf.items(), 1):
            t_dir = out_dir / "0_raw_data" / f"t{idx}"
            t_dir.mkdir(parents=True, exist_ok=True)
            single_t_conf = {key: val}
            st_path = CONFIG_DIR / f"_ai_mux_t{idx}_conf.json"
            _save_fetcher_config(st_path, single_t_conf)
            fetch_cmd = [
                sys.executable, "-m", "rational_design.cli", "pipeline",
                "--out", str(out_dir / f"_fetch_t{idx}"),
                "--email", email,
                "--target_config", str(st_path),
                "--bg_config", str(st_path),
            ]
            _run_subprocess(fetch_cmd, f"\U0001F9EC Fetching target #{idx}")
            t_paths.append(str(t_dir))

            target_folder = out_dir / f"t{idx}"
            target_folder.mkdir(parents=True, exist_ok=True)
            virtual_bg = out_dir / "0_raw_data" / f"bg_for_t{idx}"
            virtual_bg.mkdir(parents=True, exist_ok=True)
            for jdx, tp in enumerate(t_paths):
                if jdx + 1 == idx:
                    continue
                _copy_sequence_files(Path(tp), virtual_bg, prefix=f"target_bg_t{jdx+1}_")

            cmd = [
                sys.executable, "-m", "rational_design.cli", "pipeline",
                "--out", str(target_folder),
                "--local_target", str(t_dir),
                "--local_bg", str(virtual_bg),
                "--params", str(params_path),
                "--language", lang,
                "--shared_context", str(shared_context_path),
            ]
            if ai_base_url:
                cmd.extend(["--ai_base_url", ai_base_url])
            if ai_model:
                cmd.extend(["--ai_model", ai_model])
            _run_subprocess(cmd, f"\U0001F9EC Target #{idx} (AI Proposal)")
            _accumulate_shared_context(target_folder, shared_context_path)

        target_folders = [out_dir / f"t{idx}" for idx in range(1, len(targets) + 1)]
        assay_type = Prompt.ask("Assay type", choices=["qPCR", "Conventional"], default="qPCR")
        mux_cmd = [
            sys.executable, "-m", "rational_design.cli", "multiplex_analyze",
            "--folders", *[str(f) for f in target_folders],
            "--out", str(out_dir),
            "--assay_type", assay_type,
            "--language", lang,
        ]
        if ai_base_url:
            mux_cmd.extend(["--ai_base_url", ai_base_url])
        if ai_model:
            mux_cmd.extend(["--ai_model", ai_model])
        _run_subprocess(mux_cmd, "\U0001F9EA Multiplex Analysis")
        _show_run_summary(out_dir)

    else:
        console.print(f"[yellow]\u26A0\uFE0F Unknown action: {action}[/]")

    _wait_and_return()


def _items_from_proposal(items: Any, default_count: int) -> list[dict[str, Any]]:
    if not items:
        return []
    if isinstance(items, dict):
        items = [items]
    parsed = []
    for item in items:
        if not isinstance(item, dict):
            continue
        query = str(item.get("query", "")).strip()
        if not query:
            continue
        parsed.append({
            "query": query,
            "size": float(item.get("size", 0.0) or 0.0),
            "count": int(item.get("count", default_count) or default_count),
            "type": str(item.get("type", "genome") or "genome"),
        })
    return parsed


# ─────────────────────────────────────────────
#  AI EXPERT EVALUATION
# ─────────────────────────────────────────────
def _select_run_or_path() -> Path | None:
    choice = Prompt.ask("Select source", choices=["list", "path"], default="list")
    if choice == "path":
        p = _prompt_path("Run directory path", must_exist=True)
        return Path(p) if p else None

    if not RUNS_DIR.exists():
        console.print("[yellow]No runs found. Please enter a path manually.[/]")
        p = _prompt_path("Run directory path", must_exist=True)
        return Path(p) if p else None

    proj_dirs = sorted(RUNS_DIR.iterdir(), key=lambda p: p.stat().st_mtime, reverse=True)
    proj_dirs = [d for d in proj_dirs if d.is_dir()][:20]
    if not proj_dirs:
        console.print("[yellow]No runs found.[/]")
        return None

    console.print("\n[bold]Recent projects:[/]")
    for idx, d in enumerate(proj_dirs, 1):
        run_folders = sorted(d.iterdir(), key=lambda p: p.stat().st_mtime, reverse=True)
        latest = ""
        if run_folders:
            latest = datetime.fromtimestamp(run_folders[0].stat().st_mtime).strftime("%Y-%m-%d %H:%M")
        console.print(f"  [{idx}] {d.name}  [dim]{latest}[/]")

    idx_str = Prompt.ask("Select project number", default="1")
    try:
        proj = proj_dirs[int(idx_str) - 1]
    except (ValueError, IndexError):
        return None

    sub_dirs = sorted(proj.iterdir(), key=lambda p: p.stat().st_mtime, reverse=True)
    sub_dirs = [d for d in sub_dirs if d.is_dir()]
    if not sub_dirs:
        return proj

    if len(sub_dirs) == 1:
        return sub_dirs[0]

    console.print(f"\n[bold]Runs in {proj.name}:[/]")
    for idx, d in enumerate(sub_dirs, 1):
        ts = datetime.fromtimestamp(d.stat().st_mtime).strftime("%Y-%m-%d %H:%M:%S")
        console.print(f"  [{idx}] {d.name}  [dim]{ts}[/]")
    idx_str = Prompt.ask("Select run number", default="1")
    try:
        return sub_dirs[int(idx_str) - 1]
    except (ValueError, IndexError):
        return sub_dirs[0]


def _display_run_results(run_dir: Path) -> None:
    console.print(f"\n[bold cyan]\U0001F4C2 {run_dir}[/]")

    final_csv = run_dir / "FINAL_ASSAY.csv"
    audit_path = run_dir / "audit_trail.json"
    report_path = run_dir / "ai_expert_report.json"

    if final_csv.exists():
        try:
            import pandas as pd
            df = pd.read_csv(final_csv)
            console.print(f"\n[green]\u2705 FINAL_ASSAY.csv[/] \u2014 {len(df)} assay(s)")
            table = Table(box=box.SIMPLE, border_style="green")
            cols = [c for c in df.columns[:6]]
            for c in cols:
                table.add_column(str(c), style="cyan")
            for _, row in df.head(8).iterrows():
                table.add_row(*[str(row[c])[:35] for c in cols])
            console.print(table)
        except Exception:
            console.print("[yellow]   Could not parse FINAL_ASSAY.csv[/]")
    else:
        console.print("[yellow]   No FINAL_ASSAY.csv found[/]")

    if audit_path.exists():
        console.print("   [green]\u2705 audit_trail.json found[/]")

    if report_path.exists():
        try:
            with open(report_path, "r", encoding="utf-8") as f:
                report = json.load(f)
            console.print(f"\n[bold]\U0001F916 Existing AI Report[/]")
            console.print(f"   Verdict: [bold]{report.get('overall_verdict', 'N/A')}[/]")
            console.print(f"   Best assay: [green]{report.get('best_assay_name', 'N/A')}[/]")
            rec = report.get('clinical_recommendation', '')
            if rec:
                console.print(f"   Recommendation: [dim]{rec[:300]}[/dim]")
            if not Confirm.ask("Re-run AI evaluation?", default=False):
                return
        except Exception:
            pass


def _run_ai_evaluation(run_dir: Path, ai_base_url: str, ai_model: str) -> None:
    from .ai_core import LocalBackend, AssayEvaluator

    val_dir = run_dir / "3_validation"
    stats_csv = val_dir / "master_pcr_results_summary.csv"
    if not stats_csv.exists():
        stats_csv = val_dir / "pcr_results_summary.csv"
    report_path = run_dir / "ai_expert_report.json"

    if not stats_csv.exists():
        console.print("[yellow]   No validation results found. Cannot run AI evaluation.[/]")
        return

    console.print(f"\n[bold]\U0001F916 Consulting AI Expert...[/]")
    params_path = run_dir.parent.parent / "parameters.json"
    current_params = {}
    try:
        with open(params_path, "r", encoding="utf-8") as f:
            current_params = json.load(f)
    except Exception:
        pass

    params_defaults = {
        "min_sensitivity": 90.0, "product_size_min": 100, "product_size_max": 350,
        "primer_tm_min": 55.0, "primer_tm_max": 65.0,
    }
    current_params = {**params_defaults, **current_params}

    with Progress(SpinnerColumn(), TextColumn("[progress.description]{task.description}"), console=console) as progress:
        progress.add_task("[cyan]AI evaluating primers...", total=None)

        analytical_report = generate_batch_analytical_summary(str(stats_csv), current_params, language=_lang())
        backend = LocalBackend(base_url=ai_base_url, model_name=ai_model)
        evaluator = AssayEvaluator(backend)
        report = evaluator.evaluate_candidates(analytical_report, language=_lang())

    if "error" in report:
        console.print(f"[red]\u274C AI Error: {report['error']}[/]")
        return

    with open(report_path, "w", encoding="utf-8") as f:
        json.dump(report, f, ensure_ascii=False, indent=4)

    console.print(f"\n[bold green]\u2705 AI Report saved to {report_path}[/]")
    console.print(f"   Verdict: [bold]{report.get('overall_verdict', 'N/A')}[/]")
    console.print(f"   Best assay: [green]{report.get('best_assay_name', 'N/A')}[/]")
    console.print(f"   Recommendation: [dim]{report.get('clinical_recommendation', '')[:400]}[/dim]")


# ─────────────────────────────────────────────
#  SUBPROCESS RUNNER
# ─────────────────────────────────────────────
def _run_subprocess(cmd: list[str], title: str) -> None:
    console.print(f"\n[bold]{title}[/]")
    console.print(f"   [dim]{' '.join(cmd)}[/]\n")

    my_env = os.environ.copy()
    my_env["PYTHONIOENCODING"] = "utf-8"
    my_env["PYTHONUNBUFFERED"] = "1"

    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(bar_width=None),
        TimeElapsedColumn(),
        console=console,
        transient=False,
    ) as progress:
        task = progress.add_task(f"[cyan]Running...", total=None)

        process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            shell=False,
            text=True,
            bufsize=1,
            encoding="utf-8",
            env=my_env,
        )
        ACTIVE_PROCESSES.append(process)

        lines: list[str] = []
        stages = {
            "STAGE 1": "\U0001F4E6 Building datasets",
            "CYCLE": "\U0001F504 Auto-optimization",
            "STAGE 2": "\U0001F9EC Designing primers",
            "STAGE 3": "\U0001F4BB In-silico validation",
            "STAGE 4": "\U0001F50D Probe design",
            "STAGE 5": "\U0001F310 BLAST annotation",
            "STAGE 6": "\U0001F916 AI consultation",
        }

        while True:
            line = process.stdout.readline()
            if not line and process.poll() is not None:
                break
            if line:
                clean = line.strip()
                lines.append(clean)
                for keyword, desc in stages.items():
                    if keyword in clean:
                        progress.update(task, description=f"[cyan]{desc}")
                        break

        if process in ACTIVE_PROCESSES:
            ACTIVE_PROCESSES.remove(process)

        return_code = process.poll() or 0

    if return_code == 0:
        console.print("[bold green]\u2705 Pipeline completed successfully![/]")
    else:
        console.print(f"[bold red]\u274C Pipeline failed (exit code {return_code})[/]")

    show = Confirm.ask("Show full log?", default=False)
    if show:
        console.print("\n[dim]" + "\n".join(lines[-100:]) + "[/]")


def _wait_and_return() -> None:
    Prompt.ask("\nPress Enter to return to menu", default="")


if __name__ == "__main__":
    main()
