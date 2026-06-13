from __future__ import annotations

import csv
import json
import os
import re
import shutil
import signal
import subprocess
import sys
import time
import uuid
import zipfile
from io import BytesIO
from urllib.error import URLError
from urllib.parse import urlparse
from urllib.request import Request, urlopen
from pathlib import Path
from typing import Any

import pandas as pd
import psutil
from fastapi import FastAPI, File, Form, HTTPException, UploadFile
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import Response
from pydantic import BaseModel, Field


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_PARAMS_PATH = REPO_ROOT / "config" / "parameters.json"
RUN_ROOT = Path(os.environ.get("RPD_RUN_ROOT", REPO_ROOT / "runs" / "web_jobs"))
PROCESS_REGISTRY: dict[str, subprocess.Popen] = {}
LAST_HEARTBEAT = time.time()
HAS_HEARTBEAT_STARTED = False



class LocalPipelineRequest(BaseModel):
    target_path: str = Field(..., min_length=1)
    background_path: str = Field(..., min_length=1)
    output_name: str | None = None
    parameters: dict[str, Any] = Field(default_factory=dict)
    language: str = "en"
    ai_base_url: str | None = None
    ai_model: str | None = None


class QueryItem(BaseModel):
    query: str = Field(..., min_length=1)
    size: float = 0.0
    count: int = 50
    type: str = "genome"


class AutoPipelineRequest(BaseModel):
    email: str = Field(..., min_length=3)
    targets: list[QueryItem] = Field(..., min_length=1)
    backgrounds: list[QueryItem] = Field(..., min_length=1)
    output_name: str | None = None
    parameters: dict[str, Any] = Field(default_factory=dict)
    language: str = "en"
    ai_base_url: str | None = None
    ai_model: str | None = None


class PrimerPair(BaseModel):
    name: str = Field(..., min_length=1)
    fwd: str = Field(..., min_length=8)
    rev: str = Field(..., min_length=8)
    template: str | None = None
    probe: str | None = None


class ValidationRequest(BaseModel):
    primers: list[PrimerPair] = Field(..., min_length=1)
    target_path: str | None = None
    background_path: str | None = None
    email: str | None = None
    targets: list[QueryItem] = Field(default_factory=list)
    backgrounds: list[QueryItem] = Field(default_factory=list)
    output_name: str | None = None
    extract_sequence: bool = False
    max_mismatch: int = Field(default=2, ge=0, le=12)
    workers: int = Field(default=0, ge=0, le=256)
    max_len: int = Field(default=1500, ge=50, le=10000)
    action: str = "validate"


class LocalMultiplexRequest(BaseModel):
    target_paths: list[str] = Field(..., min_length=2)
    background_path: str = Field(..., min_length=1)
    output_name: str | None = None
    parameters: dict[str, Any] = Field(default_factory=dict)
    language: str = "en"
    assay_type: str = "qPCR"
    ai_base_url: str | None = None
    ai_model: str | None = None


class AutoMultiplexRequest(BaseModel):
    email: str = Field(..., min_length=3)
    targets: list[QueryItem] = Field(..., min_length=2)
    backgrounds: list[QueryItem] = Field(default_factory=list)
    output_name: str | None = None
    parameters: dict[str, Any] = Field(default_factory=dict)
    language: str = "en"
    assay_type: str = "qPCR"
    ai_base_url: str | None = None
    ai_model: str | None = None


class MultiplexAnalyzeRequest(BaseModel):
    folders: list[str] = Field(..., min_length=2)
    output_name: str | None = None
    language: str = "en"
    assay_type: str = "qPCR"
    ai_base_url: str | None = None
    ai_model: str | None = None


class ChatMessage(BaseModel):
    role: str
    content: str


class ChatRequest(BaseModel):
    messages: list[ChatMessage] = Field(..., min_length=1)
    language: str = "en"
    ai_base_url: str = "http://localhost:11434/v1"
    ai_model: str = "llama3"
    email: str | None = None
    output_name: str | None = None
    parameters: dict[str, Any] = Field(default_factory=dict)
    auto_run: bool = True
    expert_context: dict[str, Any] | None = None


class AiModelsResponse(BaseModel):
    base_url: str
    provider: str
    available: bool
    models: list[str] = Field(default_factory=list)
    error: str | None = None


class FileBrowserEntry(BaseModel):
    name: str
    path: str
    is_dir: bool
    size_bytes: int | None = None


class FileBrowserResponse(BaseModel):
    path: str
    parent: str | None
    entries: list[FileBrowserEntry]
    roots: list[str] = Field(default_factory=list)


class JobResponse(BaseModel):
    id: str
    status: str
    source: str | None = None
    output_dir: str
    created_at: float
    updated_at: float
    command: list[str]
    return_code: int | None = None


def create_app() -> FastAPI:
    app = FastAPI(
        title="Rational Primer Design API",
        version="1.0.3",
        description="FastAPI orchestration layer for the primer design pipeline.",
    )
    app.add_middleware(
        CORSMiddleware,
        allow_origins=[
            "http://localhost:3000",
            "http://127.0.0.1:3000",
            "http://localhost:3001",
            "http://127.0.0.1:3001",
        ],
        allow_credentials=True,
        allow_methods=["*"],
        allow_headers=["*"],
    )

    @app.get("/")
    def root() -> dict[str, str]:
        return {"service": "Rational Primer Design API", "health": "/api/health"}

    @app.post("/api/heartbeat")
    def post_heartbeat() -> dict[str, str]:
        global LAST_HEARTBEAT, HAS_HEARTBEAT_STARTED
        LAST_HEARTBEAT = time.time()
        HAS_HEARTBEAT_STARTED = True
        return {"status": "alive"}

    def monitor_heartbeat():
        global LAST_HEARTBEAT, HAS_HEARTBEAT_STARTED
        import threading
        # Initial startup grace period
        time.sleep(40)
        while True:
            time.sleep(2)
            # If there are active running subprocesses, reset heartbeat to avoid timeout
            active_jobs = any(proc.poll() is None for proc in PROCESS_REGISTRY.values())
            if active_jobs:
                LAST_HEARTBEAT = time.time()
                continue

            if HAS_HEARTBEAT_STARTED and (time.time() - LAST_HEARTBEAT > 300):
                print("⚠️ [SYSTEM] Heartbeat lost. Closing application and cleaning resources...")
                # 1. Kill any running subprocesses in PROCESS_REGISTRY
                for job_id, proc in list(PROCESS_REGISTRY.items()):
                    if proc.poll() is None:
                        try:
                            _terminate_process(proc)
                        except Exception:
                            pass
                
                # 2. Terminate backend API server itself
                os.kill(os.getpid(), signal.SIGTERM)
                break

    if "pytest" not in sys.modules:
        import threading
        threading.Thread(target=monitor_heartbeat, daemon=True).start()

    @app.get("/api/health")
    def health() -> dict[str, Any]:
        return {"status": "ok", "run_root": str(_ensure_run_root())}

    @app.get("/api/status")
    def status(ai_base_url: str = "http://localhost:11434/v1") -> dict[str, Any]:
        jobs = list_jobs()
        counts: dict[str, int] = {}
        for job in jobs:
            counts[job["status"]] = counts.get(job["status"], 0) + 1
        return {
            "status": "ok",
            "run_root": str(_ensure_run_root()),
            "file_browser_roots": [str(root) for root in _file_browser_roots()],
            "job_counts": counts,
            "job_count": len(jobs),
            "run_root_size_bytes": _dir_size(_ensure_run_root()),
            "disk": _disk_usage(_ensure_run_root()),
            "ai": _detect_offline_models(ai_base_url),
        }

    @app.get("/api/default-parameters")
    def default_parameters() -> dict[str, Any]:
        return _read_default_parameters()

    @app.get("/api/ai/models", response_model=AiModelsResponse)
    def list_ai_models(base_url: str = "http://localhost:11434/v1") -> dict[str, Any]:
        return _detect_offline_models(base_url)

    @app.get("/api/files/browse", response_model=FileBrowserResponse)
    def browse_files(path: str | None = None, kind: str = "any") -> dict[str, Any]:
        return _browse_files(path, kind)

    @app.get("/api/jobs", response_model=list[JobResponse])
    def list_jobs() -> list[dict[str, Any]]:
        jobs = []
        for meta_path in sorted(_ensure_run_root().glob("*/job.json"), reverse=True):
            try:
                jobs.append(_refresh_status(_read_json(meta_path)))
            except Exception:
                continue
        return jobs

    @app.get("/api/history")
    def list_history() -> dict[str, Any]:
        return {"legacy_runs": _legacy_history(), "web_jobs": list_jobs()}

    @app.post("/api/jobs/local", response_model=JobResponse)
    def create_local_job(request: LocalPipelineRequest) -> dict[str, Any]:
        target_path = Path(request.target_path).expanduser().resolve()
        background_path = Path(request.background_path).expanduser().resolve()
        _validate_input_path(target_path, "target_path")
        _validate_input_path(background_path, "background_path")

        job_id, job_dir, output_dir = _new_job_dirs(request.output_name)
        params_path = _write_parameters(job_dir, request.parameters)
        command = _build_pipeline_command(
            output_dir=output_dir,
            params_path=params_path,
            language=request.language,
            local_target=target_path,
            local_bg=background_path,
            ai_base_url=request.ai_base_url,
            ai_model=request.ai_model,
        )
        return _start_job(job_id, job_dir, output_dir, command, source="local")

    @app.post("/api/jobs/auto", response_model=JobResponse)
    def create_auto_job(request: AutoPipelineRequest) -> dict[str, Any]:
        _validate_email(request.email)
        job_id, job_dir, output_dir = _new_job_dirs(request.output_name or request.targets[0].query)
        params_path = _write_parameters(job_dir, request.parameters)
        target_config = _write_query_config(job_dir / "target_config.json", request.targets, "t")
        background_config = _write_query_config(job_dir / "background_config.json", request.backgrounds, "b")
        command = _build_auto_pipeline_command(
            output_dir=output_dir,
            params_path=params_path,
            language=request.language,
            email=request.email,
            target_config=target_config,
            background_config=background_config,
            ai_base_url=request.ai_base_url,
            ai_model=request.ai_model,
        )
        return _start_job(job_id, job_dir, output_dir, command, source="auto")

    @app.post("/api/jobs/validate", response_model=JobResponse)
    def create_validation_job(request: ValidationRequest) -> dict[str, Any]:
        return _create_validation_job(request)

    @app.post("/api/jobs/multiplex/local", response_model=JobResponse)
    def create_local_multiplex_job(request: LocalMultiplexRequest) -> dict[str, Any]:
        return _create_local_multiplex_job(request)

    @app.post("/api/jobs/multiplex/auto", response_model=JobResponse)
    def create_auto_multiplex_job(request: AutoMultiplexRequest) -> dict[str, Any]:
        return _create_auto_multiplex_job(request)

    @app.post("/api/jobs/multiplex/analyze", response_model=JobResponse)
    def create_multiplex_analysis_job(request: MultiplexAnalyzeRequest) -> dict[str, Any]:
        return _create_multiplex_analysis_job(request)

    @app.post("/api/jobs/upload", response_model=JobResponse)
    async def create_upload_job(
        target_files: list[UploadFile] = File(...),
        background_files: list[UploadFile] = File(...),
        output_name: str | None = Form(None),
        parameters: str = Form("{}"),
        language: str = Form("en"),
        ai_base_url: str | None = Form(None),
        ai_model: str | None = Form(None),
    ) -> dict[str, Any]:
        try:
            parsed_params = json.loads(parameters) if parameters else {}
        except json.JSONDecodeError as exc:
            raise HTTPException(status_code=400, detail=f"parameters is not valid JSON: {exc}") from exc

        job_id, job_dir, output_dir = _new_job_dirs(output_name)
        upload_root = job_dir / "input"
        target_dir = upload_root / "target"
        background_dir = upload_root / "background"
        target_dir.mkdir(parents=True, exist_ok=True)
        background_dir.mkdir(parents=True, exist_ok=True)

        await _save_uploads(target_files, target_dir)
        await _save_uploads(background_files, background_dir)

        params_path = _write_parameters(job_dir, parsed_params)
        command = _build_pipeline_command(
            output_dir=output_dir,
            params_path=params_path,
            language=language,
            local_target=target_dir,
            local_bg=background_dir,
            ai_base_url=ai_base_url,
            ai_model=ai_model,
        )
        return _start_job(job_id, job_dir, output_dir, command, source="upload")

    @app.get("/api/jobs/{job_id}", response_model=JobResponse)
    def get_job(job_id: str) -> dict[str, Any]:
        return _get_job_or_404(job_id)

    @app.get("/api/jobs/{job_id}/logs")
    def get_logs(job_id: str, tail: int = 12000) -> dict[str, Any]:
        job = _get_job_or_404(job_id)
        process_log = _job_dir(job_id) / "process.log"
        text = _tail_file(process_log, max(1000, min(tail, 200000)))
        return {"job": job, "logs": text or ""}

    @app.post("/api/ai/chat")
    def ai_chat(request: ChatRequest) -> dict[str, Any]:
        try:
            from .ai_core import AIExpertAgent, LocalBackend
        except Exception as exc:
            raise HTTPException(status_code=500, detail=f"AI module unavailable: {exc}") from exc

        backend = LocalBackend(base_url=request.ai_base_url, model_name=request.ai_model)
        agent = AIExpertAgent(backend)
        messages = [{"role": m.role, "content": m.content} for m in request.messages]
        expert_context = request.expert_context
        raw_reply = agent.chat(messages, expert_report=expert_context, language=request.language)
        proposal = _extract_proposal(raw_reply)
        display_reply = _strip_proposal(raw_reply) if proposal else raw_reply

        job = None
        blocked_reason = None
        if request.auto_run and proposal and proposal.get("run_immediately", False):
            try:
                job = _start_job_from_ai_proposal(proposal, request)
            except HTTPException as exc:
                blocked_reason = str(exc.detail)

        return {
            "reply": display_reply.strip(),
            "raw_reply": raw_reply,
            "proposal": proposal,
            "job": job,
            "blocked_reason": blocked_reason,
        }

    @app.get("/api/jobs/{job_id}/results")
    def get_results(job_id: str) -> dict[str, Any]:
        job = _get_job_or_404(job_id)
        output_dir = Path(job["output_dir"])
        report_dir = output_dir / "4_validation_report"
        designed_probes_csv = output_dir / "designed_probes.csv"
        report_csv = designed_probes_csv if designed_probes_csv.exists() else (output_dir / "PCR_Advanced_Report.csv")
        return {
            "job": job,
            "summary": _read_json_if_exists(report_dir / "validation_summary.json"),
            "audit": _read_json_if_exists(output_dir / "audit_trail.json"),
            "cross_contamination": _read_json_if_exists(report_dir / "cross_contamination_report.json"),
            "ai_report": _read_json_if_exists(output_dir / "ai_expert_report.json"),
            "final_assays": _read_csv_preview(output_dir / "FINAL_ASSAY.csv"),
            "candidate_summary": _read_csv_preview(output_dir / "3_validation" / "master_pcr_results_summary.csv"),
            "known_primer_validation": _read_csv_preview(report_csv),
            "known_primer_summary": _read_json_if_exists(report_dir / "per_primer_summary.json"),
            "multiplex_kits": _read_csv_preview(output_dir / "MULTIPLEX_KITS.csv"),
            "multiplex_details": _read_json_if_exists(output_dir / "multiplex_details.json"),
            "files": _result_files(output_dir),
        }

    @app.post("/api/jobs/{job_id}/cancel", response_model=JobResponse)
    def cancel_job(job_id: str) -> dict[str, Any]:
        job = _get_job_or_404(job_id)
        proc = PROCESS_REGISTRY.get(job_id)
        if proc and proc.poll() is None:
            _terminate_process(proc)
            job["status"] = "cancelled"
            job["return_code"] = proc.poll()
            job["updated_at"] = time.time()
            _write_json(_job_dir(job_id) / "job.json", job)
        return job

    @app.get("/api/jobs/{job_id}/archive")
    def archive_job(job_id: str) -> Response:
        job = _get_job_or_404(job_id)
        job_dir = _job_dir(job_id)
        buffer = BytesIO()
        with zipfile.ZipFile(buffer, "w", zipfile.ZIP_DEFLATED) as archive:
            for path in job_dir.rglob("*"):
                if path.is_file():
                    archive.write(path, path.relative_to(job_dir))
        buffer.seek(0)
        filename = f"{_safe_name(Path(job['output_dir']).name or job_id)}-{job_id}.zip"
        return Response(
            content=buffer.getvalue(),
            media_type="application/zip",
            headers={"Content-Disposition": f'attachment; filename="{filename}"'},
        )

    @app.delete("/api/jobs/{job_id}")
    def delete_job(job_id: str) -> dict[str, Any]:
        job = _get_job_or_404(job_id)
        if job["status"] == "running":
            proc = PROCESS_REGISTRY.get(job_id)
            if proc and proc.poll() is None:
                _terminate_process(proc)
        shutil.rmtree(_job_dir(job_id), ignore_errors=True)
        PROCESS_REGISTRY.pop(job_id, None)
        return {"deleted": True, "job_id": job_id}

    return app


def _ensure_run_root() -> Path:
    RUN_ROOT.mkdir(parents=True, exist_ok=True)
    return RUN_ROOT


def _read_default_parameters() -> dict[str, Any]:
    if not DEFAULT_PARAMS_PATH.exists():
        return {}
    return _read_json(DEFAULT_PARAMS_PATH)


def _browse_files(path: str | None, kind: str) -> dict[str, Any]:
    roots = _file_browser_roots()
    current = Path(path or roots[0]).expanduser()
    if not current.exists():
        current = current.parent if current.parent.exists() else roots[0]
    if current.is_file():
        current = current.parent
    current = current.resolve()
    if not _is_under_any_root(current, roots):
        raise HTTPException(status_code=403, detail=f"Path is outside allowed browser roots: {current}")

    entries: list[dict[str, Any]] = []
    try:
        children = sorted(current.iterdir(), key=lambda item: (not item.is_dir(), item.name.lower()))
    except OSError as exc:
        raise HTTPException(status_code=400, detail=f"Cannot browse path: {exc}") from exc

    include_files = kind in {"any", "file"}
    include_dirs = kind in {"any", "directory"}
    for child in children[:400]:
        try:
            is_dir = child.is_dir()
            if not ((is_dir and include_dirs) or ((not is_dir) and include_files) or is_dir):
                continue
            stat = child.stat()
            entries.append(
                {
                    "name": child.name,
                    "path": str(child),
                    "is_dir": is_dir,
                    "size_bytes": None if is_dir else stat.st_size,
                }
            )
        except OSError:
            continue

    parent = str(current.parent) if current.parent != current and _is_under_any_root(current.parent, roots) else None
    return {"path": str(current), "parent": parent, "entries": entries, "roots": [str(root) for root in roots]}


def _file_browser_roots() -> list[Path]:
    configured = os.environ.get("RPD_FILE_BROWSER_ROOTS") or os.environ.get("RPD_FILE_BROWSER_ROOT")
    if configured:
        raw_roots = [item for item in re.split(r"[:;]", configured) if item.strip()]
    else:
        raw_roots = [str(REPO_ROOT), str(Path.home())]
    roots = []
    for raw in raw_roots:
        path = Path(raw).expanduser().resolve()
        if path.exists() and path.is_dir() and path not in roots:
            roots.append(path)
    return roots or [REPO_ROOT]


def _is_under_any_root(path: Path, roots: list[Path]) -> bool:
    resolved = path.resolve()
    return any(resolved == root or root in resolved.parents for root in roots)


def _dir_size(path: Path) -> int:
    total = 0
    if not path.exists():
        return 0
    for item in path.rglob("*"):
        try:
            if item.is_file():
                total += item.stat().st_size
        except OSError:
            continue
    return total


def _disk_usage(path: Path) -> dict[str, int]:
    usage = shutil.disk_usage(path)
    return {"total": usage.total, "used": usage.used, "free": usage.free}


def _detect_offline_models(base_url: str) -> dict[str, Any]:
    normalized = (base_url or "http://localhost:11434/v1").rstrip("/")
    parsed = urlparse(normalized)
    if parsed.scheme not in {"http", "https"} or not parsed.netloc:
        return {
            "base_url": normalized,
            "provider": "unknown",
            "available": False,
            "models": [],
            "error": "Invalid model server URL.",
        }

    errors = []
    try:
        from openai import OpenAI
        client = OpenAI(api_key="sk-local", base_url=normalized, timeout=5.0)
        models = sorted({m.id for m in client.models.list().data})
        return {
            "base_url": normalized,
            "provider": "openai-compatible",
            "available": True,
            "models": models,
            "error": None,
        }
    except Exception as exc:
        errors.append(f"openai: {exc}")

    root = normalized[:-3] if normalized.endswith("/v1") else normalized
    try:
        from urllib.request import Request, urlopen
        from urllib.error import URLError
        request = Request(f"{root}/api/tags", headers={"Accept": "application/json"})
        with urlopen(request, timeout=5.0) as response:
            payload = json.loads(response.read().decode("utf-8"))
        models = sorted({str(item.get("name") or item.get("model")) for item in payload.get("models", []) if item.get("name") or item.get("model")})
        return {
            "base_url": normalized,
            "provider": "ollama",
            "available": True,
            "models": models,
            "error": None,
        }
    except Exception as exc:
        errors.append(f"ollama: {exc}")

    return {
        "base_url": normalized,
        "provider": "unknown",
        "available": False,
        "models": [],
        "error": "; ".join(errors[-2:]) if errors else "No local model server responded.",
    }


def _new_job_dirs(output_name: str | None) -> tuple[str, Path, Path]:
    stem = _safe_name(output_name or f"design-{time.strftime('%Y%m%d-%H%M%S')}")
    job_id = f"{time.strftime('%Y%m%d-%H%M%S')}-{uuid.uuid4().hex[:8]}"
    job_dir = _ensure_run_root() / job_id
    output_dir = job_dir / stem
    output_dir.mkdir(parents=True, exist_ok=True)
    return job_id, job_dir, output_dir


def _job_dir(job_id: str) -> Path:
    if "/" in job_id or "\\" in job_id or ".." in job_id:
        raise HTTPException(status_code=404, detail="Job not found")
    return _ensure_run_root() / job_id


def _safe_name(value: str) -> str:
    cleaned = "".join(ch if ch.isalnum() or ch in ("-", "_", ".") else "-" for ch in value.strip())
    cleaned = cleaned.strip(".-_")
    return cleaned[:80] or "design"


def _validate_input_path(path: Path, field_name: str) -> None:
    if not path.exists():
        raise HTTPException(status_code=400, detail=f"{field_name} does not exist: {path}")
    if not (path.is_dir() or path.suffix.lower() in {".fasta", ".fa", ".fna", ".fas"}):
        raise HTTPException(status_code=400, detail=f"{field_name} must be a FASTA file or directory: {path}")


def _validate_email(email: str) -> None:
    if "@" not in email or "." not in email:
        raise HTTPException(status_code=400, detail="A valid NCBI email address is required.")


def _write_parameters(job_dir: Path, overrides: dict[str, Any]) -> Path:
    params = _read_default_parameters()
    params.update(overrides or {})
    params_path = job_dir / "parameters.json"
    _write_json(params_path, params)
    return params_path


def _write_query_config(path: Path, items: list[QueryItem], prefix: str) -> Path:
    config = {}
    for index, item in enumerate(items, start=1):
        if not item.query.strip():
            continue
        config[f"{prefix}{index}"] = [
            item.query.strip(),
            float(item.size or 0.0),
            int(item.count or 0),
            item.type or "genome",
        ]
    if not config:
        raise HTTPException(status_code=400, detail="At least one non-empty query is required.")
    _write_json(path, config)
    return path


def _build_pipeline_command(
    *,
    output_dir: Path,
    params_path: Path,
    language: str,
    local_target: Path,
    local_bg: Path,
    ai_base_url: str | None,
    ai_model: str | None,
) -> list[str]:
    command = [
        sys.executable,
        "-m",
        "rational_design.cli",
        "pipeline",
        "--out",
        str(output_dir),
        "--local_target",
        str(local_target),
        "--local_bg",
        str(local_bg),
        "--params",
        str(params_path),
        "--language",
        language if language in {"vi", "en"} else "en",
    ]
    if ai_base_url:
        command.extend(["--ai_base_url", ai_base_url])
    if ai_model:
        command.extend(["--ai_model", ai_model])
    return command


def _build_auto_pipeline_command(
    *,
    output_dir: Path,
    params_path: Path,
    language: str,
    email: str,
    target_config: Path,
    background_config: Path,
    ai_base_url: str | None,
    ai_model: str | None,
) -> list[str]:
    command = [
        sys.executable,
        "-m",
        "rational_design.cli",
        "pipeline",
        "--out",
        str(output_dir),
        "--email",
        email,
        "--target_config",
        str(target_config),
        "--bg_config",
        str(background_config),
        "--params",
        str(params_path),
        "--language",
        language if language in {"vi", "en"} else "en",
    ]
    if ai_base_url:
        command.extend(["--ai_base_url", ai_base_url])
    if ai_model:
        command.extend(["--ai_model", ai_model])
    return command


def _create_validation_job(request: ValidationRequest) -> dict[str, Any]:
    _validate_primers(request.primers)
    job_id, job_dir, output_dir = _new_job_dirs(request.output_name or "primer-validation")
    primers_path = _write_primers_csv(job_dir / "primers.csv", request.primers)

    if request.target_path:
        target_path = Path(request.target_path).expanduser().resolve()
        _validate_input_path(target_path, "target_path")
        background_path = None
        if request.background_path:
            background_path = Path(request.background_path).expanduser().resolve()
            _validate_input_path(background_path, "background_path")

        if request.action == "design-probes":
            command = [
                sys.executable,
                "-u",
                "-m",
                "rational_design.cli",
                "design_probes",
                "-c",
                str(primers_path),
                "-t",
                str(target_path),
                "-o",
                str(output_dir / "designed_probes.csv"),
                "-e",
                str(request.max_mismatch),
                "--max_len",
                str(request.max_len),
            ]
            if background_path:
                command.extend(["-b", str(background_path)])
        else:
            command = _build_validate_primers_command(
                primers_path=primers_path,
                target_path=target_path,
                background_path=background_path,
                output_csv=output_dir / "PCR_Advanced_Report.csv",
                extract_sequence=request.extract_sequence,
                max_mismatch=request.max_mismatch,
                workers=request.workers,
                max_len=request.max_len,
            )
        return _start_job(job_id, job_dir, output_dir, command, source="validate-local")

    if not request.email:
        raise HTTPException(status_code=400, detail="Validation needs either local target_path or NCBI email + target/background queries.")
    _validate_email(request.email)
    if not request.targets:
        raise HTTPException(status_code=400, detail="At least one target query is required for online validation.")

    request_path = _write_web_job_request(
        job_dir,
        {
            "output_dir": str(output_dir),
            "primers_path": str(primers_path),
            "email": request.email,
            "targets": [_query_item_payload(item) for item in request.targets],
            "backgrounds": [_query_item_payload(item) for item in request.backgrounds],
            "extract_sequence": request.extract_sequence,
            "max_mismatch": request.max_mismatch,
            "workers": request.workers,
            "max_len": request.max_len,
            "action": request.action,
        },
    )
    command = _build_web_job_command("validate_online", request_path)
    return _start_job(job_id, job_dir, output_dir, command, source="validate-auto")


def _create_local_multiplex_job(request: LocalMultiplexRequest) -> dict[str, Any]:
    target_paths = [Path(path).expanduser().resolve() for path in request.target_paths]
    for index, target_path in enumerate(target_paths, start=1):
        _validate_input_path(target_path, f"target_paths[{index}]")
    background_path = Path(request.background_path).expanduser().resolve()
    _validate_input_path(background_path, "background_path")

    job_id, job_dir, output_dir = _new_job_dirs(request.output_name or "local-multiplex")
    params_path = _write_parameters(job_dir, request.parameters)
    request_path = _write_web_job_request(
        job_dir,
        {
            "output_dir": str(output_dir),
            "target_paths": [str(path) for path in target_paths],
            "background_path": str(background_path),
            "params_path": str(params_path),
            "language": request.language,
            "assay_type": _normalize_assay_type(request.assay_type),
            "ai_base_url": request.ai_base_url,
            "ai_model": request.ai_model,
        },
    )
    return _start_job(job_id, job_dir, output_dir, _build_web_job_command("local_multiplex", request_path), source="multiplex-local")


def _create_auto_multiplex_job(request: AutoMultiplexRequest) -> dict[str, Any]:
    _validate_email(request.email)
    job_id, job_dir, output_dir = _new_job_dirs(request.output_name or "auto-multiplex")
    params_path = _write_parameters(job_dir, request.parameters)
    request_path = _write_web_job_request(
        job_dir,
        {
            "output_dir": str(output_dir),
            "email": request.email,
            "targets": [_query_item_payload(item) for item in request.targets],
            "backgrounds": [_query_item_payload(item) for item in request.backgrounds],
            "params_path": str(params_path),
            "language": request.language,
            "assay_type": _normalize_assay_type(request.assay_type),
            "ai_base_url": request.ai_base_url,
            "ai_model": request.ai_model,
        },
    )
    return _start_job(job_id, job_dir, output_dir, _build_web_job_command("auto_multiplex", request_path), source="multiplex-auto")


def _create_multiplex_analysis_job(request: MultiplexAnalyzeRequest) -> dict[str, Any]:
    folders = [Path(folder).expanduser().resolve() for folder in request.folders]
    for index, folder in enumerate(folders, start=1):
        if not folder.exists() or not folder.is_dir():
            raise HTTPException(status_code=400, detail=f"folders[{index}] must be an existing result folder: {folder}")

    job_id, job_dir, output_dir = _new_job_dirs(request.output_name or "multiplex-analysis")
    command = [
        sys.executable,
        "-m",
        "rational_design.cli",
        "multiplex_analyze",
        "--folders",
        *[str(folder) for folder in folders],
        "--out",
        str(output_dir),
        "--language",
        request.language if request.language in {"vi", "en"} else "vi",
        "--assay_type",
        _normalize_assay_type(request.assay_type),
    ]
    if request.ai_base_url:
        command.extend(["--ai_base_url", request.ai_base_url])
    if request.ai_model:
        command.extend(["--ai_model", request.ai_model])
    return _start_job(job_id, job_dir, output_dir, command, source="multiplex-analyze")


def _build_validate_primers_command(
    *,
    primers_path: Path,
    target_path: Path,
    background_path: Path | None,
    output_csv: Path,
    extract_sequence: bool,
    max_mismatch: int,
    workers: int,
    max_len: int,
) -> list[str]:
    command = [
        sys.executable,
        "-u",
        "-m",
        "rational_design.cli",
        "validate_primers",
        "-c",
        str(primers_path),
        "-t",
        str(target_path),
        "-o",
        str(output_csv),
        "-e",
        str(max_mismatch),
        "-w",
        str(workers),
        "--max_len",
        str(max_len),
    ]
    if background_path:
        command.extend(["-b", str(background_path)])
    if extract_sequence:
        command.append("-s")
    return command


def _write_primers_csv(path: Path, primers: list[PrimerPair]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["name", "fwd", "rev", "template", "probe"])
        writer.writeheader()
        for primer in primers:
            writer.writerow(
                {
                    "name": primer.name.strip(),
                    "fwd": primer.fwd.strip().upper(),
                    "rev": primer.rev.strip().upper(),
                    "template": (primer.template or "").strip(),
                    "probe": (primer.probe or "").strip().upper(),
                }
            )
    return path


def _validate_primers(primers: list[PrimerPair]) -> None:
    allowed = set("ATGCRYSWKMNBDHV")
    for primer in primers:
        for field_name, sequence in (("fwd", primer.fwd), ("rev", primer.rev)):
            clean = sequence.strip().upper()
            if len(clean) < 8:
                raise HTTPException(status_code=400, detail=f"{primer.name}.{field_name} must be at least 8 bases.")
            bad = sorted({base for base in clean if base not in allowed})
            if bad:
                raise HTTPException(status_code=400, detail=f"{primer.name}.{field_name} contains unsupported bases: {''.join(bad)}")
        
        if primer.probe:
            clean_probe = primer.probe.strip().upper()
            if clean_probe:
                if len(clean_probe) < 8:
                    raise HTTPException(status_code=400, detail=f"{primer.name}.probe must be at least 8 bases.")
                bad_probe = sorted({base for base in clean_probe if base not in allowed})
                if bad_probe:
                    raise HTTPException(status_code=400, detail=f"{primer.name}.probe contains unsupported bases: {''.join(bad_probe)}")


def _write_web_job_request(job_dir: Path, payload: dict[str, Any]) -> Path:
    request_path = job_dir / "request.json"
    _write_json(request_path, payload)
    return request_path


def _build_web_job_command(command: str, request_path: Path) -> list[str]:
    return [sys.executable, "-m", "rational_design.web_jobs", command, "--request", str(request_path)]


def _query_item_payload(item: QueryItem) -> dict[str, Any]:
    return {"query": item.query, "size": item.size, "count": item.count, "type": item.type}


def _normalize_assay_type(value: str) -> str:
    return value if value in {"qPCR", "Conventional"} else "qPCR"


def _start_job_from_ai_proposal(proposal: dict[str, Any], request: ChatRequest) -> dict[str, Any]:
    action = proposal.get("action")
    params = dict(request.parameters or {})
    for key in ("design_target_sampling_size", "design_background_sampling_size", "validation_target_sampling_size", "validation_background_sampling_size", "degenerate_primers", "max_iupac_per_primer"):
        if key in proposal:
            params[key] = proposal[key]

    if action == "propose_local_design":
        target_path = Path(str(proposal.get("local_target", ""))).expanduser().resolve()
        background_path = Path(str(proposal.get("local_background", ""))).expanduser().resolve()
        _validate_input_path(target_path, "local_target")
        _validate_input_path(background_path, "local_background")
        job_id, job_dir, output_dir = _new_job_dirs(request.output_name or target_path.name)
        params_path = _write_parameters(job_dir, params)
        command = _build_pipeline_command(
            output_dir=output_dir,
            params_path=params_path,
            language=request.language,
            local_target=target_path,
            local_bg=background_path,
            ai_base_url=request.ai_base_url,
            ai_model=request.ai_model,
        )
        return _start_job(job_id, job_dir, output_dir, command, source="ai-local")

    if action == "propose_design":
        if not request.email:
            raise HTTPException(status_code=400, detail="AI proposed an NCBI design, but no email was provided.")
        _validate_email(request.email)
        targets = _items_from_proposal(proposal.get("target") or proposal.get("targets") or [], default_count=500)
        backgrounds = _items_from_proposal(proposal.get("background") or proposal.get("backgrounds") or [], default_count=50)
        job_id, job_dir, output_dir = _new_job_dirs(request.output_name or (targets[0].query if targets else "ai-design"))
        params_path = _write_parameters(job_dir, params)
        target_config = _write_query_config(job_dir / "target_config.json", targets, "t")
        background_config = _write_query_config(job_dir / "background_config.json", backgrounds, "b")
        command = _build_auto_pipeline_command(
            output_dir=output_dir,
            params_path=params_path,
            language=request.language,
            email=request.email,
            target_config=target_config,
            background_config=background_config,
            ai_base_url=request.ai_base_url,
            ai_model=request.ai_model,
        )
        return _start_job(job_id, job_dir, output_dir, command, source="ai-auto")

    if action == "propose_validation":
        primers = _primers_from_proposal(proposal)
        targets = _items_from_proposal(proposal.get("target") or proposal.get("targets") or [], default_count=50)
        backgrounds = _items_from_proposal(proposal.get("background") or proposal.get("backgrounds") or [], default_count=10)
        validation_request = ValidationRequest(
            primers=primers,
            email=request.email,
            targets=targets,
            backgrounds=backgrounds,
            output_name=request.output_name or (targets[0].query if targets else "ai-validation"),
            max_mismatch=int(proposal.get("max_mismatch", 2) or 2),
            max_len=int(proposal.get("max_len", 1500) or 1500),
        )
        return _create_validation_job(validation_request)

    if action == "propose_multiplex":
        if not request.email:
            raise HTTPException(status_code=400, detail="AI proposed a multiplex NCBI design, but no email was provided.")
        targets = _items_from_proposal(proposal.get("targets") or proposal.get("target") or [], default_count=500)
        backgrounds = _items_from_proposal(proposal.get("background") or proposal.get("backgrounds") or [], default_count=50)
        multiplex_request = AutoMultiplexRequest(
            email=request.email,
            targets=targets,
            backgrounds=backgrounds,
            output_name=request.output_name or "ai-multiplex",
            parameters=params,
            language=request.language,
            assay_type=str(proposal.get("assay_type", "qPCR") or "qPCR"),
            ai_base_url=request.ai_base_url,
            ai_model=request.ai_model,
        )
        return _create_auto_multiplex_job(multiplex_request)

    raise HTTPException(status_code=400, detail=f"AI proposal action is not runnable in web mode yet: {action}")


def _items_from_proposal(items: list[dict[str, Any]], default_count: int) -> list[QueryItem]:
    parsed = []
    for item in items:
        query = str(item.get("query", "")).strip()
        if not query:
            continue
        parsed.append(
            QueryItem(
                query=query,
                size=float(item.get("size", 0.0) or 0.0),
                count=int(item.get("count", default_count) or default_count),
                type=str(item.get("type", "genome") or "genome"),
            )
        )
    if not parsed:
        raise HTTPException(status_code=400, detail="AI proposal did not contain runnable target/background queries.")
    return parsed


def _primers_from_proposal(proposal: dict[str, Any]) -> list[PrimerPair]:
    primers = []
    for item in proposal.get("primers", []):
        if not isinstance(item, dict):
            continue
        primers.append(
            PrimerPair(
                name=str(item.get("name", "Primer")).strip() or "Primer",
                fwd=str(item.get("fwd", "")).strip(),
                rev=str(item.get("rev", "")).strip(),
                template=str(item.get("template", "")).strip() or None,
            )
        )
    if not primers:
        raise HTTPException(status_code=400, detail="AI proposal did not include primers to validate.")
    return primers


def _extract_proposal(text: str) -> dict[str, Any] | None:
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


def _strip_proposal(text: str) -> str:
    def replace_block(match: re.Match) -> str:
        return "" if '"action"' in match.group(1) else match.group(0)

    return re.sub(r"```json\s*([\s\S]*?)\s*```", replace_block, text).strip()


def _start_job(job_id: str, job_dir: Path, output_dir: Path, command: list[str], source: str) -> dict[str, Any]:
    process_log = open(job_dir / "process.log", "a", encoding="utf-8")
    env = os.environ.copy()
    env["PYTHONUNBUFFERED"] = "1"
    try:
        proc = subprocess.Popen(
            command,
            cwd=str(REPO_ROOT),
            stdout=process_log,
            stderr=subprocess.STDOUT,
            text=True,
            start_new_session=os.name != "nt",
            env=env,
        )
    except Exception:
        process_log.close()
        raise
    process_log.close()

    PROCESS_REGISTRY[job_id] = proc
    now = time.time()
    meta = {
        "id": job_id,
        "status": "running",
        "source": source,
        "pid": proc.pid,
        "output_dir": str(output_dir),
        "created_at": now,
        "updated_at": now,
        "command": command,
        "return_code": None,
    }
    _write_json(job_dir / "job.json", meta)
    return meta


async def _save_uploads(files: list[UploadFile], destination: Path) -> None:
    if not files:
        raise HTTPException(status_code=400, detail="At least one FASTA file is required.")
    seen: dict[str, int] = {}
    for upload in files:
        filename = _safe_name(Path(upload.filename or "sequence.fasta").name)
        if not filename.lower().endswith((".fasta", ".fa", ".fna", ".fas")):
            raise HTTPException(status_code=400, detail=f"Unsupported FASTA extension: {filename}")
        count = seen.get(filename, 0)
        seen[filename] = count + 1
        if count:
            stem, suffix = os.path.splitext(filename)
            filename = f"{stem}-{count + 1}{suffix}"
        target_path = destination / filename
        with open(target_path, "wb") as handle:
            while True:
                chunk = await upload.read(1024 * 1024)
                if not chunk:
                    break
                handle.write(chunk)


def _get_job_or_404(job_id: str) -> dict[str, Any]:
    meta_path = _job_dir(job_id) / "job.json"
    if not meta_path.exists():
        raise HTTPException(status_code=404, detail="Job not found")
    return _refresh_status(_read_json(meta_path))


def _refresh_status(job: dict[str, Any]) -> dict[str, Any]:
    job_id = job["id"]
    proc = PROCESS_REGISTRY.get(job_id)
    if proc:
        code = proc.poll()
        if code is None:
            return job
        job["return_code"] = code
        job["status"] = _infer_finished_status(Path(job["output_dir"]), code)
        job["updated_at"] = time.time()
        _write_json(_job_dir(job_id) / "job.json", job)
        PROCESS_REGISTRY.pop(job_id, None)
        return job

    if job.get("status") in {"pending", "running"}:
        pid = job.get("pid")
        if pid and psutil.pid_exists(int(pid)):
            return job
        job["status"] = _infer_finished_status(Path(job["output_dir"]), job.get("return_code"))
        job["updated_at"] = time.time()
        _write_json(_job_dir(job_id) / "job.json", job)
    return job


def _infer_finished_status(output_dir: Path, return_code: int | None) -> str:
    summary = _read_json_if_exists(output_dir / "4_validation_report" / "validation_summary.json")
    if summary:
        if isinstance(summary, dict) and "pipeline_success" in summary:
            return "completed" if summary.get("pipeline_success") else "failed"
        if return_code == 0 or (output_dir / "PCR_Advanced_Report.csv").exists():
            return "completed"
    if return_code is None:
        return "unknown"
    successful_outputs = [
        output_dir / "FINAL_ASSAY.csv",
        output_dir / "PCR_Advanced_Report.csv",
        output_dir / "MULTIPLEX_KITS.csv",
    ]
    return "completed" if return_code == 0 and any(path.exists() for path in successful_outputs) else "failed"


def _terminate_process(proc: subprocess.Popen) -> None:
    try:
        if os.name != "nt":
            os.killpg(proc.pid, signal.SIGTERM)
        else:
            proc.terminate()
        try:
            proc.wait(timeout=5)
        except subprocess.TimeoutExpired:
            if os.name != "nt":
                os.killpg(proc.pid, signal.SIGKILL)
            else:
                proc.kill()
    except ProcessLookupError:
        pass


def _tail_file(path: Path, limit: int) -> str:
    if not path.exists():
        return ""
    with open(path, "rb") as handle:
        handle.seek(0, os.SEEK_END)
        size = handle.tell()
        handle.seek(max(0, size - limit))
        return handle.read().decode("utf-8", errors="replace")


def _result_files(output_dir: Path) -> list[dict[str, Any]]:
    wanted = [
        output_dir / "FINAL_ASSAY.csv",
        output_dir / "pipeline_log.txt",
        output_dir / "audit_trail.json",
        output_dir / "ai_expert_report.json",
        output_dir / "PCR_Advanced_Report.csv",
        output_dir / "MULTIPLEX_KITS.csv",
        output_dir / "multiplex_details.json",
        output_dir / "4_validation_report" / "validation_summary.json",
        output_dir / "4_validation_report" / "cross_contamination_report.json",
        output_dir / "4_validation_report" / "per_primer_summary.json",
        output_dir / "4_validation_report" / "pcr_results_summary.csv",
        output_dir / "4_validation_report" / "target_hits.csv",
        output_dir / "4_validation_report" / "background_hits.csv",
    ]
    files = []
    for path in wanted:
        if path.exists():
            files.append({"name": path.name, "path": str(path), "size_bytes": path.stat().st_size})
    return files


def _legacy_history() -> list[dict[str, Any]]:
    try:
        from .analytics import ResultAnalysisEngine
    except Exception:
        return []

    try:
        runs = ResultAnalysisEngine(str(REPO_ROOT / "runs")).scan_historical_runs()
    except Exception:
        return []

    serialized = []
    for run in runs:
        item = dict(run)
        timestamp = item.get("timestamp")
        if hasattr(timestamp, "isoformat"):
            item["timestamp"] = timestamp.isoformat()
        serialized.append(item)
    return serialized


def _latest_expert_context() -> dict[str, Any] | None:
    for meta_path in sorted(_ensure_run_root().glob("*/job.json"), reverse=True):
        try:
            job = _refresh_status(_read_json(meta_path))
            output_dir = Path(job["output_dir"])
            report_dir = output_dir / "4_validation_report"
            context = {
                "job": job,
                "validation_report": _read_json_if_exists(report_dir / "validation_summary.json"),
                "cross_contamination_traceback": _read_json_if_exists(report_dir / "cross_contamination_report.json"),
                "ai_report": _read_json_if_exists(output_dir / "ai_expert_report.json"),
                "designed_assays": _read_csv_preview(output_dir / "FINAL_ASSAY.csv", limit=20),
                "known_primer_validation": _read_csv_preview(output_dir / "PCR_Advanced_Report.csv", limit=20),
                "known_primer_summary": _read_json_if_exists(report_dir / "per_primer_summary.json"),
                "multiplex_kits": _read_csv_preview(output_dir / "MULTIPLEX_KITS.csv", limit=20),
                "multiplex_details": _read_json_if_exists(output_dir / "multiplex_details.json"),
            }
            if any(value for key, value in context.items() if key != "job"):
                return context
        except Exception:
            continue
    return None


def _read_csv_preview(path: Path, limit: int = 50) -> dict[str, Any]:
    if not path.exists():
        return {"exists": False, "columns": [], "rows": []}
    try:
        df = pd.read_csv(path, nrows=limit)
        row_count = 0
        with open(path, "r", encoding="utf-8", errors="ignore") as f:
            for _ in f:
                row_count += 1
        if row_count > 0:
            row_count -= 1  # subtract header row
    except Exception as exc:
        return {"exists": True, "columns": [], "rows": [], "error": str(exc)}
    return {
        "exists": True,
        "columns": list(df.columns),
        "rows": df.fillna("").to_dict(orient="records"),
        "row_count": row_count,
    }


def _read_json_if_exists(path: Path) -> Any | None:
    if not path.exists():
        return None
    try:
        return _read_json(path)
    except Exception:
        return None


def _read_json(path: Path) -> Any:
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def _write_json(path: Path, data: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(data, handle, ensure_ascii=False, indent=2)


app = create_app()
