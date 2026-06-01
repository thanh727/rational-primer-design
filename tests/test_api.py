from pathlib import Path

from fastapi.testclient import TestClient

import rational_design.api as api_module


def test_api_default_parameters():
    client = TestClient(api_module.create_app())
    response = client.get("/api/default-parameters")

    assert response.status_code == 200
    assert "min_sensitivity" in response.json()


def test_api_rejects_missing_local_paths(tmp_path, monkeypatch):
    monkeypatch.setattr(api_module, "RUN_ROOT", tmp_path)
    client = TestClient(api_module.create_app())

    response = client.post(
        "/api/jobs/local",
        json={
            "target_path": str(tmp_path / "missing-target"),
            "background_path": str(tmp_path / "missing-bg"),
            "parameters": {},
        },
    )

    assert response.status_code == 400


def test_api_creates_local_job_without_shell(tmp_path, monkeypatch):
    repo_root = Path(__file__).resolve().parents[1]
    monkeypatch.setattr(api_module, "RUN_ROOT", tmp_path)

    def fake_start_job(job_id, job_dir, output_dir, command, source):
        return {
            "id": job_id,
            "status": "running",
            "source": source,
            "pid": 123,
            "output_dir": str(output_dir),
            "created_at": 1.0,
            "updated_at": 1.0,
            "command": command,
            "return_code": None,
        }

    monkeypatch.setattr(api_module, "_start_job", fake_start_job)
    client = TestClient(api_module.create_app())

    response = client.post(
        "/api/jobs/local",
        json={
            "target_path": str(repo_root / "test_data" / "target"),
            "background_path": str(repo_root / "test_data" / "background"),
            "output_name": "api-test",
            "parameters": {"enable_blast": False},
        },
    )

    payload = response.json()
    assert response.status_code == 200
    assert payload["status"] == "running"
    assert payload["command"][1:4] == ["-m", "rational_design.cli", "pipeline"]
    assert "--local_target" in payload["command"]


def test_api_creates_auto_job_from_keyword_configs(tmp_path, monkeypatch):
    monkeypatch.setattr(api_module, "RUN_ROOT", tmp_path)

    def fake_start_job(job_id, job_dir, output_dir, command, source):
        return {
            "id": job_id,
            "status": "running",
            "source": source,
            "pid": 123,
            "output_dir": str(output_dir),
            "created_at": 1.0,
            "updated_at": 1.0,
            "command": command,
            "return_code": None,
        }

    monkeypatch.setattr(api_module, "_start_job", fake_start_job)
    client = TestClient(api_module.create_app())

    response = client.post(
        "/api/jobs/auto",
        json={
            "email": "user@example.com",
            "targets": [{"query": "Salmonella enterica", "count": 500, "size": 0.0, "type": "genome"}],
            "backgrounds": [{"query": "Escherichia coli", "count": 50, "size": 0.0, "type": "genome"}],
            "parameters": {"enable_blast": False},
        },
    )

    payload = response.json()
    assert response.status_code == 200
    assert payload["source"] == "auto"
    assert "--target_config" in payload["command"]
    assert "--bg_config" in payload["command"]


def test_api_creates_known_primer_validation_job(tmp_path, monkeypatch):
    repo_root = Path(__file__).resolve().parents[1]
    monkeypatch.setattr(api_module, "RUN_ROOT", tmp_path)

    def fake_start_job(job_id, job_dir, output_dir, command, source):
        return {
            "id": job_id,
            "status": "running",
            "source": source,
            "pid": 123,
            "output_dir": str(output_dir),
            "created_at": 1.0,
            "updated_at": 1.0,
            "command": command,
            "return_code": None,
        }

    monkeypatch.setattr(api_module, "_start_job", fake_start_job)
    client = TestClient(api_module.create_app())

    response = client.post(
        "/api/jobs/validate",
        json={
            "primers": [{"name": "M1", "fwd": "ATGCGTACGTAC", "rev": "CGTACGTACGTA"}],
            "target_path": str(repo_root / "test_data" / "target"),
            "background_path": str(repo_root / "test_data" / "background"),
        },
    )

    payload = response.json()
    assert response.status_code == 200
    assert payload["source"] == "validate-local"
    assert "validate_primers" in payload["command"]
    assert "PCR_Advanced_Report.csv" in " ".join(payload["command"])


def test_api_creates_online_validation_web_job(tmp_path, monkeypatch):
    monkeypatch.setattr(api_module, "RUN_ROOT", tmp_path)

    def fake_start_job(job_id, job_dir, output_dir, command, source):
        return {
            "id": job_id,
            "status": "running",
            "source": source,
            "pid": 123,
            "output_dir": str(output_dir),
            "created_at": 1.0,
            "updated_at": 1.0,
            "command": command,
            "return_code": None,
        }

    monkeypatch.setattr(api_module, "_start_job", fake_start_job)
    client = TestClient(api_module.create_app())

    response = client.post(
        "/api/jobs/validate",
        json={
            "email": "user@example.com",
            "primers": [{"name": "M1", "fwd": "ATGCGTACGTAC", "rev": "CGTACGTACGTA"}],
            "targets": [{"query": "Target species", "count": 50, "size": 0.0, "type": "genome"}],
            "backgrounds": [{"query": "Background species", "count": 10, "size": 0.0, "type": "genome"}],
        },
    )

    payload = response.json()
    assert response.status_code == 200
    assert payload["source"] == "validate-auto"
    assert payload["command"][1:4] == ["-m", "rational_design.web_jobs", "validate_online"]


def test_api_creates_local_multiplex_web_job(tmp_path, monkeypatch):
    repo_root = Path(__file__).resolve().parents[1]
    monkeypatch.setattr(api_module, "RUN_ROOT", tmp_path)

    def fake_start_job(job_id, job_dir, output_dir, command, source):
        return {
            "id": job_id,
            "status": "running",
            "source": source,
            "pid": 123,
            "output_dir": str(output_dir),
            "created_at": 1.0,
            "updated_at": 1.0,
            "command": command,
            "return_code": None,
        }

    monkeypatch.setattr(api_module, "_start_job", fake_start_job)
    client = TestClient(api_module.create_app())

    response = client.post(
        "/api/jobs/multiplex/local",
        json={
            "target_paths": [str(repo_root / "test_data" / "target"), str(repo_root / "test_data" / "target")],
            "background_path": str(repo_root / "test_data" / "background"),
            "parameters": {"enable_blast": False},
            "assay_type": "qPCR",
        },
    )

    payload = response.json()
    assert response.status_code == 200
    assert payload["source"] == "multiplex-local"
    assert payload["command"][1:4] == ["-m", "rational_design.web_jobs", "local_multiplex"]


def test_extract_proposal_handles_nested_json_fence():
    text = """
    Tôi sẽ cấu hình.
    ```json
    {
      "action": "propose_design",
      "run_immediately": true,
      "target": [{"query": "A", "count": 500, "size": 0.0, "type": "genome"}],
      "background": [{"query": "B", "count": 50, "size": 0.0, "type": "genome"}]
    }
    ```
    """

    proposal = api_module._extract_proposal(text)

    assert proposal is not None
    assert proposal["target"][0]["query"] == "A"
    assert "```json" not in api_module._strip_proposal(text)
