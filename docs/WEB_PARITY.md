# Web parity audit

This matrix tracks the Streamlit GUI capabilities that must remain available after the Next.js + FastAPI upgrade.

| Legacy Streamlit capability | Next.js route | FastAPI endpoint / runner | Status |
| --- | --- | --- | --- |
| Clear dashboard with run history | `/dashboard` | `/api/jobs`, `/api/history` | Preserved. Shows web jobs and legacy Streamlit runs from `runs/`. |
| Local single-target primer design | `/design/local` | `POST /api/jobs/local` | Preserved. Calls the existing `rational_design.cli pipeline` with local target/background inputs. |
| Automatic NCBI keyword design | `/design/auto` | `POST /api/jobs/auto` | Preserved. API writes target/background JSON configs and calls the existing pipeline fetcher path. |
| AI chat-driven design | `/ai` | `POST /api/ai/chat` | Preserved and restored. AI proposals can auto-run supported jobs. |
| AI proposal: `propose_design` | `/ai` | `POST /api/jobs/auto` through AI orchestration | Preserved. |
| AI proposal: `propose_local_design` | `/ai` | `POST /api/jobs/local` through AI orchestration | Preserved. |
| AI proposal: `propose_validation` | `/ai`, `/validate` | `POST /api/jobs/validate`, `rational_design.web_jobs validate_online` | Restored. Supports known primer validation with local folders or online NCBI fetch. |
| AI proposal: `propose_multiplex` | `/ai`, `/multiplex` | `POST /api/jobs/multiplex/auto`, `rational_design.web_jobs auto_multiplex` | Restored. Runs sequential per-target design with shared background and final multiplex analysis. |
| Local multiplex design | `/multiplex` | `POST /api/jobs/multiplex/local`, `rational_design.web_jobs local_multiplex` | Restored. Mirrors the old local multi-target flow with virtual cross-target backgrounds and shared primer context. |
| Multiplex analysis of existing completed folders | `/multiplex` | `POST /api/jobs/multiplex/analyze` | Preserved. Calls `rational_design.cli multiplex_analyze`. |
| Known-primer in-silico PCR validation | `/validate` | `POST /api/jobs/validate` | Preserved. Produces `PCR_Advanced_Report.csv` and `4_validation_report/`. |
| Results tables and output files | all workflow routes | `GET /api/jobs/{job_id}/results` | Expanded. Returns final assays, validation reports, multiplex kits, AI reports, and file paths. |
| Live logs and cancellation | all workflow routes | `GET /api/jobs/{job_id}/logs`, `POST /api/jobs/{job_id}/cancel` | Preserved. |
| Vietnamese / English UI | all routes | frontend i18n dictionary and pipeline `--language` | Preserved for route labels, forms, controls, status, and AI language. |

Implementation notes:

- Dashboard is now a first-class route at `/dashboard`; `/` redirects there.
- Workflow routes are separated so the UI can change without collapsing the old architecture: `/design/local`, `/design/auto`, `/ai`, `/validate`, and `/multiplex`.
- The website shell now uses top navigation for primary functions and a left shared-configuration sidebar for language, run name, NCBI email, AI backend, assay type, and core pipeline parameters.
- AI chat uses the restored legacy PCR/qPCR expert role. It can provide general PCR, qPCR, assay evaluation, wet-lab optimization, and troubleshooting advice without emitting runnable JSON unless the user requests an executable design/validation/multiplex job.
- Long multi-step web jobs are run by `rational_design.web_jobs` so FastAPI remains an orchestration layer and the existing CLI algorithms remain the source of truth.
