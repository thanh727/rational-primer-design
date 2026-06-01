#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$ROOT_DIR"

export RPD_API_HOST="${RPD_API_HOST:-127.0.0.1}"
export RPD_API_PORT="${RPD_API_PORT:-8000}"

python -m rational_design.api_server
