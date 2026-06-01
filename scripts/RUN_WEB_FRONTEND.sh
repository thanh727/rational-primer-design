#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$ROOT_DIR/frontend"

export NEXT_PUBLIC_API_BASE_URL="${NEXT_PUBLIC_API_BASE_URL:-http://127.0.0.1:8000}"

if [ ! -d "node_modules" ]; then
  npm install
fi

npm run dev -- --hostname 127.0.0.1 --port "${RPD_WEB_PORT:-3000}"
