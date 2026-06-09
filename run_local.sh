#!/usr/bin/env bash
set -euo pipefail

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$DIR"

echo "========================================================"
echo "  🧬 Running Rational Primer Design (Web App)"
echo "========================================================"

if [ ! -d "venv" ]; then
    echo "❌ Error: Virtual environment (venv) not found."
    echo "   Please run './INSTALL_UX.sh' first."
    exit 1
fi

VENV_PYTHON="$DIR/venv/bin/python"

if [ ! -d "frontend/node_modules" ]; then
    echo "[*] Installing frontend dependencies..."
    cd "$DIR/frontend"
    npm install
    cd "$DIR"
fi

cleanup() {
    echo ""
    echo "[*] Shutting down..."
    if [ -n "${API_PID:-}" ] && kill -0 "$API_PID" >/dev/null 2>&1; then
        kill "$API_PID" >/dev/null 2>&1 || true
        wait "$API_PID" >/dev/null 2>&1 || true
    fi
    echo "✅ Done."
}
trap cleanup EXIT INT TERM

echo "[*] Starting backend API server (port 8000)..."
"$VENV_PYTHON" -m uvicorn rational_design.api:app --host 127.0.0.1 --port 8000 --reload &
API_PID=$!

sleep 3

echo "[*] Starting frontend (port 3000)..."
cd "$DIR/frontend"
NEXT_PUBLIC_API_BASE_URL="${NEXT_PUBLIC_API_BASE_URL:-http://127.0.0.1:8000}" npm run dev -- --hostname 127.0.0.1 --port "${RPD_WEB_PORT:-3000}"

echo ""
echo "========================================================"
echo "  Backend:  http://127.0.0.1:8000"
echo "  Frontend: http://localhost:3000"
echo "========================================================"
