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

API_PID=""
FE_PID=""
NO_CONN_COUNT=0

cleanup() {
    echo ""
    echo "[*] Shutting down..."
    if [ -n "$FE_PID" ] && kill -0 "$FE_PID" >/dev/null 2>&1; then
        kill "$FE_PID" >/dev/null 2>&1 || true
        wait "$FE_PID" >/dev/null 2>&1 || true
    fi
    if [ -n "$API_PID" ] && kill -0 "$API_PID" >/dev/null 2>&1; then
        kill "$API_PID" >/dev/null 2>&1 || true
        wait "$API_PID" >/dev/null 2>&1 || true
    fi
    echo "✅ Done."
}
trap cleanup EXIT INT TERM

echo "[*] Starting backend API server (port 8000)..."
"$VENV_PYTHON" -m uvicorn rational_design.api:app --host 127.0.0.1 --port 8000 --reload &
API_PID=$!

# Wait for API to be ready
echo "[*] Waiting for API..."
for i in $(seq 1 30); do
    if curl -sf http://127.0.0.1:8000/api/health >/dev/null 2>&1; then
        echo "   ✅ API ready."
        break
    fi
    if [ "$i" -eq 30 ]; then
        echo "   ❌ API failed to start."
        exit 1
    fi
    sleep 1
done

echo "[*] Starting frontend (port 3000)..."
cd "$DIR/frontend"
NEXT_PUBLIC_API_BASE_URL="${NEXT_PUBLIC_API_BASE_URL:-http://127.0.0.1:8000}" \
  npm run dev -- --hostname 127.0.0.1 --port "${RPD_WEB_PORT:-3000}" &
FE_PID=$!
cd "$DIR"

# Wait for frontend to be ready
echo "[*] Waiting for frontend..."
for i in $(seq 1 60); do
    if curl -sf http://127.0.0.1:3000 >/dev/null 2>&1; then
        echo "   ✅ Frontend ready."
        break
    fi
    if [ "$i" -eq 60 ]; then
        echo "   ⚠️ Frontend did not respond in time, continuing..."
    fi
    sleep 2
done

# Open browser
echo "[*] Opening browser..."
if command -v open >/dev/null 2>&1; then
    open http://localhost:3000
elif command -v xdg-open >/dev/null 2>&1; then
    xdg-open http://localhost:3000
fi

echo ""
echo "========================================================"
echo "  Backend:  http://127.0.0.1:8000"
echo "  Frontend: http://localhost:3000"
echo "========================================================"
echo "  Close the browser tab or press Ctrl+C to stop."
echo ""

# Give the browser a moment to open and establish connections
sleep 5

# Monitor for browser activity — auto-shutdown when tab is closed.
# Checks every 10 seconds if there are any ESTABLISHED TCP
# connections to the frontend port from a browser.
echo "[*] Monitoring browser connections (close tab to stop)..."
while true; do
    if ! kill -0 "$FE_PID" >/dev/null 2>&1; then
        echo "[*] Frontend process died. Shutting down..."
        break
    fi
    if ! kill -0 "$API_PID" >/dev/null 2>&1; then
        echo "[*] Backend process died. Shutting down..."
        break
    fi

    # Count ESTABLISHED TCP connections to port 3000.
    # lsof returns header + 1 line per connection; returns 0 with no matches.
    connections=$(lsof -iTCP:"${RPD_WEB_PORT:-3000}" -sTCP:ESTABLISHED 2>/dev/null | awk 'NR>1' | wc -l | tr -d ' ' || true)
    if [ "$connections" -eq 0 ]; then
        NO_CONN_COUNT=$((NO_CONN_COUNT + 1))
        # Shut down after ~30 seconds of no browser activity
        if [ "$NO_CONN_COUNT" -ge 3 ]; then
            echo ""
            echo "[*] Browser tab closed. Shutting down servers..."
            break
        fi
    else
        NO_CONN_COUNT=0
    fi

    sleep 10
done
