#!/bin/bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$DIR"

echo "========================================================"
echo " 🧬 Starting Rational Primer Design..."
echo "========================================================"

if [ ! -d "venv" ]; then
    echo "❌ Error: Virtual environment (venv) not found."
    echo "   Please run './INSTALL_UX.sh' first."
    exit 1
fi

VENV_PYTHON="$DIR/venv/bin/python"

echo " 🚀 Launching Terminal Wizard..."
"$VENV_PYTHON" -m rational_design.cli term
