#!/bin/bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$DIR"

echo "========================================================"
echo "   🧬 SETTING UP ENVIRONMENT (Virtual Environment)"
echo "========================================================"

if [ ! -d "venv" ]; then
    echo " [1/4] Creating virtual environment (venv)..."
    python3 -m venv venv
else
    echo " [1/4] Virtual environment already exists."
fi

VENV_PYTHON="$DIR/venv/bin/python"
if [ ! -x "$VENV_PYTHON" ]; then
    echo " [ERROR] Virtual environment Python not found at $VENV_PYTHON"
    echo "         Remove the existing venv folder and run this installer again."
    exit 1
fi

echo " [2/4] Upgrading pip..."
"$VENV_PYTHON" -m pip install --upgrade pip

echo " [3/4] Installing libraries..."
if [ -f "requirements.txt" ]; then
    "$VENV_PYTHON" -m pip install -r requirements.txt
else
    echo " [WARNING] requirements.txt not found!"
fi

echo " [4/4] Registering 'rational-design' command..."
"$VENV_PYTHON" -m pip install -e .

echo "========================================================"
echo "   ✅ INSTALLATION COMPLETE!"
echo "   You can now run './RUN_APP_UX.sh'"
echo "========================================================"
