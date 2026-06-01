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

source venv/bin/activate

echo " [2/4] Upgrading pip..."
pip install --upgrade pip

echo " [3/4] Installing libraries..."
if [ -f "requirements.txt" ]; then
    pip install -r requirements.txt
else
    echo " [WARNING] requirements.txt not found!"
fi

echo " [4/4] Registering 'rational-design' command..."
pip install -e .

echo "========================================================"
echo "   ✅ INSTALLATION COMPLETE!"
echo "   You can now run './RUN_APP_UX.sh'"
echo "========================================================"