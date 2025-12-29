#!/bin/bash
echo "========================================================"
echo "   🧬 SETTING UP ENVIRONMENT (Linux / macOS)"
echo "========================================================"
echo

# 1. Check for Python 3
if ! command -v python3 &> /dev/null; then
    echo "[ERROR] Python 3 is not installed."
    echo "  - macOS: Install from python.org"
    echo "  - Ubuntu/Debian: Run 'sudo apt install python3 python3-pip'"
    echo "  - Fedora: Run 'sudo dnf install python3'"
    exit 1
fi

echo " [1/3] Upgrading pip..."
python3 -m pip install --upgrade pip

echo
echo " [2/3] Installing libraries (including Turbo Levenshtein)..."
# Note: On some minimal Linux installs, you might need 'build-essential' 
# or 'python3-dev' if a pre-compiled wheel isn't available.
pip3 install -r requirements.txt

echo
echo " [3/3] Registering 'rational-design' command..."
pip3 install -e .

echo
echo "========================================================"
echo "   ✅ INSTALLATION COMPLETE!"
echo "   You can now run './RUN_PIPELINE.sh'"
echo "========================================================"