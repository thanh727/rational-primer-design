#!/bin/bash

echo "========================================================"
echo "   🧬 SETTING UP ENVIRONMENT (Virtual Environment)"
echo "========================================================"

# 1. Create a Virtual Environment if it doesn't exist
if [ ! -d "venv" ]; then
    echo " [1/4] Creating virtual environment (venv)..."
    python3 -m venv venv
else
    echo " [1/4] Virtual environment already exists."
fi

# 2. Activate the Virtual Environment
# This tells the shell to use the 'pip' and 'python' inside the folder
source venv/bin/activate

# 3. Upgrade pip inside the venv
echo " [2/4] Upgrading pip..."
pip install --upgrade pip

# 4. Install Dependencies
echo " [3/4] Installing libraries..."
if [ -f "requirements.txt" ]; then
    pip install -r requirements.txt
else
    echo " [WARNING] requirements.txt not found!"
fi

# 5. Install the package in editable mode
echo " [4/4] Registering 'rational-design' command..."
pip install -e .

echo "========================================================"
echo "   ✅ INSTALLATION COMPLETE!"
echo "   You can now run './RUN_APP_UX.sh'"
echo "========================================================"