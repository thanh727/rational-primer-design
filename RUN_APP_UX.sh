#!/bin/bash

# Get the directory where this script is located
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$DIR"

echo "========================================================"
echo " 🧬 Starting Rational Primer Design..."
echo "========================================================"

# 1. Check if venv exists
if [ ! -d "venv" ]; then
    echo "❌ Error: Virtual environment not found."
    echo "   Please run './INSTALL_UX.sh' first."
    exit 1
fi

# 2. Activate the Virtual Environment
source venv/bin/activate

# 3. Check for Streamlit
if ! python -c "import streamlit" &> /dev/null; then
    echo " [WARNING] Streamlit not found in venv. Installing..."
    pip install streamlit
fi

# 4. Run the App
echo " 🚀 Launching GUI..."
streamlit run gui.py --server.headless true