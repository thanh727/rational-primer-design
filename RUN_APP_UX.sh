#!/bin/bash

echo "========================================================"
echo "🧬 Starting Rational Primer Design (Linux/macOS)..."
echo "========================================================"

# 1. Detect Python Executable (prefer python3, fallback to python)
if command -v python3 &>/dev/null; then
    PYTHON_EXE=python3
else
    PYTHON_EXE=python
fi

# 2. Check if Streamlit is installed
if ! $PYTHON_EXE -c "import streamlit" &> /dev/null; then
    echo "[WARNING] Streamlit not found. Installing..."
    $PYTHON_EXE -m pip install streamlit
fi

# 3. Run the App
# We use the same flags: Dark Mode + Run via Module
echo "Launching GUI..."
$PYTHON_EXE -m streamlit run gui.py --theme.base="dark"