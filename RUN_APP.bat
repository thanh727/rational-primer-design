@echo off
setlocal

echo ========================================================
echo 🧬 Starting Rational Primer Design (Dark Mode)...
echo ========================================================

:: --- CONFIGURATION ---
set PYTHON_EXE=python

:: OPTIONAL: Uncomment to force a specific path if "python" isn't found
:: set PYTHON_EXE="C:\Program Files\Python311\python.exe"

:: --- CHECK INSTALLATION ---
echo Checking for Streamlit...
%PYTHON_EXE% -c "import streamlit" >nul 2>&1
if %errorlevel% neq 0 (
    echo [WARNING] Streamlit not found. Attempting to install...
    %PYTHON_EXE% -m pip install streamlit
)

:: --- RUN THE APP IN DARK MODE ---
:: We added --theme.base="dark" to force the dark theme
echo Launching GUI...
%PYTHON_EXE% -m streamlit run gui.py --theme.base="dark"

if %errorlevel% neq 0 (
    echo.
    echo [ERROR] The app crashed or closed unexpectedly.
)

pause