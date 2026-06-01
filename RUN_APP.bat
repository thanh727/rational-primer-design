@echo off
setlocal
cd /d "%~dp0"

echo ========================================================
echo 🧬 Starting Rational Primer Design (Dark Mode)...
echo ========================================================

if not exist "venv\" (
    echo [ERROR] Virtual environment not found. Please run INSTALL.bat first.
    pause
    exit /b
)

set PYTHON_EXE=venv\Scripts\python.exe

echo Launching GUI...
%PYTHON_EXE% -m streamlit run rational_design/gui.py --theme.base="dark"

if %errorlevel% neq 0 (
    echo.
    echo [ERROR] The app crashed or closed unexpectedly.
)

pause