@echo off
setlocal
title 🧬 Installing Rational Primer Design...

echo ========================================================
echo 🧬 RATIONAL PRIMER DESIGN - INSTALLER
echo ========================================================
echo.

:: 1. Check for Python
python --version >nul 2>&1
if %errorlevel% neq 0 (
    echo [ERROR] Python is not installed or not in your PATH.
    echo Please install Python 3.9+ and try again.
    pause
    exit /b
)

:: 2. Upgrade PIP (Good practice)
echo [*] Upgrading pip...
python -m pip install --upgrade pip

:: 3. Install the Tool (Editable Mode)
:: This replaces the manual command "pip install -e ."
echo.
echo [*] Installing Rational Primer Design package...
python -m pip install -e .

:: 4. Install GUI Dependencies explicitly (just to be safe)
echo.
echo [*] Verifying GUI dependencies...
python -m pip install streamlit tkinter

echo.
echo ========================================================
echo ✅ INSTALLATION COMPLETE!
echo ========================================================
echo You can now use "run_app.bat" to start the tool.
echo.
pause