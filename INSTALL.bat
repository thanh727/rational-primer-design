@echo off
setlocal
cd /d "%~dp0"
title 🧬 Installing Rational Primer Design...

echo ========================================================
echo 🧬 RATIONAL PRIMER DESIGN - INSTALLER
echo ========================================================
echo.

python --version >nul 2>&1
if %errorlevel% neq 0 (
    echo [ERROR] Python is not installed or not in your PATH.
    echo Please install Python 3.9+ and try again.
    pause
    exit /b
)

if not exist "venv\" (
    echo [*] Creating virtual environment...
    python -m venv venv
)

echo [*] Upgrading pip...
call venv\Scripts\python.exe -m pip install --upgrade pip

echo [*] Installing dependencies and Rational Primer Design...
call venv\Scripts\python.exe -m pip install -e .

echo.
echo ========================================================
echo ✅ INSTALLATION COMPLETE!
echo ========================================================
echo You can now use "RUN_APP.bat" to start the tool.
echo.
pause