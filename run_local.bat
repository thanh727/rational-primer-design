@echo off
setlocal
cd /d "%~dp0"
title 🧬 Rational Primer Design - Web App

echo ========================================================
echo   Running Rational Primer Design (Web App)
echo ========================================================

if not exist "venv\" (
    echo [ERROR] Virtual environment not found. Please run INSTALL.bat first.
    pause
    exit /b
)

if not exist "frontend\node_modules\" (
    echo [*] Installing frontend dependencies...
    cd /d "%~dp0frontend"
    call npm install
    if %errorlevel% neq 0 (
        echo [ERROR] npm install failed.
        pause
        exit /b
    )
    cd /d "%~dp0"
)

echo [*] Starting backend API server (port 8000)...
start "RPD Backend" cmd /c "cd /d "%~dp0" && venv\Scripts\python.exe -m uvicorn rational_design.api:app --host 127.0.0.1 --port 8000 --reload && pause"

timeout /t 3 /nobreak >nul

echo [*] Starting frontend (port 3000)...
start "RPD Frontend" cmd /c "cd /d "%~dp0frontend" && npm run dev && pause"

echo ========================================================
echo   Backend:  http://127.0.0.1:8000
echo   Frontend: http://localhost:3000
echo ========================================================
echo   Close both server windows to stop.
echo.
pause
