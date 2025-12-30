@echo off
echo ========================================================
echo 🧬 Starting Rational Primer Design - Desktop App...
echo ========================================================
echo.
echo Checking for updates...
pip install -r requirements.txt >nul 2>&1

echo Launching GUI...
streamlit run gui.py

pause