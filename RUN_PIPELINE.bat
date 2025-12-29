@echo off
title Rational Primer Design Pipeline
:: Force UTF-8 Encoding for Emojis
chcp 65001 >nul
color 0A
cls

echo ========================================================
echo      🧬 RATIONAL PRIMER DESIGN PIPELINE 🧬
echo ========================================================
echo.

:: --- STEP 1: PROJECT NAME ---
set /p proj_name="> Enter Project Name (e.g. MyTest): "
echo.

:: --- STEP 2: CONFIGURATION ---
echo  [CONFIGURATION]
echo  -------------------------------------------------------
echo  Default config: config\parameters.json
echo  (Press Enter to use default, or drag a new JSON here)
echo  -------------------------------------------------------
set /p param_file="> Drag Parameter JSON here: "

:: Clean quotes and set default
if defined param_file set param_file=%param_file:"=%
if not defined param_file set param_file=config\parameters.json

echo.
echo  [INFO] Using Configuration: %param_file%
echo.

echo ========================================================
echo  Please select your operation mode:
echo.
echo    [1] FULL AUTOMATION (Download from NCBI)
echo    [2] LOCAL MODE (Use my own FASTA files)
echo.
echo ========================================================
set /p mode="Enter choice (1 or 2): "

if "%mode%"=="1" goto ONLINE
if "%mode%"=="2" goto LOCAL
goto END

:ONLINE
echo.
echo  [MODE: DOWNLOAD FROM NCBI]
echo  -------------------------------------------------------
echo  Tip: Drag and drop your JSON config files here.
echo  -------------------------------------------------------
echo.
set /p t_conf="> Drag & Drop your TARGET JSON file here: "
echo.
set /p b_conf="> Drag & Drop your BACKGROUND JSON file here: "
echo.

set t_conf=%t_conf:"=%
set b_conf=%b_conf:"=%

echo  🚀 Starting Pipeline...
rational-design pipeline --target_config "%t_conf%" --bg_config "%b_conf%" --email "nguyenthanh727@gmail.com" --out "%proj_name%" --params "%param_file%"
goto DONE

:LOCAL
echo.
echo  [MODE: LOCAL FILES]
echo  -------------------------------------------------------
echo  Default Target:     database\target
echo  Default Background: database\background
echo.
echo  (Press Enter to use defaults, or drag a specific folder)
echo  -------------------------------------------------------
echo.

:: TARGET FOLDER
set /p t_path="> Path to TARGET folder: "
if defined t_path set t_path=%t_path:"=%
if not defined t_path set t_path=database\target

:: BACKGROUND FOLDER
set /p b_path="> Path to BACKGROUND folder: "
if defined b_path set b_path=%b_path:"=%
if not defined b_path set b_path=database\background

echo.
echo  [INFO] Target Source:     %t_path%
echo  [INFO] Background Source: %b_path%
echo.

echo  🚀 Starting Pipeline using local data...
rational-design pipeline --local_target "%t_path%" --local_bg "%b_path%" --out "%proj_name%" --params "%param_file%"

:DONE
echo.
echo ========================================================
echo  ✅ DONE! Check the folder '%proj_name%' for 'FINAL_ASSAY.csv'.
echo ========================================================
pause
exit

:END
echo Invalid choice. Exiting.
pause