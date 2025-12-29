@echo off
setlocal EnableDelayedExpansion
chcp 65001 >nul

cls
echo ========================================================
echo    RATIONAL PRIMER DESIGN PIPELINE (Windows)
echo ========================================================
echo.

:: --- PROJECT NAME ---
set /p proj_name="> Enter Project Name (e.g. MyTest): "

:: --- CONFIGURATION ---
echo.
echo  [CONFIGURATION]
echo  Default: config/parameters.json
set /p param_file="> Drag Parameter JSON here (or press Enter): "

:: Clean quotes
if not defined param_file set param_file=config/parameters.json
set param_file=!param_file:"=!
if "%param_file%"=="" set param_file=config/parameters.json

:: --- MODE SELECTION ---
echo.
echo  [1] Download from NCBI
echo  [2] Local Files
echo.
set /p mode="Enter choice (1 or 2): "

if "%mode%"=="1" goto ONLINE
if "%mode%"=="2" goto LOCAL
goto END

:ONLINE
echo.
echo  [MODE: NCBI DOWNLOAD]
set /p t_conf="> Drag TARGET JSON: "
set /p b_conf="> Drag BACKGROUND JSON: "
set t_conf=!t_conf:"=!
set b_conf=!b_conf:"=!

rational-design pipeline --target_config "!t_conf!" --bg_config "!b_conf!" --out "!proj_name!" --params "!param_file!"
goto END

:LOCAL
echo.
echo  [MODE: LOCAL FILES]
echo  Default Target: database\target
set /p t_path="> Path to TARGET folder: "
set /p b_path="> Path to BACKGROUND folder: "

if not defined t_path set t_path=database\target
set t_path=!t_path:"=!
if "%t_path%"=="" set t_path=database\target

if not defined b_path set b_path=database\background
set b_path=!b_path:"=!
if "%b_path%"=="" set b_path=database\background

rational-design pipeline --local_target "!t_path!" --local_bg "!b_path!" --out "!proj_name!" --params "!param_file!"

:END
echo.
echo  DONE. Check output folder.
pause