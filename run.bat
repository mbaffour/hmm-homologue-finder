@echo off
REM HMM Homologue Finder - Windows launcher.
REM The bioinformatics tools (HMMER, MAFFT, etc.) have no native Windows builds,
REM so on Windows this pipeline runs inside WSL2 (Windows Subsystem for Linux).
setlocal

where wsl >nul 2>nul
if %errorlevel% neq 0 goto NOWSL

echo Launching the HMM Homologue Finder inside WSL2...
REM Translate this folder to a WSL path and run the Linux launcher there.
for /f "usebackq delims=" %%p in (`wsl wslpath "%~dp0"`) do set "TOOLDIR=%%p"
wsl bash -lc "cd '%TOOLDIR%' && bash run.sh"
goto END

:NOWSL
echo.
echo ============================================================
echo  WSL2 is required to run this pipeline on Windows.
echo  The bioinformatics tools are not available for native Windows.
echo.
echo  One-time setup (in an Administrator PowerShell):
echo      wsl --install -d Ubuntu
echo  Then reboot, open "Ubuntu" from the Start menu, and run:
echo      cd /mnt/c/path/to/HMM_Homologue_Finder
echo      bash run.sh
echo  See README.md for details.
echo ============================================================
echo.
pause

:END
endlocal
