@echo off
echo.
echo  ===================================
echo   Drip Rate Estimator - Starting...
echo  ===================================
echo.

:: Activate conda environment - change 'base' to your env name if needed
call conda activate base 2>nul || call conda activate drip_rate 2>nul

:: Move to the app directory
cd /d "%~dp0"

:: Open browser after a short delay
start /b cmd /c "timeout /t 3 >nul && start http://localhost:5000"

:: Start the Flask app
python app.py

pause
