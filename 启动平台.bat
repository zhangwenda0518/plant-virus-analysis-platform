@echo off
chcp 65001 >nul
title 植物病毒分析平台
cd /d "%~dp0"
echo.
echo   ╔══════════════════════════════════════════╗
echo   ║        植物病毒分析平台 正在启动...        ║
echo   ╚══════════════════════════════════════════╝
echo.
python app.py
pause
