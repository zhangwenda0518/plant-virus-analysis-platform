@echo off
chcp 65001 >nul
title 环境自检（模块 / 工具 / 数据库 / 磁盘）
cd /d "%~dp0"
python main.py selfcheck
echo.
python tests\mem_check.py
pause
