@echo off
chcp 65001 >nul
rem 将源码平台的测序数据目录挂载进 exe 平台（目录联接，不复制文件）
rem 用法：双击本文件；重复运行安全（已存在则跳过）
setlocal
set SRC=%~dp0fastq
set DST=%~dp0dist\VirusPlatform\fastq
if exist "%DST%" (
  echo 已存在: %DST%
) else (
  mklink /J "%DST%" "%SRC%" >nul
  echo 已挂载: %DST%  -^>  %SRC%
)
echo.
echo 之后在平台「分析管道」新建样品时，点 📁 进入 fastq 目录即可选择测序数据。
pause
