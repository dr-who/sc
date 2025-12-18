@echo off
REM Generate hourly billing data for Windows
REM Usage: gen_hourly_billing.bat [num_rows]
REM Output: datasets\hourlybilling.csv

setlocal

set TARGET=%1
if "%TARGET%"=="" set TARGET=100000

if not exist datasets mkdir datasets

echo Generating %TARGET% rows of SaaS billing data...

REM Use built-in sc to generate data
bin\sc.exe -e "genhourlybilling(%TARGET%)" > nul 2>&1
if exist datasets\hourlybilling.csv (
    echo Generated datasets\hourlybilling.csv
    for %%A in (datasets\hourlybilling.csv) do echo Size: %%~zA bytes
) else (
    echo Error: Generation failed
    echo Make sure bin\sc.exe exists
)
