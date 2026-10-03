# Detached supervisor — survives Cursor closing. DO NOT run regen only inside Cursor terminals.
$ErrorActionPreference = "Stop"
$PythonDir = Split-Path -Parent $MyInvocation.MyCommand.Path
Set-Location $PythonDir

$StopFlag = "experiments\terminal_weighted\supervisor_stop_n24_rational.flag"
$Log = "experiments\terminal_weighted\supervisor_n24_rational.log"
$State = "experiments\terminal_weighted\supervisor_state_n24_rational.json"

if (Test-Path $StopFlag) { Remove-Item -Force $StopFlag }

$existing = Get-CimInstance Win32_Process -Filter "Name='python.exe'" -ErrorAction SilentlyContinue |
    Where-Object { $_.CommandLine -like '*supervise_c0008_rational_regen*' }
if ($existing) {
    Write-Host "Supervisor already running PID(s):" ($existing.ProcessId -join ', ')
    exit 0
}

$proc = Start-Process -FilePath "python" `
    -ArgumentList "scripts/supervise_c0008_rational_regen.py --n 24 --rational --workers 8 --adaptive" `
    -WorkingDirectory $PythonDir `
    -WindowStyle Hidden `
    -PassThru

Write-Host "Supervisor started PID=$($proc.Id)"
Write-Host "Log: $PythonDir\$Log"
Write-Host "State: $PythonDir\$State"
Write-Host "Stop:  .\scripts\stop_c0008_rational_supervisor.ps1"
Write-Host "Status: python scripts/regenerate_c0008_full_dealias_one_inf.py --n 24 --rational --status"
