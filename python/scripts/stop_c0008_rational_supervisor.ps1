# Graceful stop: supervisor + child regen workers
$ErrorActionPreference = "SilentlyContinue"
$PythonDir = Split-Path -Parent $MyInvocation.MyCommand.Path
Set-Location $PythonDir

$Base = "experiments\terminal_weighted"
$StopFlag = Join-Path $Base "supervisor_stop_n24_rational.flag"
$ChildPidFile = Join-Path $Base "supervisor_child_pid_n24_rational.txt"

New-Item -ItemType File -Path $StopFlag -Force | Out-Null

if (Test-Path $ChildPidFile) {
    $childPid = Get-Content $ChildPidFile -Raw
    if ($childPid -match '\d+') {
        taskkill /F /T /PID $matches[0] 2>$null
    }
}

Get-CimInstance Win32_Process -Filter "Name='python.exe'" |
    Where-Object {
        $_.CommandLine -like '*supervise_c0008_rational_regen*' -or
        $_.CommandLine -like '*regenerate_c0008_full_dealias_one_inf*'
    } |
    ForEach-Object { taskkill /F /T /PID $_.ProcessId 2>$null }

Write-Host "Stop flag set. Supervisor and regen processes terminated (if running)."
