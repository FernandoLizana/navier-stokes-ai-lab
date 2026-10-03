# C-0008 full-dealias repair — max CPU with adaptive RAM/OOM resilience
param(
    [int]$N = 24,
    [int]$Prec = 128,
    [int]$Workers = 8,
    [int]$MinWorkers = 1
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"
$Root = Split-Path -Parent $MyInvocation.MyCommand.Path
$Py = Join-Path $Root "python"
Set-Location $Py

# Prevent sleep during long compute
powercfg /change standby-timeout-ac 0 | Out-Null
powercfg /change hibernate-timeout-ac 0 | Out-Null

Write-Host "== C-0008 repair full-dealias n=$N workers=$Workers adaptive =="
Write-Host "Checkpoints: repair_one_inf_checkpoint_n$N.json + shell_one_inf_ck/n$N/"
Write-Host "Watchdog state: repair_watchdog_state_n$N.json"
Write-Host ""

python scripts/regenerate_c0008_full_dealias_one_inf.py `
    --n $N `
    --prec $Prec `
    --workers $Workers `
    --min-workers $MinWorkers `
    --adaptive `
    --watchdog `
    --retry-delay 15

Write-Host "== finalize on completion =="
python scripts/finalize_c0008_repair.py --n $N --prec $Prec --workers 1

Write-Host "DONE resilient full-dealias"
