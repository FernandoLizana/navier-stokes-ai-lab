# Full C-0008 repair orchestration (Windows)
param(
    [int]$N = 24,
    [int]$Prec = 128,
    [int]$Workers = 1,
    [switch]$Watchdog,
    [switch]$FinalizeOnly,
    [switch]$Bootstrap
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"
$Root = Split-Path -Parent $MyInvocation.MyCommand.Path
$Py = Join-Path $Root "python"
Set-Location $Py

if ($Bootstrap) {
    & (Join-Path $Root "scripts\bootstrap_c0008_repair.ps1")
}

if (-not $FinalizeOnly) {
    Write-Host "== full-dealias one_inf regeneration n=$N =="
    $regArgs = @("scripts/regenerate_c0008_full_dealias_one_inf.py", "--n", "$N", "--prec", "$Prec", "--workers", "$Workers")
    if ($Watchdog) { $regArgs += "--watchdog" }
    python @regArgs
}

Write-Host "== finalize repair pipeline =="
python scripts/finalize_c0008_repair.py --n $N --prec $Prec --workers $Workers

if ($Bootstrap) {
    Write-Host "== verify band cert =="
    python ../verify_c0008_terminal_certificate.py ../certificates/CERT-C0008-repair-band6-one-inf.json ../experiments/terminal_weighted/shell_manifest_repair_band6_one_inf.json
}

Write-Host "DONE run_c0008_repair_all"
