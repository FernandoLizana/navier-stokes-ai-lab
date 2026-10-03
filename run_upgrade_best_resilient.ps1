# Resilient C-0008 Phase E upgrade — outer restart loop + anti-sleep
# Usage: .\run_upgrade_best_resilient.ps1
# Stop:  Ctrl+C (checkpoints preserved)

$ErrorActionPreference = "Continue"
Set-Location $PSScriptRoot

$logDir = "experiments/terminal_weighted"
$log = Join-Path $logDir "upgrade_best_$PID.log"
$logMain = Join-Path $logDir "upgrade_best.log"
New-Item -ItemType Directory -Force -Path $logDir | Out-Null
New-Item -ItemType Directory -Force -Path "$logDir/upgrade_shell_ck" | Out-Null

Add-Type @"
using System.Runtime.InteropServices;
public static class SleepBlock {
    [DllImport("kernel32.dll")] public static extern uint SetThreadExecutionState(uint f);
    public static void On() { SetThreadExecutionState(0x80000003); }
}
"@

function Get-ShellProgress {
    $s = python -m ns_exploration.experiments.upgrade_manifest_best_bounds --status 2>$null
    if ($s -match '^(\d+)/(\d+)$') { return [int]$Matches[1], [int]$Matches[2] }
    return 0, 87
}

function Write-LogLine([string]$Line) {
    Write-Host $Line
    Add-Content -Path $log -Value $Line -ErrorAction SilentlyContinue
    Add-Content -Path $logMain -Value $Line -ErrorAction SilentlyContinue
}

Write-Host "Phase E upgrade (double watchdog)."
Write-Host "Log sesion: $log"
Write-Host "Checkpoint: $logDir/upgrade_best_checkpoint.json"
Write-Host ""

$outer = 0
while ($true) {
    $outer++
    [SleepBlock]::On()
    $done, $total = Get-ShellProgress
    if ($done -ge $total) {
        Write-LogLine "$(Get-Date -Format 'HH:mm:ss') All $total shells done."
        break
    }
    Write-LogLine "$(Get-Date -Format 'HH:mm:ss') Launch #$outer at $done/$total shells..."

    python -m ns_exploration.experiments.upgrade_manifest_best_bounds `
        --workers 2 --prec 128 --watchdog --retry-delay 5 `
        2>&1 | ForEach-Object { Write-LogLine $_ }

    [SleepBlock]::On()
    $done, $total = Get-ShellProgress
    if ($done -ge $total) { break }

    Write-LogLine "$(Get-Date -Format 'HH:mm:ss') Process ended early ($done/$total). Outer restart in 10s..."
    Start-Sleep -Seconds 10
}

Write-Host ""
Write-Host "Finalizing certificate..."
python -m ns_exploration.experiments.finalize_phase_d_from_manifest `
    experiments/terminal_weighted/shell_manifest_best_full.json

python verify_c0008_terminal_certificate.py `
    certificates/CERT-L0072-C0008-terminal-full-dealias.json `
    experiments/terminal_weighted/shell_manifest_best_full.json

Write-Host "Done."
