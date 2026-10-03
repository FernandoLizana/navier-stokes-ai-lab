#!/usr/bin/env pwsh
# Reproduce C-0008 Phase D/E certificate pipeline
#   -verify-only (default): check precomputed manifests (~seconds)
#   -full: Frobenius sprint + Phase E upgrade (~24h+) + finalize
param(
    [switch]$Full
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"
Set-Location $PSScriptRoot

$FroManifest = "experiments/terminal_weighted/shell_manifest_frobenius_full.json"
$BestManifest = "experiments/terminal_weighted/shell_manifest_best_full.json"
$ClusterManifest = "experiments/terminal_weighted/shell_manifest_cluster_best_full.json"
$Cert = "certificates/CERT-L0072-C0008-terminal-full-dealias.json"
$CertAprime = "certificates/CERT-L0073-C0008-route-a-prime-best.json"

if ($Full) {
    Write-Host "Phase D: Frobenius full-dealias sprint (~40 min)..."
    python -m ns_exploration.experiments.sprint_c0008_phase_d `
        --mode full --prec 128 --workers 4 --bound-method frobenius --cells 8

    if (-not (Test-Path $FroManifest)) {
        Write-Host "Restoring Frobenius manifest from certificate..."
        python -m ns_exploration.experiments.restore_manifest_from_cert
    }

    Write-Host "Phase E: per-shell 1-inf upgrade (~24h, checkpointed)..."
    python -m ns_exploration.experiments.upgrade_manifest_best_bounds `
        --workers 2 --prec 128 --watchdog --retry-delay 5

    Write-Host "Finalize..."
    python -m ns_exploration.experiments.finalize_phase_d_from_manifest $BestManifest
} else {
    Write-Host "Verify-only (precomputed manifests)."
    if (-not (Test-Path $BestManifest)) {
        throw "Missing $BestManifest - run with -Full or run_upgrade_best_resilient.ps1"
    }
}

Write-Host "Verifier (Phase E best manifest)..."
python verify_c0008_terminal_certificate.py $Cert $BestManifest

if (Test-Path $CertAprime) {
    Write-Host "Verifier (Route A prime cluster best)..."
    python verify_c0008_cluster_certificate.py $CertAprime $ClusterManifest
}

if (Test-Path "conjectures/active/L-0074.json") {
    Write-Host "L-0074 ladder integrity..."
    python -c "from ns_exploration.validation.l0074_certificate import verify_l0074_ladder; ok, m = verify_l0074_ladder(); print('PASS' if ok else 'FAIL', m)"
}

Write-Host "Tests..."
Push-Location python
python -m pytest ns_exploration/tests/test_phase_d.py ns_exploration/tests/test_l0073_route_a_prime.py ns_exploration/tests/test_l0074_route_a_ladder.py -q --tb=short
Pop-Location

Write-Host "Done."
