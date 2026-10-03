# Route A full-dealias ladder - clean run + watchdog retries until complete.
param(
    [int]$Workers = 2,
    [int]$MaxAttempts = 5,
    [switch]$Clean,
    [switch]$FromScratch,
    [switch]$UseSeed,
    [switch]$NoFinalize
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"
Set-Location $PSScriptRoot

$Log = "experiments/terminal_weighted/route_a_ladder_full.log"
$LogFallback = "experiments/terminal_weighted/route_a_ladder_full_$(Get-Date -Format 'yyyyMMdd_HHmmss').log"
$Checkpoint = "experiments/terminal_weighted/route_a_ladder_full_checkpoint.json"
$Summary = "experiments/terminal_weighted/route_a_ladder_full.json"
$Labels = @("equal_12", "equal_24", "equal_48", "fine_low")
$script:ActiveLog = $Log

function Write-Log {
    param([string]$Msg)
    $line = "$(Get-Date -Format 'yyyy-MM-dd HH:mm:ss') $Msg"
    Write-Host $line
    try {
        Add-Content -Path $script:ActiveLog -Value $line -ErrorAction Stop
    } catch {
        if ($script:ActiveLog -eq $Log) {
            $script:ActiveLog = $LogFallback
            Write-Host "Log locked; using $script:ActiveLog"
            Add-Content -Path $script:ActiveLog -Value $line
        }
    }
}

function Get-DoneLabels {
    if (-not (Test-Path $Checkpoint)) { return @() }
    $ck = Get-Content $Checkpoint -Raw | ConvertFrom-Json
    return @($ck.rows | Where-Object { -not $_.error } | ForEach-Object { $_.label })
}

function Invoke-Config {
    param([string]$Label)
    $pyArgs = @(
        "-m", "ns_exploration.experiments.sprint_c0008_route_a_ladder",
        "--full-sweep",
        "--workers", "$Workers",
        "--resume",
        "--only", $Label
    )
    if ($UseSeed -and $Label -eq "equal_12") {
        $pyArgs += "--seed-l0073"
    }
    Write-Log "START $Label (workers=$Workers)"
    & python @pyArgs 2>&1 | ForEach-Object { Write-Log $_ }
    return $LASTEXITCODE
}

if ($Clean -or $FromScratch) {
    Remove-Item $Checkpoint -ErrorAction SilentlyContinue
    Remove-Item $Summary -ErrorAction SilentlyContinue
    Remove-Item experiments/terminal_weighted/route_a_manifest_equal_*.json -ErrorAction SilentlyContinue
    Remove-Item experiments/terminal_weighted/route_a_manifest_fine_low.json -ErrorAction SilentlyContinue
    Remove-Item experiments/terminal_weighted/shell_manifest_route_a_best_full.json -ErrorAction SilentlyContinue
    if ($Clean) {
        Remove-Item $Log -ErrorAction SilentlyContinue
    }
    Write-Log "Cleaned ladder artifacts."
}

Write-Log "Route A ladder watchdog - labels: $($Labels -join ', ')"

$attempt = 0
while ($true) {
    $done = Get-DoneLabels
    $pending = @($Labels | Where-Object { $_ -notin $done })
    if ($pending.Count -eq 0) {
        Write-Log "All configs complete."
        break
    }
    if ($attempt -ge $MaxAttempts) {
        Write-Log "ABORT: max attempts ($MaxAttempts) reached. Pending: $($pending -join ', ')"
        exit 1
    }
    $attempt++
    Write-Log "Watchdog attempt $attempt - pending: $($pending -join ', ')"
    foreach ($label in $pending) {
        $rc = Invoke-Config -Label $label
        if ($rc -ne 0) {
            Write-Log "WARN $label exited with code $rc - will retry on next pass"
        }
        Start-Sleep -Seconds 3
    }
}

Write-Log "Finalizing summary..."
$finalizeArgs = @(
    "-m", "ns_exploration.experiments.sprint_c0008_route_a_ladder",
    "--full-sweep", "--workers", "$Workers", "--resume"
)
if ($UseSeed) { $finalizeArgs += "--seed-l0073" }
& python @finalizeArgs 2>&1 | ForEach-Object { Write-Log $_ }

if (-not (Test-Path $Summary)) {
    Write-Log "ERROR: $Summary not created"
    exit 1
}

if (-not $NoFinalize) {
    Write-Log "Updating certificate + verify..."
    python -m ns_exploration.experiments.finalize_route_a_ladder
    powershell -NoProfile -ExecutionPolicy Bypass -File .\reproduce_c0008_certificate.ps1
}

Write-Log "DONE."
