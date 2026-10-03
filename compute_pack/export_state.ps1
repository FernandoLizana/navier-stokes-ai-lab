# Export C-0008 checkpoints from dev machine into compute_pack/data/ (portable transfer)
$ErrorActionPreference = "Stop"
$PackRoot = Split-Path -Parent $MyInvocation.MyCommand.Path
$Repo = Split-Path -Parent $PackRoot
$Dest = Join-Path $PackRoot "data\terminal_weighted"
$Logs = Join-Path $PackRoot "data\logs"
New-Item -ItemType Directory -Force -Path $Dest, $Logs | Out-Null

$Sources = @(
    (Join-Path $Repo "python\experiments\terminal_weighted"),
    (Join-Path $Repo "experiments\terminal_weighted")
)

function Copy-Rel($srcRoot, $relPath) {
    foreach ($src in $Sources) {
        $f = Join-Path $src $relPath
        if (Test-Path $f) {
            $target = Join-Path $Dest $relPath
            $dir = Split-Path -Parent $target
            if ($dir) { New-Item -ItemType Directory -Force -Path $dir | Out-Null }
            if (Test-Path $f -PathType Container) {
                Copy-Item -Path $f -Destination $target -Recurse -Force
            } else {
                Copy-Item -Force $f $target
            }
            Write-Host "copied $relPath"
            return $true
        }
    }
    return $false
}

Get-ChildItem -Path $Sources[0] -ErrorAction SilentlyContinue | Where-Object {
    $_.Name -like "*rational*" -or $_.Name -like "repair_one_inf_checkpoint_n24*" -or
    $_.Name -like "repair_watchdog_state_n24*" -or $_.Name -like "repair_one_inf_heartbeat_n24*"
} | ForEach-Object {
    Copy-Rel $Sources[0] $_.Name | Out-Null
}

Copy-Rel "" "shell_one_inf_ck\n24_rational" | Out-Null
Copy-Rel "" "shell_manifest_repair_full_one_inf_n24.json" | Out-Null

$manifest = @{
    exported_at = (Get-Date).ToUniversalTime().ToString("o")
    source_repo = $Repo
    dest = $Dest
    note = "Copy entire compute_pack/ folder to Ubuntu notebook"
} | ConvertTo-Json
Set-Content -Path (Join-Path $PackRoot "data\export_manifest.json") -Value $manifest -Encoding UTF8
Write-Host "Done. Transfer: compress compute_pack/ + python/ (or full repo) to Ubuntu."
