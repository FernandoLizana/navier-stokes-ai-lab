# C-0008 repair — bootstrap and band gate tests
Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"
$Root = Split-Path -Parent $MyInvocation.MyCommand.Path
$Py = Join-Path $Root "python"
Set-Location $Py

Write-Host "== install repair deps =="
python -m pip install -q -r requirements-c0008-repair.txt

Write-Host "== band repair tests =="
python -m pytest `
  ns_exploration/tests/test_c0008_cert_repair.py `
  ns_exploration/tests/test_exact_prototype_band123.py `
  ns_exploration/tests/test_cluster_l0073r_manifest.py `
  ns_exploration/tests/test_repair_full_manifest.py `
  -q --tb=line

Write-Host "== finalize band certificate =="
python scripts/finalize_c0008_repair.py --band --prec 128

Write-Host "== band status =="
python scripts/regenerate_c0008_full_dealias_one_inf.py --n 24 --status

Write-Host "DONE bootstrap"
