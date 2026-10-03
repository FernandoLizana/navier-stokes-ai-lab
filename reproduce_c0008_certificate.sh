#!/usr/bin/env bash
# Reproduce C-0008 Phase D/E certificate pipeline
# Usage: ./reproduce_c0008_certificate.sh          # verify-only
#        ./reproduce_c0008_certificate.sh --full   # full recompute (~24h+)
set -euo pipefail
cd "$(dirname "$0")"

FULL=0
[[ "${1:-}" == "--full" ]] && FULL=1

FRO="experiments/terminal_weighted/shell_manifest_frobenius_full.json"
BEST="experiments/terminal_weighted/shell_manifest_best_full.json"
CERT="certificates/CERT-L0072-C0008-terminal-full-dealias.json"

if [[ "$FULL" -eq 1 ]]; then
  echo "Phase D: Frobenius sprint..."
  python -m ns_exploration.experiments.sprint_c0008_phase_d \
    --mode full --prec 128 --workers 4 --bound-method frobenius --cells 8
  [[ -f "$FRO" ]] || python -m ns_exploration.experiments.restore_manifest_from_cert
  echo "Phase E: upgrade (watchdog)..."
  python -m ns_exploration.experiments.upgrade_manifest_best_bounds \
    --workers 2 --prec 128 --watchdog --retry-delay 5
  python -m ns_exploration.experiments.finalize_phase_d_from_manifest "$BEST"
else
  echo "Verify-only."
  [[ -f "$BEST" ]] || { echo "Missing $BEST"; exit 1; }
fi

python verify_c0008_terminal_certificate.py "$CERT" "$BEST"
(cd python && python -m pytest ns_exploration/tests/test_phase_d.py -q --tb=short)
echo "Done."
