#!/usr/bin/env bash
# After compute on notebook: tarball for transfer back
set -euo pipefail
PACK_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT="${1:-c0008_compute_results_$(date +%Y%m%d_%H%M).tar.gz}"
tar -czvf "$OUT" -C "$PACK_ROOT" data
echo "Created $OUT — copy to main machine and extract into compute_pack/"
