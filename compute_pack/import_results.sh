#!/usr/bin/env bash
# Import results tarball from notebook into main repo compute_pack/data/
set -euo pipefail
PACK_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TAR="${1:?usage: import_results.sh results.tar.gz}"
tar -xzvf "$TAR" -C "$PACK_ROOT"
echo "Imported into $PACK_ROOT/data"
