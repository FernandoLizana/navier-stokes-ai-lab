#!/usr/bin/env bash
set -euo pipefail
PACK_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=env.sh
source "$PACK_ROOT/env.sh"

cd "$BASTARDUS_ROOT/python"
python3 -m venv "$PACK_ROOT/.venv"
# shellcheck disable=SC1091
source "$PACK_ROOT/.venv/bin/activate"
pip install -U pip wheel
pip install -r "$PACK_ROOT/requirements.txt"
pip install -r "$BASTARDUS_ROOT/python/requirements-c0008-repair.txt" 2>/dev/null || true
echo "OK: venv at $PACK_ROOT/.venv"
echo "Next: source $PACK_ROOT/env.sh && $PACK_ROOT/run.sh list"
