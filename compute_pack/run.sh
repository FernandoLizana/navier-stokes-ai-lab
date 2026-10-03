#!/usr/bin/env bash
set -euo pipefail
PACK_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=env.sh
source "$PACK_ROOT/env.sh"

if [[ -f "$PACK_ROOT/.venv/bin/activate" ]]; then
  # shellcheck disable=SC1091
  source "$PACK_ROOT/.venv/bin/activate"
fi

cd "$BASTARDUS_ROOT/python"

_recommended_workers() {
  python3 "$PACK_ROOT/ram_workers.py" --max "${C0008_WORKERS_MAX:-32}"
}

case "${1:-}" in
  list)
    python3 "$PACK_ROOT/supervisor.py" --list
    ;;
  status|status-all)
    python3 "$PACK_ROOT/supervisor.py" --status-all
    ;;
  stop)
    python3 "$PACK_ROOT/supervisor.py" --stop "${2:-}"
    ;;
  queue)
    shift || true
    W="${1:-$(_recommended_workers)}"
    echo "queue workers=$W (RAM auto; override: ./run.sh queue 16)"
    python3 "$PACK_ROOT/supervisor.py" --queue --workers "$W" --adaptive
    ;;
  ram)
    python3 "$PACK_ROOT/ram_workers.py" --json
    ;;
  "")
    echo "Usage: $0 {list|status|ram|stop [job_id]|queue [workers]|JOB_ID [workers]}"
    echo "  workers omitted → auto from RAM (up to C0008_WORKERS_MAX, default 32)"
    exit 1
    ;;
  *)
    JOB="$1"
    W="${2:-$(_recommended_workers)}"
    echo "workers=$W ($(python3 "$PACK_ROOT/ram_workers.py" --json))"
    nohup python3 "$PACK_ROOT/supervisor.py" "$JOB" --workers "$W" --adaptive \
      >> "$C0008_LOG_DIR/nohup_${JOB}.log" 2>&1 &
    echo "Started job=$JOB workers=$W pid=$!"
    echo "  tail -f $C0008_LOG_DIR/supervisor_${JOB}.log"
    echo "  $0 status"
    ;;
esac
