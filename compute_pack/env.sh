#!/usr/bin/env bash
# Source before any compute: sets portable data dirs under compute_pack/data/
PACK_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export BASTARDUS_ROOT="$(cd "$PACK_ROOT/.." && pwd)"
export C0008_DATA_DIR="$PACK_ROOT/data/terminal_weighted"
export C0008_LOG_DIR="$PACK_ROOT/data/logs"
export C0008_CERT_DIR="$PACK_ROOT/data/certificates"
export PYTHONPATH="$BASTARDUS_ROOT/python"
export C0008_WORKERS_MAX="${C0008_WORKERS_MAX:-32}"
mkdir -p "$C0008_DATA_DIR" "$C0008_LOG_DIR" "$C0008_CERT_DIR"
