"""Portable data directory for C-0008 compute (copy-friendly single folder)."""

from __future__ import annotations

import os
from pathlib import Path


def data_root() -> Path:
    """All checkpoints, manifests, heartbeats under this tree when set."""
    env = os.environ.get("C0008_DATA_DIR")
    if env:
        p = Path(env)
    else:
        p = Path("experiments/terminal_weighted")
    p.mkdir(parents=True, exist_ok=True)
    return p


def log_root() -> Path:
    env = os.environ.get("C0008_LOG_DIR")
    if env:
        p = Path(env)
    else:
        p = data_root() / "logs"
    p.mkdir(parents=True, exist_ok=True)
    return p


def cert_root() -> Path:
    env = os.environ.get("C0008_CERT_DIR")
    if env:
        p = Path(env)
    else:
        p = Path("certificates")
    p.mkdir(parents=True, exist_ok=True)
    return p


def rational_tag(rational: bool) -> str:
    return "_rational" if rational else ""
