#!/usr/bin/env python3
"""Root entry — independent C-0008 terminal verifier."""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "python"))

from ns_exploration.validation.c0008_verify import (  # noqa: E402
    verify_terminal_certificate,
    verify_terminal_file,
)


def verify_certificate(
    cert_path: Path | str,
    manifest_path: Path | str | None = None,
) -> tuple[bool, list[str]]:
    """Legacy API: (ok, messages) with FAIL-prefixed issue strings."""
    cert_path = Path(cert_path)
    manifest_path = Path(manifest_path) if manifest_path else None
    cert = json.loads(cert_path.read_text(encoding="utf-8"))
    prec = int(cert.get("arithmetic_backend", {}).get("precision_bits", 200))
    manifest: dict | None = None
    if manifest_path and manifest_path.is_file():
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    ok, issues = verify_terminal_certificate(
        cert, manifest=manifest, manifest_path=manifest_path, prec=prec
    )
    if ok:
        return True, ["PASS"]
    return False, [f"FAIL: {i}" if not str(i).startswith("FAIL") else str(i) for i in issues]


def main(argv: list[str] | None = None) -> int:
    argv = argv or sys.argv[1:]
    if not argv:
        print("Usage: verify_c0008_terminal_certificate.py <certificate.json> [manifest.json]")
        return 2
    cert_path = Path(argv[0])
    if not cert_path.is_absolute():
        cert_path = ROOT / cert_path
    manifest_path = Path(argv[1]) if len(argv) > 1 else None
    if manifest_path and not manifest_path.is_absolute():
        manifest_path = ROOT / manifest_path
    ok, verdict = verify_terminal_file(cert_path, manifest_path=manifest_path)
    print(verdict)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
