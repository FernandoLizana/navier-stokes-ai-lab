#!/usr/bin/env python3
"""Independent adversarial verifier for C-0008 terminal certificate."""

from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT / "python") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "python"))

from ns_exploration.validation.c0008_verify import verify_terminal_file  # noqa: E402


def main(argv: list[str] | None = None) -> int:
    argv = argv or sys.argv[1:]
    if not argv:
        print("Usage: verify_c0008_terminal_certificate.py <certificate.json> [manifest.json]")
        return 2
    cert_path = Path(argv[0])
    if not cert_path.is_absolute():
        cert_path = REPO_ROOT / cert_path
    manifest_path = Path(argv[1]) if len(argv) > 1 else None
    if manifest_path and not manifest_path.is_absolute():
        manifest_path = REPO_ROOT / manifest_path
    ok, verdict = verify_terminal_file(cert_path, manifest_path=manifest_path)
    print(verdict)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
