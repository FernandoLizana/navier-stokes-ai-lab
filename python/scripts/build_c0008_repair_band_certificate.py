"""Build and verify C-0008 repair-band prototype certificate."""

from __future__ import annotations

import json
from pathlib import Path

from ns_exploration.validation.c0008_repair_band_cert import save_repair_band_certificate
from ns_exploration.validation.c0008_verify import verify_terminal_certificate

REPO = Path(__file__).resolve().parents[2]
MANIFEST = REPO / "experiments" / "terminal_weighted" / "shell_manifest_repair_band6_one_inf.json"
CERT_OUT = REPO / "certificates" / "CERT-C0008-repair-band6-one-inf.json"


def main() -> dict:
    cert = save_repair_band_certificate(MANIFEST, CERT_OUT, prec=128)
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    ok, issues = verify_terminal_certificate(
        cert, manifest=manifest, manifest_path=MANIFEST, prec=128
    )
    print({"cert": str(CERT_OUT), "verify_ok": ok, "issues": issues})
    return {"cert_path": str(CERT_OUT), "verify_ok": ok, "issues": issues}


if __name__ == "__main__":
    main()
