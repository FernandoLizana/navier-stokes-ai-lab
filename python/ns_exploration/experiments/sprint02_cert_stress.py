"""
Sprint 02 — L-0006 interval certificate + C-S-0002 bounded stress.

Run: python -m ns_exploration.experiments.sprint02_cert_stress
(from the `python/` dir, with repo root as cwd for conjecture/certificate I/O).
"""

from __future__ import annotations

import json
from pathlib import Path

from ns_exploration.conjectures.cs0002_stress import stress_cs0002
from ns_exploration.validation.l0006_certificate import build_certificate, save_certificate
from ns_exploration.validation.l0006_certificate_verify import verify_certificate


def main() -> dict:
    cert = build_certificate(n=12, k_shell=2, E0=0.5, t=0.02)
    ok, checks = verify_certificate(cert.as_dict())
    save_certificate(cert, "certificates/CERT-L0006-shell-k2.json")

    stress = stress_cs0002()

    out = {
        "certificate": cert.as_dict(),
        "certificate_verified": ok,
        "certificate_checks": {k: bool(v) for k, v in checks.items()},
        "cs0002_stress": {
            "M_target": stress["M_target"],
            "best_gated": stress["best_gated"],
            "best_label": stress["best_label"],
            "refuted": stress["refuted"],
            "n_runs": stress["n_runs"],
        },
    }
    Path("reports").mkdir(exist_ok=True)
    print(json.dumps(out, indent=2, ensure_ascii=True))
    return out


if __name__ == "__main__":
    main()
