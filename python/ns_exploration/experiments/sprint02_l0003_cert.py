"""
Sprint 02 — L-0003 N5 certificate (k<=3 fallback) + Lean skeleton status.

Run from repo root:
  python -c "import sys; sys.path.insert(0,'python'); ..."
or via Makefile target sprint2-l0003-cert.
"""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

from ns_exploration.validation.l0003_certificate import (
    build_l0003_certificate,
    save_l0003_certificate,
)
from ns_exploration.validation.l0003_certificate_verify import verify_l0003_certificate
from ns_exploration.validation.l0006_certificate import build_certificate


def lean_build_status(repo_root: Path) -> dict:
    lean_dir = repo_root / "lean" / "NSGalerkin"
    if not lean_dir.is_dir():
        return {"available": False, "reason": "lean/NSGalerkin missing"}
    try:
        proc = subprocess.run(
            ["lake", "build"],
            cwd=str(lean_dir),
            capture_output=True,
            text=True,
            timeout=180,
        )
        return {
            "available": True,
            "exit_code": proc.returncode,
            "ok": proc.returncode == 0,
            "stdout_tail": (proc.stdout or "")[-500:],
            "stderr_tail": (proc.stderr or "")[-500:],
            "lemma": "NSGalerkin.EnstrophyEnergy.enstrophy_le_K2_energy",
            "clay_implication": "None",
        }
    except FileNotFoundError:
        return {"available": False, "reason": "lake not on PATH"}
    except subprocess.TimeoutExpired:
        return {"available": True, "ok": False, "reason": "lake build timeout"}


def main() -> dict:
    repo_root = Path.cwd()
    # Document L-0006 failure at k=3
    fail6 = build_certificate(n=12, k_shell=3, E0=0.5, t=0.02)
    cert3 = build_l0003_certificate(n=12, k_shell=3, E0=0.5, nu=0.1)
    ok, checks = verify_l0003_certificate(cert3.as_dict())
    save_l0003_certificate(cert3, "certificates/CERT-L0003-shell-k3.json")

    lean = lean_build_status(repo_root)

    out = {
        "l0006_k3_algebraic": {
            "closes": fail6.closes,
            "prod_hi": fail6.prod_hi,
            "Omega_bound_hi": fail6.Omega_bound_hi,
        },
        "l0003_certificate": cert3.as_dict(),
        "l0003_verified": ok,
        "l0003_checks": {k: bool(v) for k, v in checks.items()},
        "lean_skeleton": lean,
        "not_a_clay_claim": True,
    }
    Path("reports").mkdir(exist_ok=True)
    print(json.dumps(out, indent=2, ensure_ascii=True))
    return out


if __name__ == "__main__":
    main()
