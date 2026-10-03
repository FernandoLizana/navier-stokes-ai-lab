"""
Sprint 02 — N16 L-0003 N5 certificate + Lean ℚ enstrophy–energy lift.

Documents ∼M/ν² growth of the uniform Galerkin cap and the Rat-form structural
lemma behind Ω ≤ K² E.
"""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

from ns_exploration.validation.l0003_certificate import (
    build_l0003_full_dealias_certificate,
    save_l0003_certificate,
)
from ns_exploration.validation.l0003_certificate_verify import verify_l0003_certificate


def lean_build_status(repo_root: Path) -> dict:
    lean_dir = repo_root / "lean" / "NSGalerkin"
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
            "ok": proc.returncode == 0,
            "exit_code": proc.returncode,
            "stdout_tail": (proc.stdout or "")[-400:],
            "lemmas": [
                "NSGalerkin.EnstrophyEnergy.enstrophy_le_K2_energy",
                "NSGalerkin.EnstrophyEnergyRat.enstrophy_le_K2_energy_rat",
                "NSGalerkin.Certificates.full_N16_cap = 199650",
                "NSGalerkin.Certificates.full_growth_N12_to_N16",
            ],
            "clay_implication": "None",
        }
    except FileNotFoundError:
        return {"available": False, "reason": "lake not on PATH"}


def main() -> dict:
    c12 = build_l0003_full_dealias_certificate(n=12)
    c16 = build_l0003_full_dealias_certificate(n=16)
    ok12, ch12 = verify_l0003_certificate(c12.as_dict())
    ok16, ch16 = verify_l0003_certificate(c16.as_dict())
    save_l0003_certificate(c16, "certificates/CERT-L0003-full-dealias-N16.json")

    lean = lean_build_status(Path.cwd())
    out = {
        "N12": {
            "verified": ok12,
            "K2": c12.K_squared_int,
            "M": c12.M,
            "Omega_bound_hi": c12.Omega_bound_hi,
        },
        "N16": {
            "verified": ok16,
            "K2": c16.K_squared_int,
            "M": c16.M,
            "Omega_bound_hi": c16.Omega_bound_hi,
            "checks": {k: bool(v) for k, v in ch16.items()},
        },
        "growth": {
            "M_ratio": c16.M / c12.M,
            "cap_ratio": c16.Omega_bound_hi / c12.Omega_bound_hi,
            "note": "cap_ratio == M_ratio for fixed (E0,nu); documents ~M/nu^2 blow-up",
        },
        "lean_skeleton": lean,
        "not_a_clay_claim": True,
    }
    Path("reports").mkdir(exist_ok=True)
    print(json.dumps(out, indent=2, ensure_ascii=True))
    return out


if __name__ == "__main__":
    main()
