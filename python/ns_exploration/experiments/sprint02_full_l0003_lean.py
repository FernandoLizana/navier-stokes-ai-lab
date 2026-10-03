"""
Sprint 02 — full-dealias L-0003 N5 certificate + Lean comparison-ODE skeleton.

Run from repo root (with python/ on path or editable install).
"""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

from ns_exploration.validation.l0003_certificate import (
    build_l0003_certificate,
    build_l0003_full_dealias_certificate,
    save_l0003_certificate,
)
from ns_exploration.validation.l0003_certificate_verify import verify_l0003_certificate


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
            "stdout_tail": (proc.stdout or "")[-400:],
            "lemmas": [
                "NSGalerkin.EnstrophyEnergy.enstrophy_le_K2_energy",
                "NSGalerkin.ComparisonODE.uniform_cap_closes",
                "NSGalerkin.ComparisonODE.above_eq_force_nonpos",
                "NSGalerkin.ComparisonODE.below_eq_force_nonneg",
                "NSGalerkin.Certificates.shell_k3_cap = 18450",
                "NSGalerkin.Certificates.full_N12_cap = 51450",
                "NSGalerkin.Certificates.*_omegaEq_cleared (kernel decide)",
            ],
            "clay_implication": "None",
        }
    except FileNotFoundError:
        return {"available": False, "reason": "lake not on PATH"}
    except subprocess.TimeoutExpired:
        return {"available": True, "ok": False, "reason": "lake build timeout"}


def main() -> dict:
    repo_root = Path.cwd()
    shell = build_l0003_certificate(n=12, k_shell=3)
    full = build_l0003_full_dealias_certificate(n=12)
    ok_s, ch_s = verify_l0003_certificate(shell.as_dict())
    ok_f, ch_f = verify_l0003_certificate(full.as_dict())
    save_l0003_certificate(shell, "certificates/CERT-L0003-shell-k3.json")
    save_l0003_certificate(full, "certificates/CERT-L0003-full-dealias-N12.json")

    lean = lean_build_status(repo_root)
    out = {
        "shell_k3": {
            "verified": ok_s,
            "checks": {k: bool(v) for k, v in ch_s.items()},
            "Omega_bound_hi": shell.Omega_bound_hi,
            "K2": shell.K_squared_int,
            "M": shell.M,
        },
        "full_dealias_N12": {
            "verified": ok_f,
            "checks": {k: bool(v) for k, v in ch_f.items()},
            "Omega_bound_hi": full.Omega_bound_hi,
            "K2": full.K_squared_int,
            "M": full.M,
        },
        "lean_skeleton": lean,
        "not_a_clay_claim": True,
    }
    Path("reports").mkdir(exist_ok=True)
    print(json.dumps(out, indent=2, ensure_ascii=True))
    return out


if __name__ == "__main__":
    main()
