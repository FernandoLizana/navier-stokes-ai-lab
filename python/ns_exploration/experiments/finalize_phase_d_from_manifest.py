"""Finalize Phase D certificate from an existing shell manifest (skip Frobenius pass)."""

from __future__ import annotations

import json
import sys
from pathlib import Path

from gmpy2 import mpfr

from ns_exploration.terminal_weighted.adaptive_time_certificate import (
    adaptive_lipschitz_certificate,
)
from ns_exploration.terminal_weighted.certify_shell_integral import (
    save_shell_integral_certificate,
    shell_integral_certificate,
)
from ns_exploration.terminal_weighted.certify_terminal_perturbation import (
    route_d_perturbation_certificate,
)
from ns_exploration.terminal_weighted.constants import FROZEN
from ns_exploration.terminal_weighted.intervals import i_star_interval, to_float_hi
from ns_exploration.terminal_weighted.l0048_audit import l0048_value_classification
from ns_exploration.terminal_weighted.n_coverage import n_domination_analysis
from ns_exploration.terminal_weighted.weights import stokes_floor_hi


def finalize(manifest_path: str | Path, *, prec: int = 128, cells: int = 8) -> dict:
    manifest_path = Path(manifest_path)
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    n = manifest["n"]
    mode = "full" if manifest.get("n_shells", 0) >= 87 else "band"

    route_a = shell_integral_certificate(n, prec=prec, manifest=manifest)
    save_shell_integral_certificate(
        route_a, "experiments/terminal_weighted/shell_integral_certificate.json"
    )
    route_d = route_d_perturbation_certificate(n, prec=prec, manifest=manifest)
    route_b = adaptive_lipschitz_certificate(
        n, prec=prec, manifest=manifest, initial_cells=cells
    )

    I_lo, I_hi = i_star_interval(prec=prec)
    routes = {"A": route_a, "D": route_d, "B": route_b}
    best_name = min(routes, key=lambda k: mpfr(routes[k]["integral_hi"]))
    best = routes[best_name]

    cert = {
        "claim": "C-0008",
        "status": "candidate_certificate",
        "parameters": {
            "nu": FROZEN.nu,
            "E0": FROZEN.E0,
            "T": FROZEN.T,
            "N_max": n,
            "M_target": FROZEN.c0008_target_M,
        },
        "normalization": {"fft": "fftn/N^3", "hermitian": "sprintB_physical_tensor"},
        "mode_set": {
            "n": n,
            "D": manifest["D_full"],
            "n_shells": manifest["n_shells"],
            "radii_mode": mode,
        },
        "arithmetic_backend": manifest["arithmetic_backend"],
        "stokes_floor_interval": [FROZEN.stokes_floor, to_float_hi(mpfr(stokes_floor_hi(n)))],
        "i_star_interval": [str(I_lo), str(I_hi)],
        "method": f"{best_name}_{best['method']}",
        "operator_bound_method": manifest.get("bound_method", "frobenius_mpfr_up_per_shell") + "_per_shell",
        "integral_interval": ["0", best["integral_hi"]],
        "omega_T_interval": ["0", best.get("omega_T_hi", route_a.get("omega_T_hi"))],
        "strict_margin": str(I_lo - mpfr(best["integral_hi"])),
        "closes_strict": bool(best.get("closes_strict") and mode == "full"),
        "source_hashes": {"shell_manifest": manifest.get("manifest_sha256")},
        "limitations": [
            "Triangle inequality on shell op-norms; min(Frobenius,1-inf) per shell (Phase E)",
            "N7 sketch for N<=24 domination",
        ],
        "l0048_audit": l0048_value_classification(),
        "n_coverage": n_domination_analysis(n),
        "routes": routes,
        "verification_command": (
            "python verify_c0008_terminal_certificate.py "
            "certificates/CERT-L0072-C0008-terminal-full-dealias.json "
            "experiments/terminal_weighted/shell_manifest_best_full.json"
        ),
    }

    cert_path = Path("certificates/CERT-L0072-C0008-terminal-full-dealias.json")
    cert_path.parent.mkdir(parents=True, exist_ok=True)
    cert_path.write_text(json.dumps(cert, indent=2), encoding="utf-8")

    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    from verify_c0008_terminal_certificate import verify_certificate

    vok, vmsgs = verify_certificate(cert_path, manifest_path=manifest_path)
    out = {
        "mode": mode,
        "best_route": best_name,
        "closes_strict": cert["closes_strict"],
        "integral_hi": best["integral_hi"],
        "I_star_lo": str(I_lo),
        "verifier_pass": vok,
        "verifier_msgs": vmsgs,
        "c0008_status": "unchanged_open",
    }
    Path("experiments/terminal_weighted/phase_d_run.json").write_text(
        json.dumps(out, indent=2), encoding="utf-8"
    )
    print(json.dumps(out, indent=2))
    return out


if __name__ == "__main__":
    mp = sys.argv[1] if len(sys.argv) > 1 else "experiments/terminal_weighted/shell_manifest.json"
    finalize(mp)
