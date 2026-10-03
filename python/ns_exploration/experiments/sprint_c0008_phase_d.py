"""Phase D certification sprint."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from gmpy2 import mpfr

from ns_exploration.conjectures.l0048_matrixfree_fullsym import all_dealias_radii
from ns_exploration.terminal_weighted.adaptive_time_certificate import (
    adaptive_lipschitz_certificate,
)
from ns_exploration.terminal_weighted.certify_shell_integral import (
    shell_integral_certificate,
    save_shell_integral_certificate,
)
from ns_exploration.terminal_weighted.certify_terminal_perturbation import (
    route_d_perturbation_certificate,
)
from ns_exploration.terminal_weighted.constants import FROZEN
from ns_exploration.terminal_weighted.intervals import i_star_interval, to_float_hi
from ns_exploration.terminal_weighted.l0048_audit import l0048_value_classification
from ns_exploration.terminal_weighted.n_coverage import n_domination_analysis
from ns_exploration.terminal_weighted.shell_decomposition import (
    build_shell_manifest,
    save_shell_manifest,
)
from ns_exploration.terminal_weighted.weights import stokes_floor_hi


def _pick_best_route(routes: dict[str, dict]) -> tuple[str, dict]:
    return min(routes.items(), key=lambda kv: mpfr(kv[1]["integral_hi"]))


def main(argv: list[str] | None = None) -> dict:
    config_path = Path("experiments/terminal_weighted/phase_d_config.json")
    cfg = json.loads(config_path.read_text(encoding="utf-8")) if config_path.is_file() else {}

    parser = argparse.ArgumentParser(description="C-0008 Phase D certification sprint")
    parser.add_argument(
        "--mode",
        choices=("band", "full"),
        default=cfg.get("mode", "band"),
        help="band=radii 1..6 pilot; full=all dealias shells",
    )
    parser.add_argument("--prec", type=int, default=cfg.get("precision_bits", 200))
    parser.add_argument("--cells", type=int, default=cfg.get("initial_cells", 4))
    parser.add_argument("--workers", type=int, default=cfg.get("workers", 1))
    parser.add_argument("--bound-method", choices=("frobenius", "one_inf", "best"), default=cfg.get("bound_method", "best"))
    args = parser.parse_args(argv)

    n = 24
    radii = tuple(range(1, 7)) if args.mode == "band" else all_dealias_radii(n)

    print(f"Phase D mode={args.mode} radii={len(radii)}", flush=True)
    print("L-0048 audit...", flush=True)
    l48 = l0048_value_classification()
    n_cov = n_domination_analysis(n)

    print("Building shell manifest (MPFR Frobenius per shell)...", flush=True)
    manifest = build_shell_manifest(
        n, prec=args.prec, radii=radii, workers=args.workers, bound_method=args.bound_method
    )
    save_shell_manifest(manifest, "experiments/terminal_weighted/shell_manifest.json")

    print("Route A: shell integral...", flush=True)
    route_a = shell_integral_certificate(n, prec=args.prec, manifest=manifest)
    save_shell_integral_certificate(
        route_a, "experiments/terminal_weighted/shell_integral_certificate.json"
    )

    print("Route D: terminal perturbation...", flush=True)
    route_d = route_d_perturbation_certificate(n, prec=args.prec, manifest=manifest)

    print(f"Route B: adaptive Lipschitz ({args.cells} cells)...", flush=True)
    route_b = adaptive_lipschitz_certificate(
        n, prec=args.prec, manifest=manifest, initial_cells=args.cells
    )

    I_lo, I_hi = i_star_interval(prec=args.prec)
    routes = {"A": route_a, "D": route_d, "B": route_b}
    best_name, best = _pick_best_route(routes)

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
            "radii_mode": args.mode,
            "radii_count": len(radii),
        },
        "arithmetic_backend": manifest["arithmetic_backend"],
        "stokes_floor_interval": [FROZEN.stokes_floor, to_float_hi(mpfr(stokes_floor_hi(n)))],
        "i_star_interval": [str(I_lo), str(I_hi)],
        "method": f"{best_name}_{best['method']}",
        "operator_bound_method": "frobenius_mpfr_up_per_shell",
        "integral_interval": ["0", best["integral_hi"]],
        "omega_T_interval": ["0", best.get("omega_T_hi", route_a.get("omega_T_hi", "unknown"))],
        "strict_margin": str(I_lo - mpfr(best["integral_hi"])),
        "closes_strict": bool(best.get("closes_strict", False) and args.mode == "full"),
        "source_hashes": {"shell_manifest": manifest.get("manifest_sha256")},
        "limitations": [
            "Frobenius per-shell bounds (triangle inequality on shells)",
            "N7 sketch for N<=24 domination — not formal embedding proof",
            f"L-0048 classification: {l48['classification']}",
            "band mode does NOT close all-IC" if args.mode == "band" else "full-dealias shell set",
        ],
        "l0048_audit": l48,
        "n_coverage": n_cov,
        "routes": routes,
        "verification_command": (
            "python verify_c0008_terminal_certificate.py "
            "certificates/CERT-L0072-C0008-terminal-full-dealias.json"
        ),
    }

    cert_path = Path("certificates/CERT-L0072-C0008-terminal-full-dealias.json")
    cert_path.parent.mkdir(parents=True, exist_ok=True)
    cert_path.write_text(json.dumps(cert, indent=2), encoding="utf-8")

    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    from verify_c0008_terminal_certificate import verify_certificate

    vok, vmsgs = verify_certificate(cert_path)

    out = {
        "mode": args.mode,
        "best_route": best_name,
        "closes_strict": best.get("closes_strict"),
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
    main()
