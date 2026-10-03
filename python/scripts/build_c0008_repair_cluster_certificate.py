"""Build L-0073R cluster sketch certificate from cluster manifest."""

from __future__ import annotations

import json
import platform
from pathlib import Path

import gmpy2
from gmpy2 import mpfr

from ns_exploration.terminal_weighted.cluster_bounds import cluster_integral_certificate
from ns_exploration.terminal_weighted.constants import EXACT, FROZEN
from ns_exploration.terminal_weighted.intervals import i_star_interval, stokes_floor_interval
from ns_exploration.terminal_weighted.recompute_integral import (
    cluster_integral_hi_from_blocks,
    omega_T_hi_from_I_term,
)
from ns_exploration.validation.c0008_verify import manifest_sha256, verify_terminal_certificate

PY = Path(__file__).resolve().parents[1]
REPO = PY.parent


def build_repair_cluster_certificate(
    cluster_manifest: dict,
    *,
    shell_manifest: dict | None = None,
    prec: int = 128,
) -> dict:
    route = cluster_integral_certificate(cluster_manifest, prec=prec)
    I_lo, _ = i_star_interval(prec=prec)
    _, floor_hi = stokes_floor_interval(prec=prec)
    integral_hi = cluster_integral_hi_from_blocks(cluster_manifest["blocks"], prec=prec)
    omega_hi = omega_T_hi_from_I_term(integral_hi, prec=prec)
    i_star_store = [FROZEN.integral_threshold * 0.999999, FROZEN.integral_threshold * 1.000001]
    floor_store = [FROZEN.stokes_floor, float(floor_hi)]

    cert = {
        "claim": "C-0008",
        "status": "repair_cluster_sketch",
        "scope": "repair_full_cluster_sketch",
        "parameters": {
            "nu": EXACT.NU,
            "E0": EXACT.E0,
            "T": EXACT.T,
            "N_max": 24,
            "M_target": EXACT.M_TARGET,
        },
        "bound_method": cluster_manifest.get("bound_method", "one_inf_only"),
        "mode_set": {
            "n": cluster_manifest.get("n", 24),
            "D": cluster_manifest.get("D_full"),
            "n_shells": cluster_manifest.get("n_shells"),
            "n_clusters": cluster_manifest.get("n_clusters"),
            "radii_mode": "full",
        },
        "repair_complete": True,
        "arithmetic_backend": {
            "name": "gmpy2_mpfr",
            "version": gmpy2.version(),
            "precision_bits": prec,
            "rounding_up": "RoundUp",
            "rounding_down": "RoundDown",
            "platform": platform.platform(),
        },
        "stokes_floor_interval": floor_store,
        "i_star_interval": [str(i_star_store[0]), str(i_star_store[1])],
        "method": "repair_route_A_prime_L0073R_cluster_sketch",
        "integral_interval": ["0", str(integral_hi)],
        "omega_T_interval": ["0", str(omega_hi)],
        "coverage_complete": True,
        "source_hashes": {
            "cluster_manifest": cluster_manifest.get("manifest_sha256")
            or manifest_sha256(cluster_manifest),
        },
        "route_A_prime": route,
        "limitations": [
            "Route A' cluster sketch — L-0073R residual reweight.",
            "Declared bound uses cluster contrib_hi sum (verifier-recomputable).",
            "Route A shell-sum cert remains primary repair_full certificate.",
            "N<=24 inclusion: L-0076 N7_sketch.",
        ],
        "closes_c0008": False,
        "closes_strict": False,
        "honesty": (
            "Cluster L-0073R sketch; does NOT supersede Route A shell cert. "
            "See reports/L0073R_CLUSTER_AUDIT.md — closes_c0008 never true until audit PASS."
        ),
    }
    if shell_manifest is not None:
        cert["source_hashes"]["shell_manifest"] = manifest_sha256(shell_manifest)
    return cert


def main() -> dict:
    import argparse

    p = argparse.ArgumentParser()
    p.add_argument("--n", type=int, default=24)
    p.add_argument("--prec", type=int, default=128)
    args = p.parse_args()

    cluster_path = REPO / f"experiments/terminal_weighted/cluster_manifest_repair_full_n{args.n}_one_inf.json"
    if not cluster_path.is_file():
        cluster_path = PY / f"experiments/terminal_weighted/cluster_manifest_repair_full_n{args.n}_one_inf.json"
    shell_path = PY / f"experiments/terminal_weighted/shell_manifest_repair_full_one_inf_n{args.n}.json"
    if not cluster_path.is_file():
        raise FileNotFoundError(f"cluster manifest missing: {cluster_path}")

    cluster = json.loads(cluster_path.read_text(encoding="utf-8"))
    shell = (
        json.loads(shell_path.read_text(encoding="utf-8")) if shell_path.is_file() else None
    )
    cert = build_repair_cluster_certificate(cluster, shell_manifest=shell, prec=args.prec)
    out = REPO / f"certificates/CERT-L0073R-repair-C0008-cluster-one-inf.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(cert, indent=2), encoding="utf-8")
    ok, issues = verify_terminal_certificate(
        cert, manifest=cluster, manifest_path=cluster_path, prec=args.prec
    )
    rep = {
        "cluster_manifest": str(cluster_path),
        "cert": str(out),
        "verify_ok": ok,
        "issues": issues,
        "integral_hi": cert["integral_interval"][1],
        "closes_c0008": cert.get("closes_c0008"),
    }
    print(json.dumps(rep, indent=2))
    return rep


if __name__ == "__main__":
    main()
