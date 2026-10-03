"""Finalize Route A ladder: update L-0073 cert if ladder improved bound."""

from __future__ import annotations

import json
import sys
from pathlib import Path

from gmpy2 import mpfr

ROOT = Path(__file__).resolve().parents[3]
LADDER_FULL = ROOT / "experiments/terminal_weighted/route_a_ladder_full.json"
L0073_MANIFEST = ROOT / "experiments/terminal_weighted/shell_manifest_cluster_best_full.json"


def finalize(*, prec: int = 128) -> dict:
    if not LADDER_FULL.is_file():
        raise FileNotFoundError(f"missing {LADDER_FULL}")

    ladder = json.loads(LADDER_FULL.read_text(encoding="utf-8"))
    best = ladder["best"]
    best_manifest = ROOT / ladder["best_manifest"]
    if not best_manifest.is_file():
        raise FileNotFoundError(f"missing best manifest {best_manifest}")

    l0073_I = None
    if L0073_MANIFEST.is_file():
        from ns_exploration.terminal_weighted.cluster_bounds import cluster_integral_certificate

        old = json.loads(L0073_MANIFEST.read_text(encoding="utf-8"))
        l0073_I = mpfr(cluster_integral_certificate(old, prec=prec)["integral_hi"])

    new_I = mpfr(best["best_I_hi"])
    improved = l0073_I is None or new_I < l0073_I

    out = {
        "ladder_best_label": best["label"],
        "ladder_best_I_hi": best["best_I_hi"],
        "ladder_best_omega": best["omega_T_hi"],
        "l0073_I_hi": str(l0073_I) if l0073_I is not None else None,
        "improved_vs_l0073": improved,
        "best_manifest": str(best_manifest.relative_to(ROOT)).replace("\\", "/"),
    }

    if improved and best["label"] != "equal_12":
        # Copy new best to L-0073 manifest path and re-finalize cert
        L0073_MANIFEST.write_text(best_manifest.read_text(encoding="utf-8"), encoding="utf-8")
        from ns_exploration.experiments.finalize_route_a_prime_from_manifest import finalize as fin_aprime

        cert_out = fin_aprime(L0073_MANIFEST, prec=prec)
        out["certificate_updated"] = True
        out["certificate"] = cert_out
    else:
        from ns_exploration.experiments.finalize_route_a_prime_from_manifest import finalize as fin_aprime

        cert_out = fin_aprime(L0073_MANIFEST, prec=prec)
        out["certificate_updated"] = False
        out["certificate"] = cert_out
        out["note"] = "L-0073 remains best (equal_12 or no improvement)"

    run_path = ROOT / "experiments/terminal_weighted/route_a_ladder_finalize.json"
    run_path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))
    return out


if __name__ == "__main__":
    finalize()
