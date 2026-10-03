"""Finalize Route A' cluster-best certificate from manifest."""

from __future__ import annotations

import json
import sys
from pathlib import Path

from gmpy2 import mpfr

from ns_exploration.terminal_weighted.cluster_bounds import cluster_integral_certificate
from ns_exploration.terminal_weighted.constants import FROZEN
from ns_exploration.terminal_weighted.intervals import i_star_interval, to_float_hi
from ns_exploration.terminal_weighted.l0048_audit import l0048_value_classification
from ns_exploration.terminal_weighted.n_coverage import n_domination_analysis
from ns_exploration.terminal_weighted.weights import stokes_floor_hi

ROOT = Path(__file__).resolve().parents[3]
DEFAULT_MANIFEST = ROOT / "experiments/terminal_weighted/shell_manifest_cluster_best_full.json"
DEFAULT_CERT = ROOT / "certificates/CERT-L0073-C0008-route-a-prime-best.json"


def finalize(manifest_path: str | Path = DEFAULT_MANIFEST, *, prec: int = 128) -> dict:
    manifest_path = Path(manifest_path)
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    cert_body = cluster_integral_certificate(manifest, prec=prec)
    I_lo, I_hi = i_star_interval(prec=prec)

    cert = {
        "claim": "C-0008",
        "status": "candidate_certificate",
        "lemma": "L-0073",
        "parameters": {
            "nu": FROZEN.nu,
            "E0": FROZEN.E0,
            "T": FROZEN.T,
            "N_max": manifest["n"],
            "M_target": FROZEN.c0008_target_M,
        },
        "normalization": {"fft": "fftn/N^3", "hermitian": "sprintB_physical_tensor"},
        "mode_set": {
            "n": manifest["n"],
            "D": manifest["D_full"],
            "n_shells": manifest["n_shells"],
            "n_clusters": manifest["n_clusters"],
            "radii_mode": "full",
        },
        "arithmetic_backend": manifest.get("arithmetic_backend", {"precision_bits": prec}),
        "stokes_floor_interval": [FROZEN.stokes_floor, to_float_hi(mpfr(stokes_floor_hi(24)))],
        "i_star_interval": [str(I_lo), str(I_hi)],
        "method": "route_a_prime_cluster_best",
        "operator_bound_method": manifest.get("bound_method", "best"),
        "integral_interval": ["0", cert_body["integral_hi"]],
        "omega_T_interval": ["0", cert_body["omega_T_hi"]],
        "strict_margin": str(I_lo - mpfr(cert_body["integral_hi"])),
        "closes_strict": False,
        "source_hashes": {"cluster_manifest": manifest.get("manifest_sha256")},
        "limitations": [
            "Route A': min(per-shell Phase-E best sum, B_cluster max f) per cluster",
            "B_cluster = min(Frobenius, 1-inf) on combined cluster output",
            "Triangle inequality across clusters remains",
        ],
        "l0048_audit": l0048_value_classification(),
        "n_coverage": n_domination_analysis(manifest["n"]),
        "route_a_prime": cert_body,
        "verification_command": (
            "python verify_c0008_cluster_certificate.py "
            "certificates/CERT-L0073-C0008-route-a-prime-best.json "
            "experiments/terminal_weighted/shell_manifest_cluster_best_full.json"
        ),
    }

    DEFAULT_CERT.parent.mkdir(parents=True, exist_ok=True)
    DEFAULT_CERT.write_text(json.dumps(cert, indent=2), encoding="utf-8")

    sys.path.insert(0, str(ROOT))
    from verify_c0008_cluster_certificate import verify_certificate

    vok, vmsgs = verify_certificate(DEFAULT_CERT, manifest_path=manifest_path)
    out = {
        "certificate": str(DEFAULT_CERT.relative_to(ROOT)).replace("\\", "/"),
        "manifest": str(manifest_path.relative_to(ROOT)).replace("\\", "/"),
        "integral_hi": cert_body["integral_hi"],
        "omega_T_hi": cert_body["omega_T_hi"],
        "I_star_lo": str(I_lo),
        "verifier_pass": vok,
        "verifier_msgs": vmsgs,
        "closes_c0008": cert_body["closes_c0008"],
        "c0008_status": "exploring",
    }
    run_path = ROOT / "experiments/terminal_weighted/route_a_prime_cert_run.json"
    run_path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))
    return out


if __name__ == "__main__":
    mp = sys.argv[1] if len(sys.argv) > 1 else str(DEFAULT_MANIFEST)
    finalize(mp)
