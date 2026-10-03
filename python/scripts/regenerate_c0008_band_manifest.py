"""Regenerate band repair manifest with one_inf_only + cluster metadata."""

from __future__ import annotations

import json
from pathlib import Path

from ns_exploration.terminal_weighted.cluster_residual_prototype import (
    run_band123_cluster_residual,
    run_band6_cluster_residual,
)
from ns_exploration.terminal_weighted.opnorm_certified import certified_shell_blocks

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "experiments" / "terminal_weighted" / "shell_manifest_repair_band6_one_inf.json"


def regenerate_band6_manifest(*, n: int = 24, prec: int = 128, rational_coefficients: bool = False) -> dict:
    radii = tuple(range(1, 7))
    blocks = certified_shell_blocks(
        n,
        radii,
        prec=prec,
        workers=1,
        bound_method="one_inf_only",
        rational_coefficients=rational_coefficients,
    )
    cluster = run_band6_cluster_residual(n=n, prec=prec)
    tag = "C0008_cert_repair_phase7_rational" if rational_coefficients else "C0008_cert_repair_phase7"
    payload = {
        **blocks,
        "repair_tag": tag,
        "cluster_residual_L0073R": cluster["residual_combined"],
        "B_C_combined_hi": cluster["B_C_combined_hi"],
    }
    out = OUT
    if rational_coefficients:
        out = REPO / "experiments/terminal_weighted/shell_manifest_repair_band6_one_inf_rational.json"
    OUT.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    payload["_manifest_path"] = str(out)
    return payload


def main() -> dict:
    import argparse

    p = argparse.ArgumentParser()
    p.add_argument("--rational", action="store_true", help="Use mpfr_leray production path")
    p.add_argument("--prec", type=int, default=128)
    args = p.parse_args()
    rep123 = run_band123_cluster_residual()
    rep6 = regenerate_band6_manifest(prec=args.prec, rational_coefficients=args.rational)
    print(
        {
            "band123_cluster": rep123,
            "manifest_out": rep6.get("_manifest_path", str(OUT)),
            "n_shells": rep6["n_shells"],
            "rational_coefficients": rep6.get("rational_coefficients", False),
        }
    )
    return rep6


if __name__ == "__main__":
    main()
