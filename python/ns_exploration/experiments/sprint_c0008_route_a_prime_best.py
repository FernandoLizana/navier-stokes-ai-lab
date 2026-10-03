"""Route A' with cluster min(F,1-inf) + Phase E per-shell best bounds."""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

from gmpy2 import mpfr

from ns_exploration.conjectures.l0048_matrixfree_fullsym import all_dealias_radii
from ns_exploration.terminal_weighted.cluster_bounds import (
    build_cluster_manifest,
    cluster_integral_certificate,
    save_cluster_manifest,
)
from ns_exploration.terminal_weighted.constants import FROZEN
from ns_exploration.terminal_weighted.intervals import i_star_interval
from ns_exploration.terminal_weighted.recompute_integral import shell_integral_hi_from_blocks

ROOT = Path(__file__).resolve().parents[3]


def run(
    *,
    mode: str = "full",
    n_clusters: int = 12,
    prec: int = 128,
    workers: int = 4,
) -> dict:
    n = 24
    radii = tuple(range(1, 7)) if mode == "band" else all_dealias_radii(n)
    best_path = (
        ROOT / "experiments/terminal_weighted/shell_manifest.json"
        if mode == "band"
        else ROOT / "experiments/terminal_weighted/shell_manifest_best_full.json"
    )
    per_shell = json.loads(best_path.read_text(encoding="utf-8"))

    t0 = time.perf_counter()
    print(f"Route A' BEST mode={mode} clusters={n_clusters}", flush=True)
    manifest = build_cluster_manifest(
        n,
        radii=radii,
        n_clusters=n_clusters,
        prec=prec,
        workers=workers,
        per_shell_manifest=per_shell,
        bound_method="best",
    )
    suffix = "band" if mode == "band" else "full"
    out = ROOT / f"experiments/terminal_weighted/shell_manifest_cluster_best_{suffix}.json"
    save_cluster_manifest(manifest, out)
    cert = cluster_integral_certificate(manifest, prec=prec)

    I_phase_e = shell_integral_hi_from_blocks(per_shell["blocks"], prec=prec)
    I_lo, _ = i_star_interval(prec=prec)
    elapsed = time.perf_counter() - t0

    summary = {
        "mode": mode,
        "bound_method": "best",
        "n_clusters": n_clusters,
        "integral_hi": cert["integral_hi"],
        "omega_T_hi": cert["omega_T_hi"],
        "phase_e_per_shell_I_hi": str(I_phase_e),
        "improvement_vs_phase_e": float(mpfr(cert["integral_hi"]) / I_phase_e)
        if I_phase_e > mpfr(0)
        else None,
        "closes_c0008": cert["closes_c0008"],
        "I_star_lo": str(I_lo),
        "runtime_s": round(elapsed, 2),
        "manifest": str(out.relative_to(ROOT)).replace("\\", "/"),
    }

    out_json = ROOT / "experiments/terminal_weighted/route_a_prime_best_summary.json"
    out_json.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    report = ROOT / "reports/C0008_ROUTE_A_PRIME_BEST.md"
    report.write_text(_render(summary, cert), encoding="utf-8")
    return summary


def _render(summary: dict, cert: dict) -> str:
    imp = summary.get("improvement_vs_phase_e")
    return "\n".join(
        [
            "# C-0008 Route A' — Cluster + min(F,1-inf)",
            "",
            f"Mode **{summary['mode']}**, {summary['n_clusters']} clusters, "
            f"runtime {summary['runtime_s']}s.",
            "",
            "| Quantity | Value |",
            "|----------|-------|",
            f"| I_term^hi (Route A' best) | {float(summary['integral_hi']):.4f} |",
            f"| Phase E per-shell I_hi | {float(summary['phase_e_per_shell_I_hi']):.4f} |",
            f"| Ratio A'/Phase E | {imp:.4f} |" if imp else "| Ratio | n/a |",
            f"| Omega(T)^hi | {float(summary['omega_T_hi']):.4f} |",
            f"| Closes C-0008 | {summary['closes_c0008']} |",
            "",
            cert.get("honesty", ""),
            "",
        ]
    )


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--mode", choices=("band", "full"), default="full")
    p.add_argument("--clusters", type=int, default=12)
    p.add_argument("--prec", type=int, default=128)
    p.add_argument("--workers", type=int, default=4)
    args = p.parse_args()
    print(json.dumps(run(mode=args.mode, n_clusters=args.clusters, prec=args.prec, workers=args.workers), indent=2))


if __name__ == "__main__":
    main()
