"""Route A' cluster integral sprint (full-dealias + band pilot)."""

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

ROOT = Path(__file__).resolve().parents[3]


def _load_per_shell_manifest(path: Path) -> dict | None:
    if not path.is_file():
        return None
    return json.loads(path.read_text(encoding="utf-8"))


def run_route_a_prime(
    *,
    mode: str = "full",
    n_clusters: int = 12,
    prec: int = 128,
    workers: int = 4,
) -> dict:
    n = 24
    radii = tuple(range(1, 7)) if mode == "band" else all_dealias_radii(n)
    per_shell_path = (
        ROOT / "experiments/terminal_weighted/shell_manifest.json"
        if mode == "band"
        else ROOT / "experiments/terminal_weighted/shell_manifest_frobenius_full.json"
    )
    per_shell = _load_per_shell_manifest(per_shell_path)

    t0 = time.perf_counter()
    print(f"Route A' mode={mode} clusters={n_clusters} radii={len(radii)}", flush=True)
    manifest = build_cluster_manifest(
        n,
        radii=radii,
        n_clusters=n_clusters,
        prec=prec,
        workers=workers,
        per_shell_manifest=per_shell,
    )
    out_name = (
        "shell_manifest_cluster_band.json"
        if mode == "band"
        else "shell_manifest_cluster_full.json"
    )
    save_cluster_manifest(manifest, ROOT / f"experiments/terminal_weighted/{out_name}")
    cert = cluster_integral_certificate(manifest, prec=prec)
    elapsed = time.perf_counter() - t0

    I_lo, _ = i_star_interval(prec=prec)
    per_shell_I = None
    if per_shell:
        from ns_exploration.terminal_weighted.recompute_integral import (
            shell_integral_hi_from_blocks,
        )

        per_shell_I = shell_integral_hi_from_blocks(per_shell["blocks"], prec=prec)

    summary = {
        "mode": mode,
        "n_clusters": n_clusters,
        "n_shells": manifest["n_shells"],
        "D": manifest["D_full"],
        "integral_hi": cert["integral_hi"],
        "omega_T_hi": cert["omega_T_hi"],
        "I_star_lo": str(I_lo),
        "closes_c0008": cert["closes_c0008"],
        "per_shell_integral_hi": str(per_shell_I) if per_shell_I is not None else None,
        "improvement_vs_per_shell": (
            float(mpfr(cert["integral_hi"]) / per_shell_I)
            if per_shell_I and per_shell_I > mpfr(0)
            else None
        ),
        "phase_e_best_reference": 78.804 if mode == "full" else None,
        "runtime_s": round(elapsed, 2),
        "manifest": f"experiments/terminal_weighted/{out_name}",
    }

    out_json = ROOT / "experiments/terminal_weighted/route_a_prime_summary.json"
    out_json.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    report = ROOT / "reports/C0008_ROUTE_A_PRIME.md"
    report.write_text(_render_report(summary, cert), encoding="utf-8")
    return summary


def _render_report(summary: dict, cert: dict) -> str:
    imp = summary.get("improvement_vs_per_shell")
    imp_s = f"{imp:.4f}" if imp is not None else "n/a"
    lines = [
        "# C-0008 Route A' — Cluster Integral Certificate",
        "",
        f"Mode: **{summary['mode']}**. Clusters: {summary['n_clusters']}. "
        f"Runtime: {summary['runtime_s']}s.",
        "",
        "| Quantity | Value |",
        "|----------|-------|",
        f"| I_term^hi (Route A') | {float(summary['integral_hi']):.4f} |",
        f"| Omega(T)^hi | {float(summary['omega_T_hi']):.4f} |",
        f"| I_*^lo | {float(summary['I_star_lo']):.4f} |",
        f"| Per-shell Route A I_hi | {summary['per_shell_integral_hi'] or 'n/a'} |",
        f"| Ratio A'/A | {imp_s} |",
        f"| Closes C-0008 | {summary['closes_c0008']} |",
        "",
        "## Method",
        "",
        cert.get("honesty", ""),
        "",
        "**Conservative per cluster:**",
        "`min(sum_{r in C} B_r f(r), B_C max_{r in C} f(r))`.",
        "",
    ]
    if summary.get("phase_e_best_reference"):
        lines.append(
            f"Phase E best reference I_hi ~ {summary['phase_e_best_reference']} "
            "(min(F,1-inf) per shell)."
        )
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description="C-0008 Route A' cluster sprint")
    parser.add_argument("--mode", choices=("band", "full"), default="full")
    parser.add_argument("--clusters", type=int, default=12)
    parser.add_argument("--prec", type=int, default=128)
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    s = run_route_a_prime(
        mode=args.mode,
        n_clusters=args.clusters,
        prec=args.prec,
        workers=args.workers,
    )
    print(json.dumps(s, indent=2))


if __name__ == "__main__":
    main()
