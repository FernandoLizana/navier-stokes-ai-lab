"""N-ladder pilot: terminal-weighted band Frobenius at n in {12, 16, 20, 24}."""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

from gmpy2 import mpfr

from ns_exploration.conjectures.l0048_matrixfree_fullsym import all_dealias_radii
from ns_exploration.terminal_weighted.certify_shell_integral import shell_integral_certificate
from ns_exploration.terminal_weighted.constants import FROZEN
from ns_exploration.terminal_weighted.intervals import i_star_interval
from ns_exploration.terminal_weighted.n_coverage import n_domination_analysis, wavevector_shells_subset
from ns_exploration.terminal_weighted.shell_decomposition import build_shell_manifest, save_shell_manifest
from ns_exploration.validation.l0003_certificate import full_dealias_exact_stats

ROOT = Path(__file__).resolve().parents[3]


def run_ladder(
    ns: tuple[int, ...] = (12, 16, 20, 24),
    *,
    mode: str = "band",
    prec: int = 128,
    workers: int = 1,
) -> dict:
    rows = []
    t0_all = time.perf_counter()
    for n in ns:
        radii = tuple(range(1, 7)) if mode == "band" else all_dealias_radii(n)
        K2, M = full_dealias_exact_stats(n)
        embed = wavevector_shells_subset(n, 24) if n < 24 else {"subset": True, "n_small": n, "n_large": 24}
        t0 = time.perf_counter()
        manifest = build_shell_manifest(
            n, prec=prec, radii=radii, workers=workers, bound_method="frobenius"
        )
        save_shell_manifest(
            manifest, ROOT / f"experiments/terminal_weighted/shell_manifest_n{n}_frobenius.json"
        )
        cert = shell_integral_certificate(n, prec=prec, manifest=manifest)
        I_lo, _ = i_star_interval(prec=prec)
        elapsed = time.perf_counter() - t0
        rows.append(
            {
                "n": n,
                "mode": mode,
                "D": manifest["D_full"],
                "n_shells": manifest["n_shells"],
                "K2": K2,
                "M_modes": M,
                "wavevectors_embed_in_n24": embed.get("subset"),
                "I_term_hi": cert["integral_hi"],
                "omega_T_hi": cert["omega_T_hi"],
                "closes_c0008": cert["closes_c0008"],
                "I_star_lo": str(I_lo),
                "beats_I_star": mpfr(cert["integral_hi"]) < I_lo,
                "runtime_s": round(elapsed, 2),
            }
        )
        print(
            f"n={n} D={manifest['D_full']} shells={manifest['n_shells']} "
            f"I_hi={cert['integral_hi'][:12]}... t={elapsed:.1f}s",
            flush=True,
        )

    summary = {
        "mode": mode,
        "ns": list(ns),
        "M_target": FROZEN.c0008_target_M,
        "n_domination": n_domination_analysis(24),
        "rows": rows,
        "total_runtime_s": round(time.perf_counter() - t0_all, 2),
        "note": (
            "Band pilot per n; full-dealias at n<24 deferred (cost). "
            "Wavevector embedding into n=24 verified for domination sketch."
        ),
    }
    out = ROOT / "experiments/terminal_weighted/n_ladder_summary.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    report = ROOT / "reports/C0008_N_LADDER.md"
    report.write_text(_render_report(summary), encoding="utf-8")
    return summary


def _render_report(summary: dict) -> str:
    lines = [
        "# C-0008 N-Ladder Pilot",
        "",
        f"Mode: **{summary['mode']}**. Total runtime: {summary['total_runtime_s']}s.",
        "",
        "| n | D | shells | I_term^hi | Omega(T)^hi | embed in n=24 |",
        "|---|---|--------|-----------|-------------|---------------|",
    ]
    for r in summary["rows"]:
        lines.append(
            f"| {r['n']} | {r['D']} | {r['n_shells']} | {float(r['I_term_hi']):.4f} | "
            f"{float(r['omega_T_hi']):.4f} | {r['wavevectors_embed_in_n24']} |"
        )
    lines.extend(
        [
            "",
            "## Domination sketch (n=24 majorant)",
            "",
            summary["n_domination"]["conclusion"],
            "",
            "**Note:** Band results are pilots only; do not substitute for full-dealias certification.",
            "",
        ]
    )
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description="C-0008 N-ladder band pilot")
    parser.add_argument("--mode", choices=("band", "full"), default="band")
    parser.add_argument("--prec", type=int, default=128)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--ns", type=int, nargs="+", default=[12, 16, 20, 24])
    args = parser.parse_args()
    run_ladder(tuple(args.ns), mode=args.mode, prec=args.prec, workers=args.workers)


if __name__ == "__main__":
    main()
