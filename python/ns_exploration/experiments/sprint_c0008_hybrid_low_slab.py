"""C-0008 hybrid low-slab sprint: min(triad, defect, weighted, SOS cubic ODE)."""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

from ns_exploration.conjectures.l0023_cubic_dissipation import omega_ode_bound
from ns_exploration.conjectures.l0027_low_slab_cubic import lemma_l0027, save_lemma_l0027
from ns_exploration.conjectures.l0030_triad_N_bound import lemma_l0030, save_lemma_l0030
from ns_exploration.conjectures.l0031_weighted_triad import lemma_l0031, save_lemma_l0031
from ns_exploration.conjectures.l0048_fullsym_techo import C_FULLSYM_FULL_DEALIAS
from ns_exploration.terminal_weighted.constants import FROZEN

ROOT = Path(__file__).resolve().parents[3]


def _load_best_sos() -> tuple[float | None, float | None, list[dict]]:
    rows: list[dict] = []
    best_band = None
    best_greedy = None
    for path in sorted((ROOT / "certificates").glob("CERT-L*.json")):
        name = path.name.lower()
        if "sos" not in name and "onepol" not in name:
            continue
        d = json.loads(path.read_text(encoding="utf-8"))
        c = d.get("C_ub_hi") or d.get("C_shor")
        if c is None:
            continue
        c = float(c)
        row = {"cert": path.name, "C_ub_hi": c, "D": d.get("D")}
        rows.append(row)
        if "onepol-greedy" in name or "cr0008" in name or "greedy-second" in name:
            if best_greedy is None or c < best_greedy:
                best_greedy = c
        if "band123" in name and "1234" not in name:
            if best_band is None or c < best_band:
                best_band = c
    if best_band is None and rows:
        best_band = min(r["C_ub_hi"] for r in rows if r["C_ub_hi"] > 0)
    if best_greedy is None:
        best_greedy = max(
            (r["C_ub_hi"] for r in rows if "onepol" in r["cert"].lower()),
            default=best_band or 0.0,
        )
    return best_band, best_greedy, rows


def _sos_cubic_low_slab_worst(
    C_ub: float,
    *,
    E0: float,
    nu: float,
    T: float,
    omega_lo: float,
    omega_hi: float,
    n_grid: int = 25,
) -> float:
    xs = np.linspace(omega_lo, omega_hi, n_grid)
    return max(omega_ode_bound(C_ub, E0, nu, T, float(o)) for o in xs)


def run(*, n_grid: int = 15) -> dict:
    l27 = lemma_l0027(empirical=True, n_random=8, ascent_steps=6)
    save_lemma_l0027(l27)
    l30 = lemma_l0030(n_grid=n_grid)
    save_lemma_l0030(l30)
    l31 = lemma_l0031(weight="inv_sqrt_rp_rq")
    save_lemma_l0031(l31)

    best_band_sos, best_greedy_sos, sos_rows = _load_best_sos()
    sos_band_worst = None
    sos_greedy_worst = None
    if best_band_sos:
        sos_band_worst = _sos_cubic_low_slab_worst(
            best_band_sos,
            E0=l27.E0,
            nu=l27.nu,
            T=l27.T,
            omega_lo=l27.E0,
            omega_hi=l27.Omega_star,
        )
    if best_greedy_sos:
        sos_greedy_worst = _sos_cubic_low_slab_worst(
            best_greedy_sos,
            E0=l27.E0,
            nu=l27.nu,
            T=l27.T,
            omega_lo=l27.E0,
            omega_hi=l27.Omega_star,
        )

    candidates: dict[str, float] = {
        "triad_hybrid_ode": l30.low_slab_hybrid_worst,
        "triad_min_defect_ode": l30.low_slab_min_defect_worst,
        "weighted_min_defect_ode": l31.low_slab_hybrid_weighted_worst,
    }
    if sos_band_worst is not None:
        candidates["sos_band123_cubic_ode"] = sos_band_worst
    if sos_greedy_worst is not None:
        candidates["sos_greedy_onepol_cubic_ode"] = sos_greedy_worst

    best_label = min(candidates, key=candidates.get)
    best_omega = candidates[best_label]
    best_all_ic_label = min(
        (k for k in candidates if not k.startswith("sos_")),
        key=lambda k: candidates[k],
    )
    best_all_ic_omega = candidates[best_all_ic_label]
    M_sharp = FROZEN.c0008_target_M
    M_c0007 = l27.c0007_M
    subclass_only = best_label.startswith("sos_")

    summary = {
        "sprint": "C0008_hybrid_low_slab",
        "C_dagger": l27.C_dagger,
        "C_emp_max": l27.C_emp_max,
        "C_fullsym_all_ic": C_FULLSYM_FULL_DEALIAS,
        "Omega_star": l27.Omega_star,
        "M_c0007": M_c0007,
        "M_sharp_c0008": M_sharp,
        "L0073_terminal_omega_hi": 54.231,
        "candidates": candidates,
        "best_label": best_label,
        "best_low_slab_omega_T": best_omega,
        "best_all_ic_label": best_all_ic_label,
        "best_all_ic_low_slab_omega_T": best_all_ic_omega,
        "best_is_subclass_sos": subclass_only,
        "closes_c0007_low_slab_all_ic": best_all_ic_omega <= M_c0007,
        "closes_c0008_sharp": False,
        "closes_c0008_sharp_subclass_low_slab_only": subclass_only and best_omega <= M_sharp,
        "gap_vs_M_sharp_all_ic": best_all_ic_omega - M_sharp,
        "gap_vs_M_sharp": best_omega - M_sharp,
        "best_sos_band_C_ub": best_band_sos,
        "best_sos_greedy_C_ub": best_greedy_sos,
        "sos_cert_count": len(sos_rows),
        "all_ic_closes_C_dagger": False,
        "verdict": _verdict(best_all_ic_omega, M_sharp, M_c0007),
        "subclass_sos_verdict": (
            "subclass_sos_low_slab_below_M"
            if subclass_only and best_omega <= M_sharp
            else None
        ),
        "viable_next": [
            "Prove all-IC C <= C_dagger (fullsym techo blocks naive Shor)",
            "Terminal non-partition structure (L-0074 refuted finer clusters)",
            "Cluster greedy SOS D=1772 (subclass, ~113h deferred)",
        ],
    }

    out = ROOT / "experiments/terminal_weighted/hybrid_low_slab_summary.json"
    out.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    reg = ROOT / "conjectures/active/L-0075.json"
    reg.write_text(
        json.dumps(
            {
                "lemma_id": "L-0075",
                "route": "B",
                "status": "best_low_slab_hybrid",
                "evidence_level": "N7",
                "parent": "L-0027",
                "claim": "C-0008",
                "best_label": best_label,
                "best_low_slab_omega_T": round(best_omega, 4),
                "M_sharp": M_sharp,
                "M_c0007": M_c0007,
                "best_all_ic_label": best_all_ic_label,
                "best_all_ic_low_slab_omega_T": round(best_all_ic_omega, 4),
                "closes_c0008_sharp": False,
                "closes_c0008_sharp_subclass_low_slab_only": subclass_only
                and best_omega <= M_sharp,
                "best_is_subclass_sos": subclass_only,
                "summary": str(out.relative_to(ROOT)).replace("\\", "/"),
                "notes": (
                    "Hybrid min(triad/defect/weighted/SOS cubic ODE) on low slab. "
                    "Subclass SOS; does not close sharp M or all-IC C<=C_dagger."
                ),
            },
            indent=2,
        ),
        encoding="utf-8",
    )

    report = ROOT / "reports/C0008_HYBRID_LOW_SLAB.md"
    report.write_text(_render(summary), encoding="utf-8")
    return summary


def _verdict(best_omega: float, M_sharp: float, M_c0007: float) -> str:
    if best_omega <= M_sharp:
        return "closes_c0008_sharp_all_ic"
    if best_omega <= M_c0007:
        return "closes_c0007_not_sharp"
    if best_omega < 54.0:
        return "low_slab_beats_terminal_sketch"
    return "low_slab_open"


def _render(summary: dict) -> str:
    lines = [
        "# C-0008 Hybrid Low-Slab Sprint",
        "",
        f"**Verdict:** `{summary['verdict']}`",
        "",
        f"**Best (all-IC sketch):** `{summary['best_all_ic_label']}` — Ω ≤ **{summary['best_all_ic_low_slab_omega_T']:.4f}**",
        f"**Best (subclass SOS):** `{summary['best_label']}` — Ω ≤ **{summary['best_low_slab_omega_T']:.4f}**",
        f"(M_sharp={summary['M_sharp_c0008']:.4f}, M_c0007={summary['M_c0007']:.4f}, "
        f"L-0073 terminal Ω≈{summary['L0073_terminal_omega_hi']:.3f})",
        "",
        "| Route | Ω(T) worst (low slab) |",
        "|-------|----------------------|",
    ]
    for k, v in sorted(summary["candidates"].items(), key=lambda kv: kv[1]):
        mark = " **" if k == summary["best_label"] else ""
        end = "**" if k == summary["best_label"] else ""
        lines.append(f"| {mark}{k}{end} | {v:.4f} |")
    lines.extend(
        [
            "",
            f"Gap vs M_sharp (all-IC): **{summary['gap_vs_M_sharp_all_ic']:.4f}**",
            "",
            "## SOS cubic (subclass)",
            "",
            f"- Best band SOS C_ub: {summary['best_sos_band_C_ub']}",
            f"- Best greedy one-pol C_ub: {summary['best_sos_greedy_C_ub']}",
            "",
            "## Next",
            "",
        ]
    )
    for i, v in enumerate(summary["viable_next"], 1):
        lines.append(f"{i}. {v}")
    return "\n".join(lines)


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
