"""SOS + triad-aware stretch sprint toward all-IC C <= C_dagger (C-0008)."""

from __future__ import annotations

import json
from pathlib import Path

from ns_exploration.conjectures.l0027_low_slab_cubic import lemma_l0027, save_lemma_l0027
from ns_exploration.conjectures.l0030_triad_N_bound import lemma_l0030, save_lemma_l0030
from ns_exploration.conjectures.l0048_fullsym_techo import C_FULLSYM_FULL_DEALIAS, lemma_l0048

ROOT = Path(__file__).resolve().parents[3]


def _try_sos_band(label: str, fn) -> dict | None:
    try:
        b = fn()
        return {
            "lemma": label,
            "C_ub_hi": getattr(
                b,
                "C_ub_hi",
                getattr(b, "C_ub_n5_hi", getattr(b, "best_sos_n5_hi", None)),
            ),
            "closes_vs_Cdagger": getattr(b, "closes_vs_Cdagger", getattr(b, "closes_vs_Cdagger", None)),
            "D": getattr(b, "D", None),
            "status": getattr(b, "status", "ok"),
        }
    except Exception as exc:
        return {"lemma": label, "error": str(exc)}


def _load_sos_certs() -> list[dict]:
    rows: list[dict] = []
    cert_dir = ROOT / "certificates"
    for path in sorted(cert_dir.glob("CERT-L005*-sos-*.json")):
        d = json.loads(path.read_text(encoding="utf-8"))
        c = d.get("C_ub_hi")
        rows.append(
            {
                "lemma": path.stem,
                "C_ub_hi": float(c) if c is not None else None,
                "closes_vs_Cdagger": d.get("closes_vs_Cdagger"),
                "D": d.get("D"),
                "source": str(path.name),
            }
        )
    return rows


def run_sprint() -> dict:
    l27 = lemma_l0027(empirical=True, n_random=8, ascent_steps=6)
    save_lemma_l0027(l27)
    l30 = lemma_l0030()
    save_lemma_l0030(l30)
    l48 = lemma_l0048()

    sos_rows: list[dict] = _load_sos_certs()
    for label, mod, fn_name in [
        ("L-0052", "l0052_sos_gram_n5", "lemma_l0052"),
        ("L-0053", "l0053_onepol_sos_gram_n5", "lemma_l0053"),
    ]:
        try:
            mod_obj = __import__(
                f"ns_exploration.conjectures.{mod}",
                fromlist=[fn_name],
            )
            fn = getattr(mod_obj, fn_name)
            row = _try_sos_band(label, fn)
            if row and "error" not in row:
                sos_rows.append(row)
        except Exception as exc:
            sos_rows.append({"lemma": label, "error": str(exc)})

    best_sos = None
    for row in sos_rows:
        c = row.get("C_ub_hi")
        if c is not None and float(c) > 0 and (best_sos is None or float(c) < best_sos):
            best_sos = float(c)

    summary = {
        "sprint": "C0008_SOS_TRIAD_STRETCH",
        "C_dagger": l27.C_dagger,
        "C_emp_max": l27.C_emp_max,
        "C_fullsym_all_ic": C_FULLSYM_FULL_DEALIAS,
        "C_from_triad_Rstar": l30.C_from_Rstar,
        "R_star": l30.R_star,
        "N_triad_at_equipartition": l30.N_triad_at_equipartition,
        "low_slab_triad_hybrid_worst": l30.low_slab_hybrid_worst,
        "low_slab_min_defect_worst": l30.low_slab_min_defect_worst,
        "best_sos_C_ub_hi": best_sos,
        "sos_rows": sos_rows,
        "all_ic_closes_C_dagger": False,
        "verdict": _verdict(l27, l30, l48, best_sos),
        "next_targets": [
            "Weighted triad / SOS flattening (L-0030 notes)",
            "Extend N5 SOS one-pol greedy beyond C-R-0014",
            "Hybrid min(triad, defect, SOS) on low slab ODE",
        ],
    }

    out = ROOT / "experiments/terminal_weighted/sos_triad_stretch_summary.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    report = ROOT / "reports/C0008_SOS_TRIAD_STRETCH.md"
    report.write_text(_render_report(summary, l27, l30), encoding="utf-8")
    return summary


def _verdict(l27, l30, l48, best_sos) -> str:
    if l48.C_fullsym_full <= l27.C_dagger:
        return "fullsym_closes"
    if l30.C_from_Rstar <= l27.C_dagger and l30.low_slab_hybrid_worst <= l27.c0007_M:
        return "triad_closes"
    if best_sos is not None and best_sos > 0 and best_sos <= l27.C_dagger:
        return "sos_subclass_only"
    return "all_ic_open"


def _render_report(summary: dict, l27, l30) -> str:
    lines = [
        "# C-0008 SOS / Triad Stretch Sprint",
        "",
        f"**Verdict:** `{summary['verdict']}` — all-IC C <= C_dagger remains **open**.",
        "",
        "## Thresholds",
        "",
        f"| Quantity | Value | vs C_dagger ({summary['C_dagger']:.4f}) |",
        "|----------|-------|------------------|",
        f"| C_emp (N2) | {summary['C_emp_max']:.6f} | yes (empirical) |",
        f"| C_fullsym all-IC | {summary['C_fullsym_all_ic']:.4f} | **no** |",
        f"| C from triad R_star | {summary['C_from_triad_Rstar']:.2f} | **no** |",
        f"| Best SOS C_ub_hi | {summary['best_sos_C_ub_hi'] or 'n/a'} | subclass |",
        "",
        "## L-0030 triad",
        "",
        f"- R_star = {summary['R_star']}",
        f"- Low-slab hybrid Omega(T) worst ~ {summary['low_slab_triad_hybrid_worst']:.2f}",
        f"- min(hybrid, defect) worst ~ {summary['low_slab_min_defect_worst']:.2f}",
        "",
        l30.notes,
        "",
        "## SOS bands attempted",
        "",
        "| Lemma | C_ub_hi | vs C_dagger |",
        "|-------|---------|-------------|",
    ]
    for row in summary["sos_rows"]:
        c = row.get("C_ub_hi", row.get("error", "-"))
        cd = "yes" if row.get("closes_vs_Cdagger") else ("no" if "error" not in row else "err")
        lines.append(f"| {row['lemma']} | {c} | {cd} |")
    lines.extend(["", "## Next targets", ""])
    for i, t in enumerate(summary["next_targets"], 1):
        lines.append(f"{i}. {t}")
    lines.extend(["", l27.notes, ""])
    return "\n".join(lines)


if __name__ == "__main__":
    print(json.dumps(run_sprint(), indent=2))
