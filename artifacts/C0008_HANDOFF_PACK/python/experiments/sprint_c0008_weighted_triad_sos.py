"""Weighted triad (L-0031) + full SOS ladder audit toward all-IC C <= C_dagger."""

from __future__ import annotations

import json
from pathlib import Path

from ns_exploration.conjectures.l0030_triad_N_bound import lemma_l0030, save_lemma_l0030
from ns_exploration.conjectures.l0031_weighted_triad import lemma_l0031, save_lemma_l0031
from ns_exploration.conjectures.l0027_low_slab_cubic import lemma_l0027, save_lemma_l0027
from ns_exploration.conjectures.l0048_fullsym_techo import C_FULLSYM_FULL_DEALIAS

ROOT = Path(__file__).resolve().parents[3]


def _load_all_sos_certs() -> list[dict]:
    rows = []
    for path in sorted((ROOT / "certificates").glob("CERT-L*.json")):
        if "sos" not in path.name.lower() and "onepol" not in path.name.lower():
            continue
        d = json.loads(path.read_text(encoding="utf-8"))
        c = d.get("C_ub_hi") or d.get("C_shor")
        if c is None:
            continue
        rows.append(
            {
                "cert": path.name,
                "C_ub_hi": float(c),
                "D": d.get("D"),
                "closes_vs_Cdagger": d.get("closes_vs_Cdagger"),
                "subclass": d.get("subclass") or d.get("band_radii"),
            }
        )
    return rows


def run() -> dict:
    l27 = lemma_l0027(empirical=True, n_random=8, ascent_steps=6)
    save_lemma_l0027(l27)
    l30 = lemma_l0030()
    save_lemma_l0030(l30)
    l31 = lemma_l0031()
    l31b = lemma_l0031(weight="inv_sqrt_rp_rq")
    save_lemma_l0031(l31)

    sos = _load_all_sos_certs()
    best_sos = min((r["C_ub_hi"] for r in sos if r["C_ub_hi"] > 0), default=None)
    largest_sos = max((r["C_ub_hi"] for r in sos), default=None)

    # Simple extrapolation: fit log(C) vs log(D) on SOS ladder
    pts = [(r["D"], r["C_ub_hi"]) for r in sos if r.get("D") and r["C_ub_hi"] > 0]
    extrap_full = None
    if len(pts) >= 3:
        import math

        xs = [math.log(float(d)) for d, _ in pts]
        ys = [math.log(float(c)) for _, c in pts]
        mx = sum(xs) / len(xs)
        my = sum(ys) / len(ys)
        num = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
        den = sum((x - mx) ** 2 for x in xs) or 1.0
        slope = num / den
        intercept = my - slope * mx
        extrap_full = math.exp(intercept + slope * math.log(6748.0))

    summary = {
        "sprint": "weighted_triad_sos",
        "C_dagger": l27.C_dagger,
        "C_fullsym_all_ic": C_FULLSYM_FULL_DEALIAS,
        "R_star_unweighted": l30.R_star,
        "R_star_weighted_inv_sqrt": l31.R_star_weighted,
        "R_star_weighted_inv_sqrt_rp_rq": l31b.R_star_weighted,
        "C_from_weighted_triad_rp_rq": l31b.C_from_weighted_R,
        "C_from_unweighted_triad": l30.C_from_Rstar,
        "C_from_weighted_triad": l31.C_from_weighted_R,
        "low_slab_weighted_worst": l31.low_slab_hybrid_weighted_worst,
        "best_sos_C_ub_hi": best_sos,
        "largest_sos_C_ub_hi": largest_sos,
        "sos_extrap_D6748": extrap_full,
        "sos_cert_count": len(sos),
        "sos_rows": sos,
        "all_ic_closes_C_dagger": False,
        "verdict": (
            "all_ic_open"
            if (l31.C_from_weighted_R > l27.C_dagger and C_FULLSYM_FULL_DEALIAS > l27.C_dagger)
            else "needs_review"
        ),
    }

    out = ROOT / "experiments/terminal_weighted/weighted_triad_sos_summary.json"
    out.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    report = ROOT / "reports/C0008_WEIGHTED_TRIAD_SOS.md"
    report.write_text(_render(summary, l31), encoding="utf-8")
    return summary


def _render(summary: dict, l31) -> str:
    lines = [
        "# C-0008 Weighted Triad + SOS Ladder",
        "",
        f"**Verdict:** `{summary['verdict']}`",
        "",
        "## Triad constants (N=24)",
        "",
        f"| Variant | R | C bound | vs C_dagger |",
        f"|---------|---|---------|-------------|",
        f"| Unweighted R_★ | {summary['R_star_unweighted']} | {summary['C_from_unweighted_triad']:.2f} | no |",
        f"| Weighted 1/√r_s | {summary['R_star_weighted_inv_sqrt']:.2f} | {summary['C_from_weighted_triad']:.2f} | no |",
        f"| Fullsym Shor | — | {summary['C_fullsym_all_ic']:.2f} | no |",
        "",
        f"Low-slab weighted+defect Ω(T) worst: **{summary['low_slab_weighted_worst']:.2f}**",
        "",
        "## SOS ladder",
        "",
        f"- Certificates scanned: {summary['sos_cert_count']}",
        f"- Best C_ub_hi: **{summary['best_sos_C_ub_hi']}** (subclass)",
        f"- Largest certified band: **{summary['largest_sos_C_ub_hi']}**",
        f"- log-log extrap to D=6748: **{summary['sos_extrap_D6748']:.2f}**"
        if summary.get("sos_extrap_D6748")
        else "- Extrap: n/a",
        "",
        l31.notes,
        "",
    ]
    return "\n".join(lines)


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
