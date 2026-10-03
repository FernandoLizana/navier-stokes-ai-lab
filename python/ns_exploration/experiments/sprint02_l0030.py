"""Sprint: L-0030 full-mask triad ‖N‖ bound (C_† techo)."""

from __future__ import annotations

import json
from pathlib import Path

from ns_exploration.conjectures.l0030_triad_N_bound import lemma_l0030, save_lemma_l0030
from ns_exploration.validation.l0030_certificate import (
    build_l0030_certificate,
    save_l0030_certificate,
    verify_l0030_certificate,
)


def main() -> dict:
    b = lemma_l0030(n=24, n_grid=15)
    save_lemma_l0030(b)
    cert = build_l0030_certificate(n=24)
    ok, checks = verify_l0030_certificate(cert.as_dict())
    save_l0030_certificate(cert)
    if not ok:
        raise RuntimeError(f"L-0030 certificate failed: {checks}")
    if b.closes_c0007:
        raise RuntimeError("unexpected closure")

    c7_path = Path("conjectures/active/C-0007.json")
    c7 = json.loads(c7_path.read_text(encoding="utf-8"))
    c7["notes"] = (
        f"Open all-IC. C-R-0002 high slab. L-0027 needs C≤C_†≈{b.C_dagger:.4f}. "
        f"L-0030: R_★={b.R_star}, ‖N‖≤2√(R_★EΩ) beats Young at Ω∼E "
        f"({b.N_triad_at_equipartition:.2f} vs {b.N_young_at_E0:.2f}); "
        f"cubic C≲{b.C_from_Rstar:.1f} still ≫C_†. Need weighted/SOS triad."
    )
    related = set(c7.get("related") or [])
    related.update(
        ["C-R-0002", "L-0026", "L-0027", "L-0028", "L-0029", "L-0030"]
    )
    c7["related"] = sorted(related)
    c7_path.write_text(json.dumps(c7, indent=2), encoding="utf-8")

    out = {
        "R_star": b.R_star,
        "n_triads": b.n_triads,
        "N_triad_at_equipartition": b.N_triad_at_equipartition,
        "N_young_at_E0": b.N_young_at_E0,
        "C_from_Rstar": b.C_from_Rstar,
        "C_dagger": b.C_dagger,
        "low_slab_hybrid_worst": b.low_slab_hybrid_worst,
        "closes_c0007": False,
        "certificate_ok": ok,
        "checks": {k: bool(v) for k, v in checks.items()},
        "not_a_clay_claim": True,
    }
    Path("reports").mkdir(exist_ok=True)
    Path("experiments/exploratory/sprint02_l0030").mkdir(parents=True, exist_ok=True)
    Path("experiments/exploratory/sprint02_l0030/summary.json").write_text(
        json.dumps(out, indent=2), encoding="utf-8"
    )
    Path("reports/SPRINT_02_L0030.md").write_text(
        f"""# Sprint — L-0030 full-mask triad ‖N‖

| Item | Value |
|------|-------|
| R_★ | {b.R_star} |
| n_triads | {b.n_triads} |
| ‖N‖ triad @ Ω=E | ≈{b.N_triad_at_equipartition:.4f} |
| ‖N‖ Young @ E0 | ≈{b.N_young_at_E0:.4f} |
| C≲2√2√R_★ | ≈{b.C_from_Rstar:.4f} |
| C_† | ≈{b.C_dagger:.4f} |
| Low-slab hybrid Ω(T) | ≈{b.low_slab_hybrid_worst:.4f} |
| C-0007 | still open |

FINITE Galerkin only. Not continuum. Not Clay.
""",
        encoding="utf-8",
    )
    print(json.dumps(out, indent=2, ensure_ascii=True))
    return out


if __name__ == "__main__":
    main()
