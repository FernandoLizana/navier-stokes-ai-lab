"""Sprint: L-0026 case-split techo for C-0007 (gap remains open)."""

from __future__ import annotations

import json
from pathlib import Path

from ns_exploration.conjectures.l0026_c0007_techo import lemma_l0026, save_lemma_l0026
from ns_exploration.validation.l0026_certificate import (
    build_l0026_certificate,
    save_l0026_certificate,
    verify_l0026_certificate,
)


def main() -> dict:
    b = lemma_l0026(n=24)
    save_lemma_l0026(b)
    cert = build_l0026_certificate(n=24)
    ok, checks = verify_l0026_certificate(cert.as_dict())
    save_l0026_certificate(cert)
    if not b.high_slab_closes:
        raise RuntimeError("L-0026 high-slab case-split failed")
    if b.closes_c0007:
        raise RuntimeError("unexpected: L-0026 claims C-0007 closed")
    if not ok:
        raise RuntimeError(f"L-0026 certificate failed: {checks}")

    # Refresh active C-0007 notes with L-0026 pointer (do not prove).
    c7_path = Path("conjectures/active/C-0007.json")
    c7 = json.loads(c7_path.read_text(encoding="utf-8"))
    c7["notes"] = (
        f"Open: M≈{b.c0007_M:.6f} below L-0024 majorant {b.L0024_majorant:.6f}. "
        f"L-0026: high slab Ω0≥Ω★≈{b.Omega_star:.4f} already ≤M via L-0024; "
        f"low slab Ω0∈[E0,Ω★) is the obstruction (structural techo on defect ODE). "
        f"Certificate: CERT-L0026-c0007-techo-N24."
    )
    c7["lemma_techo"] = "L-0026"
    c7["Omega_star"] = b.Omega_star
    c7_path.write_text(json.dumps(c7, indent=2), encoding="utf-8")

    out = {
        "C0007_M": b.c0007_M,
        "L0024_majorant": b.L0024_majorant,
        "Omega_star": b.Omega_star,
        "high_slab_closes": b.high_slab_closes,
        "closes_c0007": b.closes_c0007,
        "certificate_ok": ok,
        "checks": {k: bool(v) for k, v in checks.items()},
        "not_a_clay_claim": True,
    }
    Path("reports").mkdir(exist_ok=True)
    Path("experiments/exploratory/sprint02_l0026").mkdir(parents=True, exist_ok=True)
    Path("experiments/exploratory/sprint02_l0026/summary.json").write_text(
        json.dumps(out, indent=2), encoding="utf-8"
    )
    report = f"""# Sprint — L-0026 techo for C-0007

| Item | Result |
|------|--------|
| **C-0007** (N≤24, M≈{b.c0007_M:.6f}) | **still open** |
| L-0024 majorant | {b.L0024_majorant:.6f} |
| Stokes floor | {b.Stokes_floor:.6f} |
| Ω★ (high slab) | {b.Omega_star:.6f} |
| High slab Ω0≥Ω★ | **closes** (≤M via L-0024) |
| Low slab | obstruction; defect-ODE techo |

FINITE Galerkin only. Not continuum. Not Clay.
"""
    Path("reports/SPRINT_02_L0026.md").write_text(report, encoding="utf-8")
    print(json.dumps(out, indent=2, ensure_ascii=True))
    return out


if __name__ == "__main__":
    main()
