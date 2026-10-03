"""Sprint: L-0027 low-slab cubic + C-R-0002 high-slab proved."""

from __future__ import annotations

import json
from pathlib import Path

from ns_exploration.conjectures.l0026_c0007_techo import lemma_l0026
from ns_exploration.conjectures.l0027_low_slab_cubic import lemma_l0027, save_lemma_l0027
from ns_exploration.validation.l0027_certificate import (
    build_l0027_certificate,
    save_l0027_certificate,
    verify_l0027_certificate,
)


def main() -> dict:
    b26 = lemma_l0026(n=24, n_grid=61)
    b = lemma_l0027(n=24, empirical=True, n_random=12, ascent_steps=8)
    save_lemma_l0027(b)
    cert = build_l0027_certificate(n=24)
    ok, checks = verify_l0027_certificate(cert.as_dict())
    save_l0027_certificate(cert)
    if not ok:
        raise RuntimeError(f"L-0027 certificate failed: {checks}")
    if b.proved_closes_c0007:
        raise RuntimeError("unexpected proved closure")
    if b.all_ic_cubic_closes_c0007:
        raise RuntimeError("unexpected all-IC cubic closure")

    # C-R-0002: high-Ω0 restriction of C-0007 — proved by L-0024+L-0026
    cr = {
        "id": "C-R-0002",
        "route": "B",
        "status": "proved_restricted",
        "evidence_level": "N7",
        "lemma": "L-0024+L-0026",
        "domain": (
            "Galerkin/pseudospectral T^3, N<=24, 2/3 dealias, "
            "IC with Omega(0)>=Omega_star. NOT continuum."
        ),
        "nu": 0.1,
        "energy": 0.5,
        "t_end": 0.02,
        "n_max": 24,
        "proposed_bound_M": b26.c0007_M,
        "Omega_star": b26.Omega_star,
        "parent_open": "C-0007",
        "statement": (
            f"For all divergence-free Fourier fields on N≤24 with energy ≤0.5 and "
            f"Ω(0)≥Ω★≈{b26.Omega_star:.6f}, dealiased Galerkin NS (ν=0.1) satisfies "
            f"Ω(0.02)≤{b26.c0007_M:.8f}. Proved by L-0024 spectral-defect ODE "
            f"(monotone in Ω0; L-0026 case-split). FINITE-DIMENSIONAL only."
        ),
        "clay_implication": (
            "None. High-Ω0 Galerkin restriction only; not all-IC C-0007; "
            "not continuum; not Clay."
        ),
        "certificate": "CERT-L0026-c0007-techo-N24",
        "notes": (
            f"Ω★={b26.Omega_star:.6f}; L-0024 majorant at Ω★ equals M. "
            f"Low slab remains open (see L-0027 cubic conditional C_†≈{b.C_dagger:.4f})."
        ),
    }
    Path("conjectures/proved_restricted").mkdir(parents=True, exist_ok=True)
    Path("conjectures/proved_restricted/C-R-0002.json").write_text(
        json.dumps(cr, indent=2), encoding="utf-8"
    )

    # Refresh C-0007 notes
    c7_path = Path("conjectures/active/C-0007.json")
    c7 = json.loads(c7_path.read_text(encoding="utf-8"))
    c7["notes"] = (
        f"Open all-IC. High slab Ω0≥Ω★≈{b26.Omega_star:.4f} proved as C-R-0002 "
        f"(L-0024+L-0026). Low slab: L-0027 says C≤C_†≈{b.C_dagger:.4f} in the "
        f"cubic stretch would close (C_emp≈{b.C_emp_max:.4g}); proved C still open."
    )
    c7["related"] = ["C-R-0002", "L-0026", "L-0027"]
    c7_path.write_text(json.dumps(c7, indent=2), encoding="utf-8")

    out = {
        "C_R_0002_M": b26.c0007_M,
        "Omega_star": b26.Omega_star,
        "C_dagger": b.C_dagger,
        "C_emp_max": b.C_emp_max,
        "emp_closes_low_slab": b.emp_closes_low_slab,
        "all_ic_cubic_C0_floor": b.all_ic_cubic_C0_floor,
        "closes_c0007": False,
        "certificate_ok": ok,
        "checks": {k: bool(v) for k, v in checks.items()},
        "not_a_clay_claim": True,
    }
    Path("reports").mkdir(exist_ok=True)
    Path("experiments/exploratory/sprint02_l0027").mkdir(parents=True, exist_ok=True)
    Path("experiments/exploratory/sprint02_l0027/summary.json").write_text(
        json.dumps(out, indent=2), encoding="utf-8"
    )
    Path("reports/SPRINT_02_L0027.md").write_text(
        f"""# Sprint — L-0027 + C-R-0002

| Item | Result |
|------|--------|
| **C-R-0002** (Ω0≥Ω★≈{b26.Omega_star:.4f}, M≈{b26.c0007_M:.6f}) | **proved** (L-0024+L-0026) |
| **L-0027** C_† (low-slab cubic) | ≈{b.C_dagger:.6f} |
| C_emp | ≈{b.C_emp_max:.6g} (≪ C_†) |
| **C-0007** all-IC | **still open** (need proved C≤C_†) |

All-IC cubic at C=0 already floors at ≈{b.all_ic_cubic_C0_floor:.4f}>M; high slab needs defect ODE.

FINITE Galerkin only. Not continuum. Not Clay.
""",
        encoding="utf-8",
    )
    print(json.dumps(out, indent=2, ensure_ascii=True))
    return out


if __name__ == "__main__":
    main()
