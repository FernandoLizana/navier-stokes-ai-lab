"""Sprint: L-0025 @ N=32 proves C-0006; open C-0007 (N=24 tighter gap)."""

from __future__ import annotations

import json
from pathlib import Path

from ns_exploration.conjectures.l0025_multi_n_defect import (
    c0007_open_target,
    lemma_l0025,
    save_lemma_l0025,
)
from ns_exploration.validation.l0025_certificate import (
    build_l0025_certificate,
    save_l0025_certificate,
    verify_l0025_certificate,
)


def main() -> dict:
    b = lemma_l0025(n=32)
    save_lemma_l0025(b)
    cert = build_l0025_certificate(n=32)
    ok, checks = verify_l0025_certificate(cert.as_dict())
    save_l0025_certificate(cert)
    if not b.closes_c0006:
        raise RuntimeError("L-0025 did not close C-0006")

    c6 = {
        "id": "C-0006",
        "route": "B",
        "status": "proved_restricted",
        "evidence_level": "N7",
        "lemma": "L-0025",
        "domain": "Galerkin/pseudospectral T^3, N<=32, 2/3 dealias (NOT continuum PDE)",
        "nu": b.nu,
        "energy": b.E0,
        "t_end": b.T,
        "n_max": b.n,
        "proposed_bound_M": b.c0006_M,
        "proved_bound_Omega_T": b.Omega_T_worst,
        "parent": "C-0005",
        "stokes_floor": b.Stokes_floor,
        "envelope_cap": b.envelope_cap,
        "statement": (
            f"For all divergence-free Fourier fields on N≤{b.n} with energy ≤{b.E0}, "
            f"dealiased Galerkin NS (ν={b.nu}), Ω({b.T}) ≤ {b.c0006_M:.8f}. "
            f"Proved via L-0025/L-0024 spectral-defect ODE "
            f"(majorant Ω(T)≤{b.Omega_T_worst:.8f}). FINITE-DIMENSIONAL only."
        ),
        "clay_implication": b.clay_implication,
        "K2": b.K2,
        "gap": b.gap,
        "m_K": b.m_K,
        "certificate": "CERT-L0025-spectral-defect-C0006-N32",
        "notes": b.notes,
    }
    Path("conjectures/proved_restricted").mkdir(parents=True, exist_ok=True)
    Path("conjectures/proved_restricted/C-0006.json").write_text(
        json.dumps(c6, indent=2), encoding="utf-8"
    )

    g7 = c0007_open_target(n=24)
    c7 = {
        "id": "C-0007",
        "route": "B",
        "status": "exploring",
        "evidence_level": "N6",
        "domain": "Galerkin/pseudospectral T^3, N<=24, 2/3 dealias (NOT continuum PDE)",
        "nu": 0.1,
        "energy": 0.5,
        "t_end": 0.02,
        "n_max": 24,
        "proposed_bound_M": g7["proposed_bound_M"],
        "parent": "C-0005",
        "stokes_floor": g7["Stokes_floor"],
        "L0024_majorant": g7["L0024_majorant"],
        "envelope_cap": g7["envelope_cap"],
        "statement": (
            f"For all divergence-free Fourier fields on N≤24 with energy ≤0.5, "
            f"dealiased Galerkin NS (ν=0.1), Ω(0.02) ≤ {g7['proposed_bound_M']:.8f}. "
            f"Target sits in the open gap (Stokes floor {g7['Stokes_floor']:.6f}, "
            f"L-0024 majorant {g7['L0024_majorant']:.6f}). FINITE-DIMENSIONAL only."
        ),
        "clay_implication": "None. Finite-dimensional only; no continuum claim.",
        "refuted": False,
        "notes": (
            f"Open: M≈{g7['proposed_bound_M']:.6f} is below L-0024 majorant "
            f"{g7['L0024_majorant']:.6f}; needs a sharper defect/stretch estimate."
        ),
    }
    Path("conjectures/active").mkdir(parents=True, exist_ok=True)
    Path("conjectures/active/C-0007.json").write_text(
        json.dumps(c7, indent=2), encoding="utf-8"
    )

    out = {
        "C0006_M": b.c0006_M,
        "C0006_Omega_T_worst": b.Omega_T_worst,
        "C0006_proved": True,
        "C0007_M": g7["proposed_bound_M"],
        "C0007_open": True,
        "certificate_ok": ok,
        "checks": {k: bool(v) for k, v in checks.items()},
        "not_a_clay_claim": True,
    }
    Path("reports").mkdir(exist_ok=True)
    Path("experiments/exploratory/sprint02_l0025").mkdir(parents=True, exist_ok=True)
    Path("experiments/exploratory/sprint02_l0025/summary.json").write_text(
        json.dumps(out, indent=2), encoding="utf-8"
    )
    print(json.dumps(out, indent=2, ensure_ascii=True))
    return out


if __name__ == "__main__":
    main()
