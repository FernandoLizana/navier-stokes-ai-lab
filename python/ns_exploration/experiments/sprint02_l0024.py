"""Sprint: L-0024 spectral defect → prove C-0005."""

from __future__ import annotations

import json
from pathlib import Path

from ns_exploration.conjectures.l0024_spectral_defect import lemma_l0024, save_lemma_l0024
from ns_exploration.validation.l0024_certificate import (
    build_l0024_certificate,
    save_l0024_certificate,
    verify_l0024_certificate,
)


def main() -> dict:
    b = lemma_l0024()
    save_lemma_l0024(b)
    cert = build_l0024_certificate()
    ok, checks = verify_l0024_certificate(cert.as_dict())
    save_l0024_certificate(cert)

    if not b.closes_c0005:
        raise RuntimeError("L-0024 did not close C-0005")

    # Move / write proved C-0005
    active = Path("conjectures/active/C-0005.json")
    proved_dir = Path("conjectures/proved_restricted")
    proved_dir.mkdir(parents=True, exist_ok=True)
    c5 = {
        "id": "C-0005",
        "route": "B",
        "status": "proved_restricted",
        "evidence_level": "N7",
        "lemma": "L-0024",
        "domain": (
            "Galerkin/pseudospectral T^3, N<=24, 2/3 dealias (NOT continuum PDE)"
        ),
        "nu": b.nu,
        "energy": b.E0,
        "t_end": b.T,
        "n_max": b.n,
        "proposed_bound_M": b.c0005_M,
        "proved_bound_Omega_T": b.Omega_T_worst,
        "parent": "C-0004",
        "stokes_floor": 40.82462307930434,
        "envelope_cap": 73.5,
        "statement": (
            f"For all divergence-free Fourier fields on N≤{b.n} with energy ≤{b.E0}, "
            f"dealiased Galerkin NS (ν={b.nu}), Ω({b.T}) ≤ {b.c0005_M:.8f}. "
            f"Proved via L-0024 spectral-defect comparison ODE "
            f"(majorant Ω(T)≤{b.Omega_T_worst:.8f}). FINITE-DIMENSIONAL only."
        ),
        "clay_implication": b.clay_implication,
        "K2": b.K2,
        "gap": b.gap,
        "m_K": b.m_K,
        "Omega0_worst": b.Omega0_worst,
        "certificate": "CERT-L0024-spectral-defect-C0005-N24",
        "notes": b.notes,
    }
    Path(proved_dir / "C-0005.json").write_text(json.dumps(c5, indent=2), encoding="utf-8")
    if active.exists():
        active.unlink()

    out = {
        "closes_c0005": True,
        "Omega_T_worst": b.Omega_T_worst,
        "Omega0_worst": b.Omega0_worst,
        "M": b.c0005_M,
        "K2": b.K2,
        "gap": b.gap,
        "m_K": b.m_K,
        "certificate_ok": ok,
        "checks": {k: bool(v) for k, v in checks.items()},
        "not_a_clay_claim": True,
    }
    Path("reports").mkdir(exist_ok=True)
    Path("experiments/exploratory/sprint02_l0024").mkdir(parents=True, exist_ok=True)
    Path("experiments/exploratory/sprint02_l0024/summary.json").write_text(
        json.dumps(out, indent=2), encoding="utf-8"
    )
    print(json.dumps(out, indent=2, ensure_ascii=True))
    return out


if __name__ == "__main__":
    main()
