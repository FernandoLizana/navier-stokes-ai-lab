"""Sprint: L-0016 radial Young → longer C-S-0002-SHORT."""

from __future__ import annotations

import json
from pathlib import Path

from ns_exploration.conjectures.l0016_radial import lemma_l0016, save_lemma_l0016
from ns_exploration.validation.l0016_certificate import (
    build_l0016_short_certificate,
    save_l0016_short_certificate,
    verify_l0016_short_certificate,
)


def main() -> dict:
    b = lemma_l0016()
    save_lemma_l0016(b, "conjectures/proved_restricted/L-0016.json")

    short = {
        "id": "C-S-0002-SHORT",
        "route": "B",
        "status": "proved_restricted",
        "evidence_level": "N7",
        "parent": "C-S-0002",
        "lemma": "L-0016",
        "domain": (
            "Pseudospectral T^3, N=16, IC |k|_inf<=4, E=0.5, full dealias. NOT continuum."
        ),
        "statement": (
            f"For all divergence-free Fourier ICs on N=16 with |k|_∞≤4 and energy 0.5, "
            f"evolved with dealiased Galerkin NS (ν={b.nu}), one has Ω(t) ≤ {b.cs0002_M} "
            f"for all t ∈ [0, {b.T_star:.8f}], via radial Young ρ_★={b.rho_star:.0f} "
            f"plus NS cancellation ODE. Does NOT claim T=0.02."
        ),
        "T_max": b.T_star,
        "proposed_bound_M": b.cs0002_M,
        "original_conjecture_T": 0.02,
        "improvement_vs_L0015": b.improvement_vs_L0015,
        "rho_star": b.rho_star,
        "r_star": b.r_star,
        "m_star": b.m_star,
        "N_max": b.N_max,
        "U_eff": b.U_eff,
        "W": b.W,
        "alpha": b.alpha,
        "beta": b.beta,
        "clay_implication": "None. Short-time finite Galerkin only.",
        "notes": b.notes,
    }
    Path("conjectures/proved_restricted/C-S-0002-SHORT.json").write_text(
        json.dumps(short, indent=2), encoding="utf-8"
    )

    cert = build_l0016_short_certificate()
    ok, checks = verify_l0016_short_certificate(cert.as_dict())
    save_l0016_short_certificate(
        cert, "certificates/CERT-L0016-radial-CS0002-short-N16.json"
    )

    out = {
        "T_star_L0015": b.T_L0015_ref,
        "T_star_L0016": b.T_star,
        "improvement_vs_L0015": b.improvement_vs_L0015,
        "ratio_to_0.02": b.T_star / 0.02,
        "rho_star": b.rho_star,
        "U_eff": b.U_eff,
        "U_L0015": b.U_L0015,
        "Omega_at_0.02_comparison": b.Omega_at_0_02,
        "certificate_verified": ok,
        "certificate_T_cert": cert.T_cert,
        "certificate_Omega_hi": cert.Omega_bound_hi,
        "checks": {k: bool(v) for k, v in checks.items()},
        "not_a_clay_claim": True,
    }
    Path("reports").mkdir(exist_ok=True)
    print(json.dumps(out, indent=2, ensure_ascii=True))
    return out


if __name__ == "__main__":
    main()
