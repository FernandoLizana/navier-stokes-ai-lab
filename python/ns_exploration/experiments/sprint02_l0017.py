"""Sprint: L-0017 radial H×L cross → longer C-S-0002-SHORT."""

from __future__ import annotations

import json
from pathlib import Path

from ns_exploration.conjectures.l0017_cross import lemma_l0017, save_lemma_l0017
from ns_exploration.validation.l0017_certificate import (
    build_l0017_short_certificate,
    save_l0017_short_certificate,
    verify_l0017_short_certificate,
)


def main() -> dict:
    b = lemma_l0017()
    save_lemma_l0017(b, "conjectures/proved_restricted/L-0017.json")

    short = {
        "id": "C-S-0002-SHORT",
        "route": "B",
        "status": "proved_restricted",
        "evidence_level": "N7",
        "parent": "C-S-0002",
        "lemma": "L-0017",
        "domain": (
            "Pseudospectral T^3, N=16, IC |k|_inf<=4, E=0.5, full dealias. NOT continuum."
        ),
        "statement": (
            f"For all divergence-free Fourier ICs on N=16 with |k|_∞≤4 and energy 0.5, "
            f"evolved with dealiased Galerkin NS (ν={b.nu}), one has Ω(t) ≤ {b.cs0002_M} "
            f"for all t ∈ [0, {b.T_star:.8f}], via radial Young bounds on L×L and H×L "
            f"(ρ_★={b.rho_star:.0f}, γ_cross={b.gamma_cross:.4g}) plus NS cancellation ODE. "
            f"Does NOT claim T=0.02."
        ),
        "T_max": b.T_star,
        "proposed_bound_M": b.cs0002_M,
        "original_conjecture_T": 0.02,
        "improvement_vs_L0016": b.improvement_vs_L0016,
        "rho_star": b.rho_star,
        "U_eff": b.U_eff,
        "gamma_cross": b.gamma_cross,
        "W_sqrt_B": b.W_sqrt_B,
        "alpha": b.alpha,
        "beta": b.beta,
        "clay_implication": "None. Short-time finite Galerkin only.",
        "notes": b.notes,
    }
    Path("conjectures/proved_restricted/C-S-0002-SHORT.json").write_text(
        json.dumps(short, indent=2), encoding="utf-8"
    )

    cert = build_l0017_short_certificate()
    ok, checks = verify_l0017_short_certificate(cert.as_dict())
    save_l0017_short_certificate(
        cert, "certificates/CERT-L0017-cross-CS0002-short-N16.json"
    )

    out = {
        "T_star_L0016": b.T_L0016_ref,
        "T_star_L0017": b.T_star,
        "improvement_vs_L0016": b.improvement_vs_L0016,
        "ratio_to_0.02": b.T_star / 0.02,
        "gamma_cross": b.gamma_cross,
        "W_sqrt_B": b.W_sqrt_B,
        "U_eff": b.U_eff,
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
