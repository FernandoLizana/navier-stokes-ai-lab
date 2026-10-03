"""Sprint: L-0014 triad multiplicity → longer C-S-0002-SHORT."""

from __future__ import annotations

import json
from pathlib import Path

from ns_exploration.conjectures.l0014_triad import lemma_l0014, save_lemma_l0014
from ns_exploration.validation.l0014_certificate import (
    build_l0014_short_certificate,
    save_l0014_short_certificate,
    verify_l0014_short_certificate,
)


def main() -> dict:
    b = lemma_l0014()
    save_lemma_l0014(b, "conjectures/proved_restricted/L-0014.json")

    short = {
        "id": "C-S-0002-SHORT",
        "route": "B",
        "status": "proved_restricted",
        "evidence_level": "N7",
        "parent": "C-S-0002",
        "lemma": "L-0014",
        "domain": (
            "Pseudospectral T^3, N=16, IC |k|_inf<=4, E=0.5, full dealias. NOT continuum."
        ),
        "statement": (
            f"For all divergence-free Fourier ICs on N=16 with |k|_∞≤4 and energy 0.5, "
            f"evolved with dealiased Galerkin NS (ν={b.nu}), one has Ω(t) ≤ {b.cs0002_M} "
            f"for all t ∈ [0, {b.T_star:.8f}], via shell→high triad multiplicity R_H={b.R_H} "
            f"plus NS cancellation ODE. Does NOT claim T=0.02."
        ),
        "T_max": b.T_star,
        "proposed_bound_M": b.cs0002_M,
        "original_conjecture_T": 0.02,
        "improvement_vs_L0013": b.improvement_vs_L0013,
        "R_H": b.R_H,
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

    cert = build_l0014_short_certificate()
    ok, checks = verify_l0014_short_certificate(cert.as_dict())
    save_l0014_short_certificate(
        cert, "certificates/CERT-L0014-triad-CS0002-short-N16.json"
    )

    out = {
        "T_star_L0013": b.T_L0013_ref,
        "T_star_L0014": b.T_star,
        "improvement_vs_L0013": b.improvement_vs_L0013,
        "ratio_to_0.02": b.T_star / 0.02,
        "R_H": b.R_H,
        "U_eff": b.U_eff,
        "U_L_L0013": b.U_L_L0013,
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
