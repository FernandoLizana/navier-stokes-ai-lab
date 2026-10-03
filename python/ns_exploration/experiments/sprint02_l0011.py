"""Sprint: L-0011 sharpen C-S-0002-SHORT horizon."""

from __future__ import annotations

import json
from pathlib import Path

from ns_exploration.conjectures.l0011_sharpened_short import lemma_l0011, save_lemma_l0011
from ns_exploration.validation.l0011_certificate import (
    build_l0011_short_certificate,
    save_l0011_short_certificate,
    verify_l0011_short_certificate,
)


def main() -> dict:
    b = lemma_l0011()
    save_lemma_l0011(b, "conjectures/proved_restricted/L-0011.json")

    short = {
        "id": "C-S-0002-SHORT",
        "route": "B",
        "status": "proved_restricted",
        "evidence_level": "N7",
        "parent": "C-S-0002",
        "lemma": "L-0011",
        "domain": (
            "Pseudospectral T^3, N=16, IC |k|_inf<=4, E=0.5, full dealias. NOT continuum."
        ),
        "statement": (
            f"For all divergence-free Fourier ICs on N=16 with |k|_∞≤4 and energy 0.5, "
            f"evolved with dealiased Galerkin NS (ν={b.nu}), one has Ω(t) ≤ {b.cs0002_M} "
            f"for all t ∈ [0, {b.T_star_rigorous:.8f}] (div-free embedding + viscous Duhamel). "
            f"Does NOT claim T=0.02."
        ),
        "T_max": b.T_star_rigorous,
        "T_max_LxL_conditional_N6": b.T_star_LxL_conditional,
        "Omega_LxL_at_0.02_conditional_N6": b.Omega_LxL_at_T_target,
        "proposed_bound_M": b.cs0002_M,
        "original_conjecture_T": 0.02,
        "improvement_vs_L0010": b.improvement_vs_L0010,
        "clay_implication": "None. Short-time finite Galerkin only.",
        "notes": b.notes,
    }
    Path("conjectures/proved_restricted/C-S-0002-SHORT.json").write_text(
        json.dumps(short, indent=2), encoding="utf-8"
    )

    cert = build_l0011_short_certificate()
    ok, checks = verify_l0011_short_certificate(cert.as_dict())
    save_l0011_short_certificate(
        cert, "certificates/CERT-L0011-duhamel-CS0002-short-N16.json"
    )

    out = {
        "T_star_L0010": 0.0010112385074348496,
        "T_star_L0011_rigorous": b.T_star_rigorous,
        "improvement": b.improvement_vs_L0010,
        "T_star_LxL_conditional": b.T_star_LxL_conditional,
        "Omega_LxL_at_0.02_conditional": b.Omega_LxL_at_T_target,
        "certificate_verified": ok,
        "certificate_T_cert": cert.T_cert,
        "checks": {k: bool(v) for k, v in checks.items()},
        "not_a_clay_claim": True,
    }
    Path("reports").mkdir(exist_ok=True)
    print(json.dumps(out, indent=2, ensure_ascii=True))
    return out


if __name__ == "__main__":
    main()
