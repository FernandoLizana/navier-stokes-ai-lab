"""Sprint: L-0010 short-time Duhamel proof of C-S-0002 bound on [0,T*]."""

from __future__ import annotations

import json
from pathlib import Path

from ns_exploration.conjectures.l0010_duhamel_short import lemma_l0010, save_lemma_l0010
from ns_exploration.validation.l0010_certificate import (
    build_l0010_short_certificate,
    save_l0010_short_certificate,
    verify_l0010_short_certificate,
)


def main() -> dict:
    b = lemma_l0010()
    save_lemma_l0010(b, "conjectures/proved_restricted/L-0010.json")

    # Record short-time proved corollary of C-S-0002
    short = {
        "id": "C-S-0002-SHORT",
        "route": "B",
        "status": "proved_restricted",
        "evidence_level": "N7",
        "parent": "C-S-0002",
        "domain": (
            "Pseudospectral T^3, N=16, IC |k|_inf<=4, E=0.5, full dealias evolution. "
            "NOT continuum PDE."
        ),
        "statement": (
            f"For all divergence-free Fourier ICs on N=16 with |k|_∞≤4 and energy 0.5, "
            f"evolved with the dealiased Galerkin NS (ν={b.nu}), one has Ω(t) ≤ {b.cs0002_M} "
            f"for all t ∈ [0, {b.T_duhamel_for_cs0002:.8f}]. "
            f"Does NOT claim the original C-S-0002 horizon T=0.02."
        ),
        "T_max": b.T_duhamel_for_cs0002,
        "proposed_bound_M": b.cs0002_M,
        "original_conjecture_T": 0.02,
        "clay_implication": "None. Short-time finite Galerkin only.",
        "notes": b.notes,
    }
    Path("conjectures/proved_restricted").mkdir(parents=True, exist_ok=True)
    Path("conjectures/proved_restricted/C-S-0002-SHORT.json").write_text(
        json.dumps(short, indent=2), encoding="utf-8"
    )

    cert = build_l0010_short_certificate()
    ok, checks = verify_l0010_short_certificate(cert.as_dict())
    save_l0010_short_certificate(
        cert, "certificates/CERT-L0010-duhamel-CS0002-short-N16.json"
    )

    out = {
        "L0010": {
            "T_duhamel_cs0002": b.T_duhamel_for_cs0002,
            "T_twoscale_max": b.T_twoscale_max,
            "closes_at_0.02": b.duhamel_closes_at_T_target,
            "ratio_Tstar_over_0.02": b.T_duhamel_for_cs0002 / 0.02,
        },
        "C_S_0002_SHORT": short,
        "certificate_verified": ok,
        "certificate_T_cert": cert.T_cert,
        "certificate_checks": {k: bool(v) for k, v in checks.items()},
        "not_a_clay_claim": True,
    }
    Path("reports").mkdir(exist_ok=True)
    print(json.dumps(out, indent=2, ensure_ascii=True))
    return out


if __name__ == "__main__":
    main()
