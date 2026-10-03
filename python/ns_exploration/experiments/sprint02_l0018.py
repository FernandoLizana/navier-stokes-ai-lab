"""Sprint: L-0018 spectral envelope → prove C-S-0002 on N≤16."""

from __future__ import annotations

import json
from pathlib import Path

from ns_exploration.conjectures.l0018_envelope import lemma_l0018, save_lemma_l0018
from ns_exploration.validation.l0018_certificate import (
    build_l0018_envelope_certificate,
    save_l0018_envelope_certificate,
    verify_l0018_envelope_certificate,
)


def main() -> dict:
    b = lemma_l0018()
    save_lemma_l0018(b, "conjectures/proved_restricted/L-0018.json")

    proved = {
        "id": "C-S-0002",
        "route": "B",
        "status": "proved_restricted",
        "evidence_level": "N7",
        "lemma": "L-0018",
        "domain": (
            "Pseudospectral T^3, N≤16, IC |k|_inf≤4, E=0.5, 2/3 dealias Galerkin. "
            "NOT continuum."
        ),
        "k_shell": 4,
        "nu": b.nu,
        "energy": b.E0,
        "t_end": 0.02,
        "n_max": b.n_max,
        "proposed_bound_M": b.cs0002_M,
        "statement": (
            f"For all divergence-free Fourier ICs on N³ (N≤{b.n_max}) with |k|_∞≤4 "
            f"and kinetic energy E0={b.E0}, evolved under dealiased Galerkin NS (ν={b.nu}), "
            f"one has Ω(t) ≤ {b.cs0002_M} for all t ≥ 0 (hence on [0, 0.02]), via the "
            f"spectral envelope Ω ≤ K² E ≤ K²({b.n_worst}) E0 = {b.Omega_cap} "
            f"(L-0018). FINITE-DIMENSIONAL only."
        ),
        "Omega_cap": b.Omega_cap,
        "K2_worst": b.K2_worst,
        "n_worst": b.n_worst,
        "margin": b.margin,
        "parent": "C-S-0001",
        "clay_implication": b.clay_implication,
        "notes": b.notes,
    }
    Path("conjectures/proved_restricted/C-S-0002.json").write_text(
        json.dumps(proved, indent=2), encoding="utf-8"
    )

    # Remove from active (proved).
    active = Path("conjectures/active/C-S-0002.json")
    if active.exists():
        active.unlink()

    cert = build_l0018_envelope_certificate()
    ok, checks = verify_l0018_envelope_certificate(cert.as_dict())
    save_l0018_envelope_certificate(
        cert, "certificates/CERT-L0018-envelope-CS0002-N16.json"
    )

    out = {
        "closes_cs0002": b.closes_cs0002,
        "Omega_cap": b.Omega_cap,
        "K2_worst": b.K2_worst,
        "margin": b.margin,
        "certificate_verified": ok,
        "checks": {k: bool(v) for k, v in checks.items()},
        "not_a_clay_claim": True,
        "does_not_close_C0002": True,
    }
    Path("reports").mkdir(exist_ok=True)
    print(json.dumps(out, indent=2, ensure_ascii=True))
    return out


if __name__ == "__main__":
    main()
