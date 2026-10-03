"""
C-0004 / C-0005: all-IC enstrophy line on N≤24.

C-0004: M = K²(24)·E0 = 73.5 — proved by L-0018 envelope (all t≥0).
C-0005: M = 1.5 · Stokes_floor(T) ≈ 61.24 — proved by L-0024 spectral-defect ODE.
  Stokes floor at T=0.02: Ω = 147·E0·exp(-2ν·147·T) ≈ 40.825 (exact single-mode).
  Any all-IC claim with M < Stokes floor is false; C-0005 sits above that floor
  and below the envelope.

FINITE Galerkin / dealias only. Not continuum. Not Clay.
"""

from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

from ns_exploration.conjectures.c0002_stokes_refuter import stokes_Omega_exact
from ns_exploration.conjectures.l0018_envelope import lemma_l0018, save_lemma_l0018
from ns_exploration.validation.l0003_certificate import full_dealias_exact_stats
from ns_exploration.validation.l0018_certificate import (
    build_l0018_envelope_certificate,
    save_l0018_envelope_certificate,
    verify_l0018_envelope_certificate,
)


@dataclass
class ConjectureC0005:
    id: str = "C-0005"
    route: str = "B"
    status: str = "exploring"
    evidence_level: str = "N6"
    domain: str = "Galerkin/pseudospectral T^3, N<=24, 2/3 dealias (NOT continuum PDE)"
    nu: float = 0.1
    energy: float = 0.5
    t_end: float = 0.02
    n_max: int = 24
    proposed_bound_M: float = 0.0
    parent: str = "C-0004"
    stokes_floor: float = 0.0
    envelope_cap: float = 0.0
    statement: str = ""
    clay_implication: str = (
        "None. Finite-dimensional only; no continuum regularity or blow-up claim."
    )
    numerical_support_max: float | None = None
    refuted: bool = False
    refutation_value: float | None = None
    notes: str = ""

    def as_dict(self) -> dict:
        return asdict(self)


def stokes_floor_N(n: int, E0: float, nu: float, t: float) -> tuple[float, int]:
    K2, _ = full_dealias_exact_stats(n)
    return stokes_Omega_exact(K2, E0, nu, t), int(K2)


def build_c0004_and_c0005(
    n_max: int = 24,
    E0: float = 0.5,
    nu: float = 0.1,
    t_end: float = 0.02,
    cushion: float = 1.5,
) -> dict:
    b18 = lemma_l0018(n_max=n_max, E0=E0, nu=nu, cs0002_M=1e100)
    M4 = b18.Omega_cap
    floor, K2 = stokes_floor_N(n_max, E0, nu, t_end)
    M5 = float(cushion * floor)
    assert M5 < M4, "C-0005 must sit strictly below the envelope to be nontrivial"
    assert M5 > floor

    c4 = {
        "id": "C-0004",
        "route": "B",
        "status": "proved_restricted",
        "evidence_level": "N7",
        "lemma": "L-0018",
        "domain": f"Galerkin/pseudospectral T^3, N≤{n_max}, 2/3 dealias (NOT continuum PDE)",
        "nu": nu,
        "energy": E0,
        "t_end": t_end,
        "n_max": n_max,
        "proposed_bound_M": M4,
        "parent": "C-0003",
        "statement": (
            f"For all divergence-free Fourier fields on the N³ grid with N≤{n_max}, "
            f"dealiased Galerkin NS (ν={nu}), kinetic energy ≤ {E0}, one has "
            f"Ω(t) ≤ {M4} for all t ≥ 0, via Ω ≤ K² E ≤ K²({b18.n_worst}) E0 = {M4} "
            f"(L-0018). FINITE-DIMENSIONAL only."
        ),
        "Omega_cap": M4,
        "K2_worst": b18.K2_worst,
        "clay_implication": b18.clay_implication,
        "notes": (
            f"Envelope extension of C-0003 to N≤{n_max}. "
            f"Stokes floor at T={t_end} is ≈{floor:.6f} (does not threaten M={M4})."
        ),
    }

    Path("conjectures/proved_restricted").mkdir(parents=True, exist_ok=True)
    Path("conjectures/proved_restricted/C-0004.json").write_text(
        json.dumps(c4, indent=2), encoding="utf-8"
    )
    # C-0005 is proved by L-0024; do not recreate an active exploring record.
    active_c5 = Path("conjectures/active/C-0005.json")
    if active_c5.exists():
        active_c5.unlink()

    # L-0018 artifact for N=24
    b24 = lemma_l0018(n_max=n_max, E0=E0, nu=nu, cs0002_M=M4)
    b24.notes = (
        f"Ω ≤ K² E0 with K²({b24.n_worst})={b24.K2_worst} ⇒ Ω_cap={b24.Omega_cap}. "
        f"Closes C-0004 (M={M4}) on N≤{n_max}."
    )
    save_lemma_l0018(b24, "conjectures/proved_restricted/L-0018-N24.json")

    cert = build_l0018_envelope_certificate(
        n_max=n_max,
        E0=E0,
        R=M4,
        certificate_id="CERT-L0018-envelope-C0004-N24",
    )
    ok, checks = verify_l0018_envelope_certificate(cert.as_dict())
    save_l0018_envelope_certificate(
        cert, "certificates/CERT-L0018-envelope-C0004-N24.json"
    )

    return {
        "C0004_M": M4,
        "C0004_proved": ok and cert.closes,
        "C0005_M": M5,
        "stokes_floor": floor,
        "K2": K2,
        "gap_envelope_minus_M5": M4 - M5,
        "checks": {k: bool(v) for k, v in checks.items()},
        "not_a_clay_claim": True,
    }
