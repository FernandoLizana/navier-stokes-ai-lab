"""Certificate for L-0030 full-mask triad R_★ / ‖N‖ bound."""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path

from ns_exploration.conjectures.l0020_duhamel_h1 import full_mask_rho_star
from ns_exploration.conjectures.l0030_triad_N_bound import (
    N_triad_bound,
    N_young_bound,
    cubic_C_from_Rstar,
    full_mask_triad_Rstar,
    lemma_l0030,
)


@dataclass
class CertificateL0030:
    cert_id: str
    lemma_id: str
    n: int
    M_modes: int
    R_star: int
    n_triads: int
    rho_star: float
    C_dagger: float
    C_from_Rstar: float
    N_triad_at_equipartition: float
    N_young_at_E0: float
    low_slab_hybrid_worst: float
    c0007_M: float
    closes_c0007: bool
    clay_implication: str

    def as_dict(self) -> dict:
        return asdict(self)


def build_l0030_certificate(n: int = 24) -> CertificateL0030:
    b = lemma_l0030(n=n, n_grid=11)
    return CertificateL0030(
        cert_id="CERT-L0030-triad-N-C0007-N24",
        lemma_id="L-0030",
        n=b.n,
        M_modes=b.M_modes,
        R_star=b.R_star,
        n_triads=b.n_triads,
        rho_star=b.rho_star,
        C_dagger=b.C_dagger,
        C_from_Rstar=b.C_from_Rstar,
        N_triad_at_equipartition=b.N_triad_at_equipartition,
        N_young_at_E0=b.N_young_at_E0,
        low_slab_hybrid_worst=b.low_slab_hybrid_worst,
        c0007_M=b.c0007_M,
        closes_c0007=b.closes_c0007,
        clay_implication=b.clay_implication,
    )


def verify_l0030_certificate(d: dict) -> tuple[bool, dict]:
    n = int(d["n"])
    R_star, n_tri = full_mask_triad_Rstar(n)
    rho, _, _ = full_mask_rho_star(n)
    E0 = 0.5
    checks = {
        "R_star_matches": R_star == int(d["R_star"]),
        "n_triads_matches": n_tri == int(d["n_triads"]),
        "R_star_lt_M": R_star < int(d["M_modes"]),
        "R_star_positive": R_star > 1000,
        "C_from_R_matches": abs(cubic_C_from_Rstar(R_star) - float(d["C_from_Rstar"]))
        < 1e-6,
        "N_triad_beats_young_at_eq": N_triad_bound(E0, E0, R_star)
        < N_young_bound(E0, rho),
        "N_values_match": abs(N_triad_bound(E0, E0, R_star) - float(d["N_triad_at_equipartition"]))
        < 1e-6,
        "hybrid_above_M": float(d["low_slab_hybrid_worst"]) > float(d["c0007_M"]),
        "C_above_Cdagger": float(d["C_from_Rstar"]) > float(d["C_dagger"]),
        "closes_false": d["closes_c0007"] is False,
    }
    return all(checks.values()), checks


def save_l0030_certificate(
    cert: CertificateL0030,
    path: str | Path = "certificates/CERT-L0030-triad-N-C0007-N24.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(cert.as_dict(), indent=2), encoding="utf-8")
