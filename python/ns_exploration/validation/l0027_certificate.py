"""Certificate for L-0027 low-slab cubic conditional."""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path

from ns_exploration.conjectures.l0023_cubic_dissipation import omega_ode_bound
from ns_exploration.conjectures.l0026_c0007_techo import lemma_l0026
from ns_exploration.conjectures.l0027_low_slab_cubic import (
    C_star_low_slab,
    all_ic_cubic_floor,
    lemma_l0027,
)
from ns_exploration.validation.l0003_certificate import full_dealias_exact_stats


@dataclass
class CertificateL0027:
    cert_id: str
    lemma_id: str
    n: int
    E0: float
    nu: float
    T: float
    c0007_M: float
    Omega_star: float
    C_dagger: float
    Omega_T_at_Cdagger: float
    all_ic_cubic_C0_floor: float
    proved_closes_c0007: bool
    clay_implication: str

    def as_dict(self) -> dict:
        return asdict(self)


def build_l0027_certificate(
    n: int = 24,
    E0: float = 0.5,
    nu: float = 0.1,
    T: float = 0.02,
) -> CertificateL0027:
    # Skip slow empirical search in cert build
    b = lemma_l0027(n=n, E0=E0, nu=nu, T=T, empirical=False)
    return CertificateL0027(
        cert_id="CERT-L0027-low-slab-cubic-C0007-N24",
        lemma_id="L-0027",
        n=b.n,
        E0=b.E0,
        nu=b.nu,
        T=b.T,
        c0007_M=b.c0007_M,
        Omega_star=b.Omega_star,
        C_dagger=b.C_dagger,
        Omega_T_at_Cdagger=b.Omega_T_at_Cdagger,
        all_ic_cubic_C0_floor=b.all_ic_cubic_C0_floor,
        proved_closes_c0007=b.proved_closes_c0007,
        clay_implication=b.clay_implication,
    )


def verify_l0027_certificate(d: dict) -> tuple[bool, dict]:
    n = int(d["n"])
    E0 = float(d["E0"])
    nu = float(d["nu"])
    T = float(d["T"])
    b26 = lemma_l0026(n=n, E0=E0, nu=nu, T=T, n_grid=41)
    M = b26.c0007_M
    om_star = b26.Omega_star
    Cd, omT = C_star_low_slab(E0, nu, T, M, om_star, n_grid=21)
    K2, _ = full_dealias_exact_stats(n)
    floor0 = all_ic_cubic_floor(E0, nu, T, K2, n_grid=25)
    # Spot-check low slab at C_dagger
    spots = [
        omega_ode_bound(Cd, E0, nu, T, float(o))
        for o in (E0, 0.5 * (E0 + om_star), om_star)
    ]
    checks = {
        "M_matches": abs(M - float(d["c0007_M"])) < 1e-9,
        "Omega_star_matches": abs(om_star - float(d["Omega_star"])) < 1e-5,
        "C_dagger_matches": abs(Cd - float(d["C_dagger"])) < 1e-4,
        "Omega_T_le_M": omT <= M + 1e-9,
        "spots_le_M": max(spots) <= M + 1e-8,
        "all_ic_C0_exceeds_M": floor0 > M,
        "proved_closes_false": d["proved_closes_c0007"] is False,
        "C_dagger_positive": Cd > 1.0,
    }
    return all(checks.values()), checks


def save_l0027_certificate(
    cert: CertificateL0027,
    path: str | Path = "certificates/CERT-L0027-low-slab-cubic-C0007-N24.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(cert.as_dict(), indent=2), encoding="utf-8")
