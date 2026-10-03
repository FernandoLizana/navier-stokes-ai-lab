"""N5-lite certificate wrapper for L-0026 (recomputes case-split constants)."""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path

from ns_exploration.conjectures.l0024_spectral_defect import (
    dealias_shell_stats,
    defect_ode_omega_T,
)
from ns_exploration.conjectures.l0025_multi_n_defect import c0007_open_target
from ns_exploration.conjectures.l0026_c0007_techo import lemma_l0026, omega_star_for_M


@dataclass
class CertificateL0026:
    cert_id: str
    lemma_id: str
    n: int
    E0: float
    nu: float
    T: float
    c0007_M: float
    Stokes_floor: float
    L0024_majorant: float
    Omega_star: float
    Omega_T_at_star: float
    high_slab_closes: bool
    closes_c0007: bool
    monotone_in_Omega0: bool
    clay_implication: str

    def as_dict(self) -> dict:
        return asdict(self)


def build_l0026_certificate(
    n: int = 24,
    E0: float = 0.5,
    nu: float = 0.1,
    T: float = 0.02,
) -> CertificateL0026:
    b = lemma_l0026(n=n, E0=E0, nu=nu, T=T, n_grid=61)
    return CertificateL0026(
        cert_id="CERT-L0026-c0007-techo-N24",
        lemma_id="L-0026",
        n=b.n,
        E0=b.E0,
        nu=b.nu,
        T=b.T,
        c0007_M=b.c0007_M,
        Stokes_floor=b.Stokes_floor,
        L0024_majorant=b.L0024_majorant,
        Omega_star=b.Omega_star,
        Omega_T_at_star=b.Omega_T_at_star,
        high_slab_closes=b.high_slab_closes,
        closes_c0007=b.closes_c0007,
        monotone_in_Omega0=b.monotone_in_Omega0,
        clay_implication=b.clay_implication,
    )


def verify_l0026_certificate(d: dict) -> tuple[bool, dict]:
    n = int(d["n"])
    E0 = float(d["E0"])
    nu = float(d["nu"])
    T = float(d["T"])
    gap = c0007_open_target(n=n, E0=E0, nu=nu, T=T)
    stats = dealias_shell_stats(n)
    M = float(gap["proposed_bound_M"])
    om_star, omT_star = omega_star_for_M(M, E0, nu, T, stats)
    # High-slab spot checks
    K2 = stats["K2"]
    highs = [
        defect_ode_omega_T(E0, float(o), nu, T, stats)
        for o in (om_star, 0.5 * (om_star + K2 * E0), K2 * E0)
    ]
    checks = {
        "M_matches": abs(M - float(d["c0007_M"])) < 1e-9,
        "floor_matches": abs(gap["Stokes_floor"] - float(d["Stokes_floor"])) < 1e-9,
        "majorant_matches": abs(gap["L0024_majorant"] - float(d["L0024_majorant"]))
        < 1e-8,
        "Omega_star_matches": abs(om_star - float(d["Omega_star"])) < 1e-6,
        "Omega_T_at_star_le_M": omT_star <= M + 1e-9,
        "high_slab_le_M": max(highs) <= M + 1e-9,
        "closes_c0007_false": d["closes_c0007"] is False,
        "high_slab_closes_true": d["high_slab_closes"] is True,
        "gap_nontrivial": gap["Stokes_floor"] < M < gap["L0024_majorant"],
    }
    return all(checks.values()), checks


def save_l0026_certificate(
    cert: CertificateL0026,
    path: str | Path = "certificates/CERT-L0026-c0007-techo-N24.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(cert.as_dict(), indent=2), encoding="utf-8")
