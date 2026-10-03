"""Certificate for L-0034 shell-factorized ‖N‖ bound."""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path

from ns_exploration.conjectures.l0020_duhamel_h1 import full_mask_rho_star
from ns_exploration.conjectures.l0034_shell_N_bound import (
    N_shell_bound,
    lemma_l0034,
    shell_multiplicities,
)


@dataclass
class CertificateL0034:
    cert_id: str
    lemma_id: str
    n: int
    c0007_M: float
    C_dagger: float
    n_shells: int
    N_at_equipartition: float
    N_young: float
    C_eff_at_equipartition: float
    low_slab_ode_worst: float
    closes_c0007: bool
    clay_implication: str

    def as_dict(self) -> dict:
        return asdict(self)


def build_l0034_certificate(n: int = 24) -> CertificateL0034:
    b = lemma_l0034(n=n, n_grid=9)
    return CertificateL0034(
        cert_id="CERT-L0034-shell-N-C0007-N24",
        lemma_id="L-0034",
        n=b.n,
        c0007_M=b.c0007_M,
        C_dagger=b.C_dagger,
        n_shells=b.n_shells,
        N_at_equipartition=b.N_at_equipartition,
        N_young=b.N_young,
        C_eff_at_equipartition=b.C_eff_at_equipartition,
        low_slab_ode_worst=b.low_slab_ode_worst,
        closes_c0007=b.closes_c0007,
        clay_implication=b.clay_implication,
    )


def verify_l0034_certificate(d: dict) -> tuple[bool, dict]:
    n = int(d["n"])
    radii, masses = shell_multiplicities(n)
    rho, _, _ = full_mask_rho_star(n)
    N_eq = N_shell_bound(0.5, 0.5, radii, masses, rho)
    N_y = 2.0 * 0.5 * (rho**0.5)
    checks = {
        "n_shells_match": len(radii) == int(d["n_shells"]),
        "N_eq_matches": abs(N_eq - float(d["N_at_equipartition"])) < 1e-6,
        "beats_young_at_eq": N_eq < N_y,
        "C_eff_above_Cdagger": float(d["C_eff_at_equipartition"]) > float(d["C_dagger"]),
        "ode_above_M": float(d["low_slab_ode_worst"]) > float(d["c0007_M"]),
        "closes_false": d["closes_c0007"] is False,
    }
    return all(checks.values()), checks


def save_l0034_certificate(
    cert: CertificateL0034,
    path: str | Path = "certificates/CERT-L0034-shell-N-C0007-N24.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(cert.as_dict(), indent=2), encoding="utf-8")
