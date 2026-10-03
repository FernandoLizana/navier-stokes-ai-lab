"""Certificate for L-0031/32/33 techos + α=0 mono-radial subclass."""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path

from ns_exploration.conjectures.l0020_duhamel_h1 import full_mask_rho_star, lemma_l0020, omega_duhamel_H1
from ns_exploration.conjectures.l0031_0033_techos import (
    absolute_cubic_frobenius,
    full_mask_Smax,
    lemma_l0031_33,
    zero_alpha_shells_from_l0021,
)


@dataclass
class CertificateL0031_33:
    cert_id: str
    lemma_id: str
    n: int
    c0007_M: float
    C_dagger: float
    compatible_worst_OmT: float
    L0024_majorant: float
    S_max_full: float
    N_from_Smax: float
    N_young: float
    cubic_F: float
    C_from_cubic_F: float
    zero_alpha_shells: list[int]
    zero_alpha_duhamel_OmT: float
    closes_c0007: bool
    clay_implication: str

    def as_dict(self) -> dict:
        return asdict(self)


def build_l0031_33_certificate(n: int = 24) -> CertificateL0031_33:
    b = lemma_l0031_33(n=n, n_grid=21)
    return CertificateL0031_33(
        cert_id="CERT-L0031-0033-techos-C0007-N24",
        lemma_id="L-0031+L-0032+L-0033",
        n=b.n,
        c0007_M=b.c0007_M,
        C_dagger=b.C_dagger,
        compatible_worst_OmT=b.compatible_worst_OmT,
        L0024_majorant=b.L0024_majorant,
        S_max_full=b.S_max_full,
        N_from_Smax=b.N_from_Smax,
        N_young=b.N_young,
        cubic_F=b.cubic_F,
        C_from_cubic_F=b.C_from_cubic_F,
        zero_alpha_shells=list(b.zero_alpha_shells or []),
        zero_alpha_duhamel_OmT=b.zero_alpha_duhamel_OmT,
        closes_c0007=b.closes_c0007,
        clay_implication=b.clay_implication,
    )


def verify_l0031_33_certificate(d: dict) -> tuple[bool, dict]:
    n = int(d["n"])
    Smax = full_mask_Smax(n)
    rho, _, _ = full_mask_rho_star(n)
    F, _ = absolute_cubic_frobenius(n)
    zshells = zero_alpha_shells_from_l0021()
    b20 = lemma_l0020(n=n)
    z_om = omega_duhamel_H1(b20.S_T, 0.5, 0.0, b20.I_sigma)
    checks = {
        "compatible_above_M": float(d["compatible_worst_OmT"]) > float(d["c0007_M"]),
        "compatible_above_L0024": float(d["compatible_worst_OmT"])
        >= float(d["L0024_majorant"]) - 1e-3,
        "Smax_matches": abs(Smax - float(d["S_max_full"])) < 1e-6,
        "Smax_worse_than_young": 2.0 * 0.5 * (Smax**0.5) > 2.0 * 0.5 * (rho**0.5),
        "F_matches": abs(F - float(d["cubic_F"])) < 1e-2,
        "C_cubic_above_Cdagger": float(d["C_from_cubic_F"]) > float(d["C_dagger"]),
        "zero_alpha_nonempty": len(zshells) >= 3,
        "zero_alpha_le_M": z_om <= float(d["c0007_M"]) + 1e-9,
        "closes_false": d["closes_c0007"] is False,
    }
    return all(checks.values()), checks


def save_l0031_33_certificate(
    cert: CertificateL0031_33,
    path: str | Path = "certificates/CERT-L0031-0033-techos-C0007-N24.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(cert.as_dict(), indent=2), encoding="utf-8")
