"""
Interval-arithmetic certificates for the L-0003 uniform-in-time Galerkin bound.

Supports two finite supports on the N³ grid:
  - shell:        Euclidean |k| ≤ k_shell (hard truncation; L-0005)
  - full_dealias: 2/3-rule dealias mask (standard pseudospectral Galerkin)

What is certified (N5 for the CONSTANT):
  Given the hand differential inequality (N7, L-0003)
      dΩ/dt ≤ -ν (2 Ω² / E0) + C Ω^{3/2},
  the comparison ODE yields, for all t ≥ 0,
      Ω(t) ≤ max( K² E0,  6 M E0² / ν² ) =: B.
  This module computes an OUTWARD interval enclosure of B and exports B_hi.

Scope: finite Galerkin ODE ONLY. Fixed N,M — not uniform as N→∞.
NOT continuum NS. NOT Clay. Evidence: N5 (constant) on top of N7 (inequality).
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from ns_exploration.spectral.dealias import dealias_mask
from ns_exploration.spectral.fourier_conventions import k_squared
from ns_exploration.validation.intervals import Interval
from ns_exploration.validation.l0006_certificate import shell_exact_stats


def full_dealias_exact_stats(n: int) -> tuple[int, int]:
    """Return (K_squared_int, M_count) for the 2/3 dealias mask."""
    mask = dealias_mask(n)
    k2 = k_squared(n)
    k2_ret = np.rint(k2[mask]).astype(np.int64)
    K2 = int(k2_ret.max()) if k2_ret.size else 0
    M = int(mask.sum())
    return K2, M


@dataclass
class UniformIntervalCertificate:
    certificate_id: str
    lemma_ref: str
    route: str
    evidence_level: str
    support: str  # "shell" | "full_dealias"
    n: int
    k_shell: int | None
    K_squared_int: int
    M: int
    E0: float
    nu: float
    Omega0_hi: float
    Omega_eq_hi: float
    Omega_bound_hi: float
    closes: bool
    clay_implication: str
    payload_sha256: str
    notes: str = ""

    def as_dict(self) -> dict:
        return asdict(self)


def _exact_stats(support: str, n: int, k_shell: int | None) -> tuple[int, int]:
    if support == "shell":
        if k_shell is None:
            raise ValueError("k_shell required for support='shell'")
        return shell_exact_stats(n, int(k_shell))
    if support == "full_dealias":
        return full_dealias_exact_stats(n)
    raise ValueError(f"unknown support: {support}")


def build_l0003_certificate(
    n: int = 12,
    k_shell: int | None = 3,
    E0: float = 0.5,
    nu: float = 0.1,
    support: str = "shell",
) -> UniformIntervalCertificate:
    if nu <= 0.0 or E0 <= 0.0:
        raise ValueError("require nu>0, E0>0")
    if support == "full_dealias":
        k_shell = None

    K2, M = _exact_stats(support, n, k_shell)
    IE0 = Interval.from_float(E0)
    Inu = Interval.from_float(nu)
    IK2 = Interval(float(K2), float(K2))
    IM = Interval(float(M), float(M))
    I6 = Interval(6.0, 6.0)

    Omega0 = IK2 * IE0
    Omega_eq = (I6 * IM * IE0.sqr()) / Inu.sqr()
    Omega_bound_hi = max(float(Omega0.hi), float(Omega_eq.hi))

    payload = json.dumps(
        {
            "n": n,
            "support": support,
            "k_shell": k_shell,
            "K2": K2,
            "M": M,
            "E0": E0,
            "nu": nu,
        },
        sort_keys=True,
    ).encode()
    digest = hashlib.sha256(payload).hexdigest()

    if support == "shell":
        cert_id = f"CERT-L0003-shell-k{k_shell}"
        support_note = f"hard shell |k|≤{k_shell}"
    else:
        cert_id = f"CERT-L0003-full-dealias-N{n}"
        support_note = f"full 2/3 dealias mask on N={n}"

    return UniformIntervalCertificate(
        certificate_id=cert_id,
        lemma_ref=(
            f"L-0003 (uniform-in-time), truncated Galerkin ({support_note})"
        ),
        route="B",
        evidence_level="N5",
        support=support,
        n=n,
        k_shell=k_shell,
        K_squared_int=K2,
        M=M,
        E0=E0,
        nu=nu,
        Omega0_hi=float(Omega0.hi),
        Omega_eq_hi=float(Omega_eq.hi),
        Omega_bound_hi=Omega_bound_hi,
        closes=True,
        clay_implication=(
            "None. Uniform in time for fixed finite support only; "
            "Ω_eq ∼ M/ν² → ∞ as resolution → ∞. Not continuum NS."
        ),
        payload_sha256=digest,
        notes=(
            f"Outward enclosure: Ω(t)≤max(K²E0, 6ME0²/ν²)≤{Omega_bound_hi:.6e} "
            f"for all t≥0 under {support_note}. "
            f"Ω0_hi={float(Omega0.hi):.6e}, Ω_eq_hi={float(Omega_eq.hi):.6e}."
        ),
    )


def build_l0003_full_dealias_certificate(
    n: int = 12,
    E0: float = 0.5,
    nu: float = 0.1,
) -> UniformIntervalCertificate:
    return build_l0003_certificate(n=n, k_shell=None, E0=E0, nu=nu, support="full_dealias")


def save_l0003_certificate(cert: UniformIntervalCertificate, path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(cert.as_dict(), indent=2), encoding="utf-8")
