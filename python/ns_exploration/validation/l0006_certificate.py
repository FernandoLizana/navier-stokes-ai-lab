"""
Interval-arithmetic certificate for the L-0006 algebraic short-time bound,
restricted to the hard Galerkin truncation |k| <= 2 on the N=12 grid.

What is certified (rigorously, N5 for the CONSTANT):
  Given the hand differential inequality (N7)
      dΩ/dt ≤ 2 a Ω^{3/2},   a = sqrt(2)*sqrt(3M),
  with Ω0 ≤ K² E0, the comparison solution gives, while a·√Ω0·t < 1,
      Ω(t) ≤ Ω0 / (1 - a√Ω0 t)²  =: B.
  This module computes an OUTWARD interval enclosure of B using directed
  rounding, and exports B_hi (a guaranteed upper bound) with the exact
  integer data (K²=4, M) needed to re-verify independently.

Scope: shell-truncated Galerkin ODE ONLY (see L-0005). Inviscid upper bound,
so it holds for every ν ≥ 0. NOT continuum NS. NOT Clay. Evidence: N5 (constant)
on top of N7 (inequality).
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from ns_exploration.conjectures.l0005_shell_galerkin import shell_mask
from ns_exploration.spectral.fourier_conventions import k_squared
from ns_exploration.validation.intervals import Interval


def shell_exact_stats(n: int, k_shell: int) -> tuple[int, int]:
    """Return (K_squared_int, M_count) exactly for the Euclidean shell on the dealiased grid."""
    mask = shell_mask(n, float(k_shell), dealias=True)
    k2 = k_squared(n)
    k2_ret = np.rint(k2[mask]).astype(np.int64)
    K2 = int(k2_ret.max()) if k2_ret.size else 0
    M = int(mask.sum())
    return K2, M


@dataclass
class IntervalCertificate:
    certificate_id: str
    lemma_ref: str
    route: str
    evidence_level: str
    n: int
    k_shell: int
    K_squared_int: int
    M: int
    E0: float
    t: float
    a_lo: float
    a_hi: float
    prod_hi: float  # (a * sqrt(Ω0) * t) upper bound; must be < 1
    closes: bool
    Omega_bound_hi: float
    clay_implication: str
    payload_sha256: str

    def as_dict(self) -> dict:
        return asdict(self)


def build_certificate(
    n: int = 12,
    k_shell: int = 2,
    E0: float = 0.5,
    t: float = 0.02,
) -> IntervalCertificate:
    K2, M = shell_exact_stats(n, k_shell)
    IE0 = Interval.from_float(E0)
    It = Interval.from_float(t)
    IK2 = Interval(float(K2), float(K2))  # exact integer
    I3M = Interval(float(3 * M), float(3 * M))
    I2 = Interval(2.0, 2.0)

    Omega0 = IK2 * IE0
    a = (I2 * I3M).sqrt()  # sqrt(2 * 3M) = sqrt(6M) = sqrt(2)*sqrt(3M)
    z0 = Omega0.sqrt()
    prod = a * z0 * It
    closes = bool(prod.hi < 1.0)
    if not closes:
        Omega_bound_hi = float("inf")
    else:
        one = Interval(1.0, 1.0)
        denom = one - prod  # positive interval since prod.hi < 1
        denom_sq = denom.sqr()
        cap = Omega0 / denom_sq
        Omega_bound_hi = float(cap.hi)

    payload = json.dumps(
        {"n": n, "k_shell": k_shell, "K2": K2, "M": M, "E0": E0, "t": t}, sort_keys=True
    ).encode()
    digest = hashlib.sha256(payload).hexdigest()

    return IntervalCertificate(
        certificate_id="CERT-L0006-shell-k2",
        lemma_ref="L-0006 (algebraic short-time), truncated Galerkin |k|<=k_shell (L-0005 support)",
        route="B",
        evidence_level="N5",
        n=n,
        k_shell=k_shell,
        K_squared_int=K2,
        M=M,
        E0=E0,
        t=t,
        a_lo=float(a.lo),
        a_hi=float(a.hi),
        prod_hi=float(prod.hi),
        closes=closes,
        Omega_bound_hi=Omega_bound_hi,
        clay_implication="None. Shell-truncated Galerkin ODE, finite M; not continuum NS.",
        payload_sha256=digest,
    )


def save_certificate(cert: IntervalCertificate, path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(cert.as_dict(), indent=2), encoding="utf-8")
