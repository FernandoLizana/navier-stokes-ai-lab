"""
N5 certificate for L-0009 two-scale bound on N=12, euclidean k≤2:
  Ω(t) ≤ B cosh²(u*),  F(u*) = F(0) + K_full t.

Finds u_hi by float search, then verifies with outward interval arithmetic that
F(u_hi) - F(0) ≥ K t (so u* ≤ u_hi) and exports B_hi · cosh(u_hi)_hi².
"""

from __future__ import annotations

import hashlib
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from ns_exploration.conjectures.l0009_twoscale import lemma_l0009, shell_mode_count
from ns_exploration.validation.intervals import Interval
from ns_exploration.validation.l0003_certificate import full_dealias_exact_stats
from ns_exploration.conjectures.l0007_shell_ic_fullgrid import euclidean_shell_K2


def _tanh_outward(x: Interval) -> Interval:
    return Interval(
        np.nextafter(math.tanh(x.lo), -np.inf),
        np.nextafter(math.tanh(x.hi), np.inf),
    )


def _log_abs_outward(x: Interval) -> Interval:
    """Outward log|·| on an interval that does not contain 0."""
    if x.lo <= 0.0 <= x.hi:
        raise ValueError("log_abs interval contains 0")
    # |x| on interval away from 0
    if x.lo > 0:
        abs_iv = x
    else:
        abs_iv = Interval(-x.hi, -x.lo)
    return Interval(
        np.nextafter(math.log(abs_iv.lo), -np.inf),
        np.nextafter(math.log(abs_iv.hi), np.inf),
    )


def _F_interval(u: float, U_L: Interval, lam: Interval) -> Interval:
    """Interval enclosure of F(u) for fixed float u ≥ 0."""
    Iu = Interval.from_float(u)
    half = Iu * Interval.from_float(0.5)
    th = _tanh_outward(half)
    s = (U_L.sqr() + lam.sqr()).sqrt()
    # num = th - lam/U_L + s/U_L ; den = th - lam/U_L - s/U_L
    ratio_l = lam / U_L
    ratio_s = s / U_L
    base = th - ratio_l
    num = base + ratio_s
    den = base - ratio_s
    return _log_abs_outward(num / den) / s


@dataclass
class L0009Certificate:
    certificate_id: str
    lemma_ref: str
    route: str
    evidence_level: str
    n: int
    ic_kind: str
    k_ic: int
    K_ic_squared: int
    K_full_squared: int
    M_L: int
    M_H: int
    E0: float
    t: float
    B_hi: float
    U_L_hi: float
    W_hi: float
    u_hi: float
    F_gap_lo: float  # verified F(u_hi)-F(0) lower bound
    Kt_hi: float
    Omega_bound_hi: float
    closes: bool
    clay_implication: str
    payload_sha256: str
    notes: str = ""

    def as_dict(self) -> dict:
        return asdict(self)


def build_l0009_certificate(
    n: int = 12,
    k_ic: int = 2,
    ic_kind: str = "euclidean",
    E0: float = 0.5,
    t: float = 0.02,
) -> L0009Certificate:
    if ic_kind != "euclidean":
        raise ValueError("N5 certificate currently implemented for euclidean shells")

    K_ic2 = euclidean_shell_K2(n, k_ic)
    K2_full, M_full = full_dealias_exact_stats(n)
    M_L = shell_mode_count(n, k_ic, ic_kind)
    M_H = M_full - M_L
    K_full = math.sqrt(float(K2_full))

    IE0 = Interval.from_float(E0)
    It = Interval.from_float(t)
    IB = Interval(float(K_ic2), float(K_ic2)) * IE0
    IK = Interval(float(K2_full), float(K2_full)).sqrt()
    IM_L = Interval(float(M_L), float(M_L))
    IM_H = Interval(float(M_H), float(M_H))
    I3 = Interval(3.0, 3.0)
    I2 = Interval(2.0, 2.0)

    U_L = (I3 * IM_L).sqrt() * (I2 * IE0).sqrt()
    W = (I3 * IM_H).sqrt() * I2.sqrt() if M_H > 0 else Interval(0.0, 0.0)
    lam = (W * IB.sqrt()) / IK

    # Float search for a candidate u_hi (slightly inflated)
    b = lemma_l0009(n=n, k_ic=k_ic, ic_kind=ic_kind, E0=E0, t=t, cs0002_M=None)
    if not b.two_scale_closes or b.u_star is None:
        raise RuntimeError("two-scale branch does not close; cannot certify")
    u_hi = float(b.u_star) * (1.0 + 1e-6) + 1e-9

    # Inflate until interval verification succeeds
    verified = False
    F_gap_lo = float("-inf")
    Kt_hi = float((IK * It).hi)
    for _ in range(40):
        F_u = _F_interval(u_hi, U_L, lam)
        F_0 = _F_interval(0.0, U_L, lam)
        gap = F_u - F_0
        F_gap_lo = float(gap.lo)
        if F_gap_lo >= Kt_hi:
            verified = True
            break
        u_hi = u_hi * 1.000001 + 1e-8

    if not verified:
        raise RuntimeError("failed to certify F(u_hi)-F(0) >= K t")

    Iu = Interval.from_float(u_hi)
    cosh_u = Iu.cosh()
    cap = IB * cosh_u.sqr()

    payload = json.dumps(
        {"n": n, "ic_kind": ic_kind, "k_ic": k_ic, "K_ic2": K_ic2,
         "M_L": M_L, "M_H": M_H, "E0": E0, "t": t},
        sort_keys=True,
    ).encode()
    digest = hashlib.sha256(payload).hexdigest()

    return L0009Certificate(
        certificate_id=f"CERT-L0009-twoscale-{ic_kind}-k{k_ic}-N{n}",
        lemma_ref="L-0009 two-scale embedding (shell IC, full dealias)",
        route="B",
        evidence_level="N5",
        n=n,
        ic_kind=ic_kind,
        k_ic=k_ic,
        K_ic_squared=K_ic2,
        K_full_squared=K2_full,
        M_L=M_L,
        M_H=M_H,
        E0=E0,
        t=t,
        B_hi=float(IB.hi),
        U_L_hi=float(U_L.hi),
        W_hi=float(W.hi),
        u_hi=u_hi,
        F_gap_lo=F_gap_lo,
        Kt_hi=Kt_hi,
        Omega_bound_hi=float(cap.hi),
        closes=True,
        clay_implication="None. Finite Galerkin two-scale bound; not continuum NS.",
        payload_sha256=digest,
        notes=f"Verified F(u_hi)-F(0) >= Kt with u_hi={u_hi:.8f}; Ω_hi={float(cap.hi):.6e}.",
    )


def verify_l0009_certificate(cert: dict) -> tuple[bool, dict]:
    rebuilt = build_l0009_certificate(
        n=int(cert["n"]),
        k_ic=int(cert["k_ic"]),
        ic_kind=str(cert["ic_kind"]),
        E0=float(cert["E0"]),
        t=float(cert["t"]),
    )
    checks = {
        "K_ic_matches": rebuilt.K_ic_squared == int(cert["K_ic_squared"]),
        "M_L_matches": rebuilt.M_L == int(cert["M_L"]),
        "M_H_matches": rebuilt.M_H == int(cert["M_H"]),
        "bound_is_valid_upper": float(cert["Omega_bound_hi"]) >= rebuilt.Omega_bound_hi - 1e-3,
        "F_gap_certified": float(cert["F_gap_lo"]) >= float(cert["Kt_hi"]) - 1e-12,
        "evidence_level_N5": cert.get("evidence_level") == "N5",
        "no_clay_claim": "None" in cert.get("clay_implication", ""),
    }
    return all(checks.values()), checks


def save_l0009_certificate(cert: L0009Certificate, path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(cert.as_dict(), indent=2), encoding="utf-8")
