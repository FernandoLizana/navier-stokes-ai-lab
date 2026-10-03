"""
N5 certificate for L-0015: S_max energy bound + L-0012 cancellation ODE.
"""

from __future__ import annotations

import hashlib
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from ns_exploration.conjectures.l0007_shell_ic_fullgrid import linf_shell_K2
from ns_exploration.conjectures.l0010_duhamel_short import kappa_H_min
from ns_exploration.conjectures.l0012_cancellation import max_T_cancellation
from ns_exploration.conjectures.l0013_embedding import (
    embedding_W_vector,
    nonzero_high_mode_count,
    nonzero_shell_mode_count,
)
from ns_exploration.conjectures.l0015_smax import (
    Nmax_from_Smax,
    effective_U_from_Smax,
    shell_high_Smax,
)
from ns_exploration.validation.intervals import Interval
from ns_exploration.validation.l0003_certificate import full_dealias_exact_stats


@dataclass
class L0015ShortCertificate:
    certificate_id: str
    lemma_ref: str
    route: str
    evidence_level: str
    n: int
    k_ic: int
    K_ic_squared: int
    K_full_squared: int
    M_L: int
    M_H: int
    S_max: float
    N_max: float
    kappa_H: float
    E0: float
    nu: float
    R: float
    T_cert: float
    T_star: float
    U_eff_hi: float
    W_hi: float
    alpha_hi: float
    beta_hi: float
    z_hi: float
    Omega_bound_hi: float
    closes: bool
    clay_implication: str
    payload_sha256: str
    notes: str = ""

    def as_dict(self) -> dict:
        return asdict(self)


def _exp_outward(x: Interval) -> Interval:
    return Interval(
        float(np.nextafter(math.exp(x.lo), -np.inf)),
        float(np.nextafter(math.exp(x.hi), np.inf)),
    )


def build_l0015_short_certificate(
    n: int = 16,
    k_ic: int = 4,
    E0: float = 0.5,
    nu: float = 0.1,
    R: float = 48.15928260103266,
    safety: float = 0.99,
) -> L0015ShortCertificate:
    K_ic2 = linf_shell_K2(k_ic)
    K2_full, _ = full_dealias_exact_stats(n)
    M_L = nonzero_shell_mode_count(n, k_ic, "linf")
    M_H = nonzero_high_mode_count(n, k_ic, "linf")
    S_max = shell_high_Smax(n, k_ic, "linf")
    N_max = Nmax_from_Smax(S_max, E0)
    kappa = float(kappa_H_min(n, k_ic, "linf") or 0.0)
    K_full = math.sqrt(float(K2_full))
    U_eff = effective_U_from_Smax(S_max, E0, K_ic2)

    IE0 = Interval.from_float(E0)
    Inu = Interval.from_float(nu)
    Ikap = Interval.from_float(kappa)
    IB = Interval(float(K_ic2), float(K_ic2)) * IE0
    IK2 = Interval(float(K2_full), float(K2_full))
    # U_eff = N_max / (√2 √B) with N_max = 2 E0 √S_max, B = K_ic2 E0
    IS = Interval(float(S_max), float(S_max))
    INmax = Interval(2.0, 2.0) * IE0 * IS.sqrt()
    IUL = INmax / ((Interval(2.0, 2.0).sqrt()) * IB.sqrt())
    IW = (Interval(2.0, 2.0) * Interval(float(M_H), float(M_H))).sqrt()

    sqrtB = IB.sqrt()
    beta = IUL * sqrtB
    alpha = IW * sqrtB - Inu * Ikap

    B = float(K_ic2) * E0
    Wf = embedding_W_vector(M_H)
    T_star = max_T_cancellation(B, K_full, U_eff, Wf, R, nu, kappa)
    T_cert = safety * T_star
    IT = Interval.from_float(T_cert)

    if alpha.lo <= 0:
        raise RuntimeError("expected alpha > 0 for C-S N=16 regime")
    exp_at = _exp_outward(alpha * IT)
    z = (beta / alpha) * (exp_at - Interval(1.0, 1.0))
    Omega = IB + IK2 * z.sqr()

    closes = bool(Omega.hi <= R + 1e-9)
    if not closes:
        raise RuntimeError(f"Omega.hi={Omega.hi} exceeds R={R} at T_cert")

    payload = json.dumps(
        {
            "n": n,
            "k_ic": k_ic,
            "R": R,
            "T_cert": T_cert,
            "S_max": S_max,
            "lemma": "L-0015",
        },
        sort_keys=True,
    ).encode()
    digest = hashlib.sha256(payload).hexdigest()

    return L0015ShortCertificate(
        certificate_id="CERT-L0015-smax-CS0002-short-N16",
        lemma_ref="L-0015 S_max energy bound + L-0012 cancellation ODE",
        route="B",
        evidence_level="N5",
        n=n,
        k_ic=k_ic,
        K_ic_squared=K_ic2,
        K_full_squared=K2_full,
        M_L=M_L,
        M_H=M_H,
        S_max=S_max,
        N_max=N_max,
        kappa_H=kappa,
        E0=E0,
        nu=nu,
        R=R,
        T_cert=T_cert,
        T_star=T_star,
        U_eff_hi=float(IUL.hi),
        W_hi=float(IW.hi),
        alpha_hi=float(alpha.hi),
        beta_hi=float(beta.hi),
        z_hi=float(z.hi),
        Omega_bound_hi=float(Omega.hi),
        closes=closes,
        clay_implication="None. Short-time finite Galerkin only; not C-S-0002 at T=0.02.",
        payload_sha256=digest,
        notes=(
            f"S_max={S_max:.8f}; N_max={N_max:.8f}; T_cert={T_cert:.8f}; "
            f"Omega_hi={float(Omega.hi):.8f} <= R={R:.8f}."
        ),
    )


def verify_l0015_short_certificate(cert: dict) -> tuple[bool, dict]:
    rebuilt = build_l0015_short_certificate(
        n=int(cert["n"]),
        k_ic=int(cert["k_ic"]),
        E0=float(cert["E0"]),
        nu=float(cert["nu"]),
        R=float(cert["R"]),
        safety=float(cert["T_cert"]) / float(cert["T_star"])
        if float(cert["T_star"]) > 0
        else 0.99,
    )
    checks = {
        "K_ic_matches": rebuilt.K_ic_squared == int(cert["K_ic_squared"]),
        "M_L_matches": rebuilt.M_L == int(cert["M_L"]),
        "M_H_matches": rebuilt.M_H == int(cert["M_H"]),
        "S_max_matches": abs(rebuilt.S_max - float(cert["S_max"])) < 1e-9,
        "closes": rebuilt.closes and bool(cert.get("closes")),
        "Omega_valid_upper": float(cert["Omega_bound_hi"]) >= rebuilt.Omega_bound_hi - 1e-6,
        "Omega_le_R": float(cert["Omega_bound_hi"]) <= float(cert["R"]) + 1e-6,
        "short_horizon": float(cert["T_cert"]) < 0.02,
        "evidence_level_N5": cert.get("evidence_level") == "N5",
        "no_clay_claim": "None" in cert.get("clay_implication", ""),
        "sharper_than_L0014_T": float(cert["T_star"]) > 0.0043669005451596726,
    }
    return all(checks.values()), checks


def save_l0015_short_certificate(cert: L0015ShortCertificate, path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(cert.as_dict(), indent=2), encoding="utf-8")
