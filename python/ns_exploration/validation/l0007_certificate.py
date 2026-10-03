"""
Interval certificate for the L-0007 exponential short-time bound
(shell IC → full dealias grid):

  Ω(t) ≤ Ω0 · exp(α t),
  Ω0 = K_IC² E0,
  α = 2 K_full √(3 M_full) √(2 E0).

Uses outward rounding. Finite Galerkin only. Not Clay. Evidence: N5.
"""

from __future__ import annotations

import hashlib
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

from ns_exploration.conjectures.l0007_shell_ic_fullgrid import linf_shell_K2, euclidean_shell_K2
from ns_exploration.validation.intervals import Interval
from ns_exploration.validation.l0003_certificate import full_dealias_exact_stats


@dataclass
class L0007ExpCertificate:
    certificate_id: str
    lemma_ref: str
    route: str
    evidence_level: str
    support: str
    n: int
    ic_kind: str
    k_ic: int
    K_ic_squared: int
    K_full_squared: int
    M_full: int
    E0: float
    t: float
    Omega0_hi: float
    alpha_hi: float
    Omega_bound_hi: float
    closes: bool
    clay_implication: str
    payload_sha256: str
    notes: str = ""

    def as_dict(self) -> dict:
        return asdict(self)


def _exp_outward(x: Interval) -> Interval:
    """Outward enclosure of exp on a bounded interval (via math.exp + nextafter)."""
    import numpy as np

    lo = float(np.nextafter(math.exp(x.lo), -np.inf))
    hi = float(np.nextafter(math.exp(x.hi), np.inf))
    return Interval(lo, hi)


def build_l0007_exp_certificate(
    n: int = 12,
    k_ic: int = 2,
    ic_kind: str = "euclidean",
    E0: float = 0.5,
    t: float = 0.02,
) -> L0007ExpCertificate:
    if ic_kind == "linf":
        K_ic2 = linf_shell_K2(k_ic)
    elif ic_kind == "euclidean":
        K_ic2 = euclidean_shell_K2(n, k_ic)
    else:
        raise ValueError(ic_kind)

    K2_full, M_full = full_dealias_exact_stats(n)
    IE0 = Interval.from_float(E0)
    It = Interval.from_float(t)
    IK_ic2 = Interval(float(K_ic2), float(K_ic2))
    IK2f = Interval(float(K2_full), float(K2_full))
    IM = Interval(float(M_full), float(M_full))

    Omega0 = IK_ic2 * IE0
    # α = 2 * √K2_full * √(3M) * √(2 E0)
    two = Interval(2.0, 2.0)
    three = Interval(3.0, 3.0)
    K_full = IK2f.sqrt()
    a_part = (three * IM).sqrt()  # √(3M)
    e_part = (two * IE0).sqrt()  # √(2 E0)
    alpha = two * K_full * a_part * e_part
    prod = alpha * It
    exp_prod = _exp_outward(prod)
    cap = Omega0 * exp_prod

    payload = json.dumps(
        {"n": n, "ic_kind": ic_kind, "k_ic": k_ic, "K_ic2": K_ic2, "K2_full": K2_full,
         "M": M_full, "E0": E0, "t": t},
        sort_keys=True,
    ).encode()
    digest = hashlib.sha256(payload).hexdigest()

    return L0007ExpCertificate(
        certificate_id=f"CERT-L0007-exp-{ic_kind}-k{k_ic}-N{n}",
        lemma_ref="L-0007 exponential (shell IC, full dealias stretch)",
        route="B",
        evidence_level="N5",
        support="shell_ic_full_dealias",
        n=n,
        ic_kind=ic_kind,
        k_ic=k_ic,
        K_ic_squared=K_ic2,
        K_full_squared=K2_full,
        M_full=M_full,
        E0=E0,
        t=t,
        Omega0_hi=float(Omega0.hi),
        alpha_hi=float(alpha.hi),
        Omega_bound_hi=float(cap.hi),
        closes=True,
        clay_implication="None. Finite Galerkin shell-IC/full-grid; not continuum NS.",
        payload_sha256=digest,
        notes=(
            f"Outward: Ω(t)≤Ω0 exp(αt) ≤ {float(cap.hi):.6e} "
            f"(Ω0_hi={float(Omega0.hi):.6e}, α_hi={float(alpha.hi):.6e})."
        ),
    )


def verify_l0007_exp_certificate(cert: dict) -> tuple[bool, dict]:
    checks: dict[str, bool] = {}
    rebuilt = build_l0007_exp_certificate(
        n=int(cert["n"]),
        k_ic=int(cert["k_ic"]),
        ic_kind=str(cert["ic_kind"]),
        E0=float(cert["E0"]),
        t=float(cert["t"]),
    )
    checks["K_ic_matches"] = rebuilt.K_ic_squared == int(cert["K_ic_squared"])
    checks["K_full_matches"] = rebuilt.K_full_squared == int(cert["K_full_squared"])
    checks["M_matches"] = rebuilt.M_full == int(cert["M_full"])
    checks["bound_is_valid_upper"] = float(cert["Omega_bound_hi"]) >= rebuilt.Omega_bound_hi - 1e-6
    checks["evidence_level_N5"] = cert.get("evidence_level") == "N5"
    checks["no_clay_claim"] = "None" in cert.get("clay_implication", "")
    return all(checks.values()), checks


def save_l0007_exp_certificate(cert: L0007ExpCertificate, path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(cert.as_dict(), indent=2), encoding="utf-8")
