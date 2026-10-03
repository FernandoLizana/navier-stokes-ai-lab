"""
L-0031: Shell-weighted triad multiplicity (toward C_† all-IC).

Refines L-0030 by weighting each triad p+q=s with w(p,q)=1/sqrt(r_s)
where r_s = |s|^2. This down-weights high-shell triads in the counting
constant R_w = max_s sum_{p+q=s} w(p,q).

Still finite Galerkin; not continuum; not Clay.
"""

from __future__ import annotations

import json
import math
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from ns_exploration.conjectures.l0020_duhamel_h1 import full_mask_rho_star
from ns_exploration.conjectures.l0024_spectral_defect import N_defect_bound, dealias_shell_stats
from ns_exploration.conjectures.l0027_low_slab_cubic import lemma_l0027
from ns_exploration.conjectures.l0030_triad_N_bound import (
    N_hybrid_bound,
    N_triad_bound,
    N_young_bound,
    cubic_C_from_Rstar,
    dealias_nonzero_modes,
    full_mask_triad_Rstar,
    low_slab_ode_worst,
)
from ns_exploration.validation.l0003_certificate import full_dealias_exact_stats


def _radius(k: tuple[int, int, int]) -> int:
    return int(k[0] * k[0] + k[1] * k[1] + k[2] * k[2])


def weighted_triad_Rstar(
    n: int = 24,
    *,
    weight: str = "inv_sqrt_rs",
) -> tuple[float, int, int]:
    """
    Return (R_w, R_unweighted, n_triads).
    R_w = max_s sum_{p+q=s} w(p,q) with w = 1/sqrt(r_s) by default.
    """
    modes = dealias_nonzero_modes(n)
    mset = set(modes)
    ctr: Counter[tuple[int, int, int]] = Counter()
    wsum: defaultdict[tuple[int, int, int], float] = defaultdict(float)
    for p in modes:
        px, py, pz = p
        for qx, qy, qz in modes:
            q = (qx, qy, qz)
            s = (px + qx, py + qy, pz + qz)
            if s not in mset:
                continue
            ctr[s] += 1
            rs = max(_radius(s), 1)
            rp = max(_radius(p), 1)
            rq = max(_radius(q), 1)
            if weight == "inv_sqrt_rs":
                w = 1.0 / math.sqrt(float(rs))
            elif weight == "inv_rs":
                w = 1.0 / float(rs)
            elif weight == "inv_sqrt_rp_rq":
                w = 1.0 / math.sqrt(float(rp) * float(rq))
            else:
                w = 1.0
            wsum[s] += w
    if not ctr:
        return 0.0, 0, 0
    R_u = int(max(ctr.values()))
    R_w = float(max(wsum.values()))
    return R_w, R_u, int(sum(ctr.values()))


def N_weighted_triad_bound(E: float, Omega: float, R_w: float) -> float:
    return 2.0 * math.sqrt(float(R_w) * max(float(E), 0.0) * max(float(Omega), 0.0))


def cubic_C_from_weighted_R(R_w: float) -> float:
    return 2.0 * math.sqrt(2.0) * math.sqrt(float(R_w))


@dataclass
class GalerkinBoundL0031:
    lemma_id: str = "L-0031"
    route: str = "B"
    status: str = "exploring"
    evidence_level: str = "N7"
    n: int = 24
    R_star_unweighted: int = 0
    R_star_weighted: float = 0.0
    weight_scheme: str = "inv_sqrt_rs"
    C_dagger: float = 0.0
    C_from_Rstar: float = 0.0
    C_from_weighted_R: float = 0.0
    low_slab_weighted_worst: float = 0.0
    low_slab_hybrid_weighted_worst: float = 0.0
    closes_c0007: bool = False
    clay_implication: str = "None. Weighted triad sketch; not all-IC proof; not Clay."
    notes: str = ""

    def as_dict(self) -> dict:
        return asdict(self)


def lemma_l0031(n: int = 24, *, weight: str = "inv_sqrt_rs") -> GalerkinBoundL0031:
    b7 = lemma_l0027(n=n, empirical=False)
    R_w, R_u, n_tri = weighted_triad_Rstar(n, weight=weight)
    rho, _, _ = full_mask_rho_star(n)
    K2, _ = full_dealias_exact_stats(n)
    stats = dealias_shell_stats(n)
    C_w = cubic_C_from_weighted_R(R_w)
    C_u = cubic_C_from_Rstar(R_u)

    def weighted_triad(E: float, Om: float) -> float:
        return N_weighted_triad_bound(E, Om, R_w)

    def hybrid_w(E: float, Om: float) -> float:
        return min(
            N_young_bound(E, rho),
            N_weighted_triad_bound(E, Om, R_w),
            N_triad_bound(E, Om, R_u),
        )

    def min_def_w(E: float, Om: float) -> float:
        return min(hybrid_w(E, Om), N_defect_bound(E, Om, stats))

    w_worst = low_slab_ode_worst(
        weighted_triad, b7.E0, b7.nu, b7.T, int(K2), b7.Omega_star, n_grid=15
    )
    hy_worst = low_slab_ode_worst(
        min_def_w, b7.E0, b7.nu, b7.T, int(K2), b7.Omega_star, n_grid=15
    )
    notes = (
        f"Weighted triad ({weight}): R_w={R_w:.4g} vs R_★={R_u}. "
        f"C≤2√2√R_w={C_w:.4g} (unweighted {C_u:.4g}, C_†={b7.C_dagger:.4g}). "
        f"Low-slab weighted Ω(T)≤{w_worst:.4g}, min(hybrid_w,defect)≤{hy_worst:.4g} "
        f"(M={b7.c0007_M:.4g}); closes: False."
    )
    return GalerkinBoundL0031(
        n=n,
        R_star_unweighted=R_u,
        R_star_weighted=R_w,
        weight_scheme=weight,
        C_dagger=b7.C_dagger,
        C_from_Rstar=C_u,
        C_from_weighted_R=C_w,
        low_slab_weighted_worst=w_worst,
        low_slab_hybrid_weighted_worst=hy_worst,
        closes_c0007=False,
        notes=notes,
    )


def save_lemma_l0031(
    bound: GalerkinBoundL0031,
    path: str | Path = "conjectures/active/L-0031.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(bound.as_dict(), indent=2), encoding="utf-8")
