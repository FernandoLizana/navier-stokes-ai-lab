"""
L-0031…L-0033: Exhaustive techos on the C-0007 low-slab / C_† route.

L-0031 — Shell-moment compatible (E,Ω) majorant.
  Replace L-0024's split (dissipation S≥Ω²/E, stretch √(2K²Ω)‖N‖) by a
  feasible-S envelope: for each state take max_{S∈[Smin,Smax]} of
  -2νS + √(2S) N_defect(E,Ω). Integrates to worst Ω(T)≈41.825 > M≈41.284
  and *worse* than L-0024's 41.743. Techo for autonomous (E,Ω)-only
  compatible-S refinements of the defect ODE.

L-0032 — Full-mask L-0015 S_max.
  S_max := max_p Σ_{q: p+q∈S} |q|² = 169574 on nonzero dealias S ⇒
  ‖N‖≤2E√S_max≈411.8 ≫ Young 2E√ρ_★≈78.7. Naive full-mask S_max is
  strictly worse than radial Young; the shell→high split was essential.

L-0033 — Absolute weighted cubic Frobenius for stretch.
  With a_k=|k|‖û_k‖/√(2Ω) and |stretch|≤(2Ω)^{3/2} Σ (|s|/|p|) a_p a_q a_s
  over p+q=s∈S, Cauchy gives C≤2√2 ‖c‖_F with ‖c‖_F≈3179 ⇒ C≲8992≪useful.
  Absolute triad majorants without cancellation geometry cannot reach C_†≈9.56.

Together with L-0026…L-0030: embeddings, raw R_★, compatible-S, full S_max,
and absolute cubic Frobenius are all blocked. Remaining path: cancellation-aware
weighted triad / SOS flattening of the *signed* stretch form.

Also records mono-radial shells with α_r=0 (no self-triads): Duhamel collapses
to the Stokes floor ≤M, proving C-0007 on that IC subclass (C-R-0003).

FINITE dealias Galerkin only. Not continuum. Not Clay.
Evidence: N7 (+ N2 ODE probes for L-0031).
"""

from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from ns_exploration.conjectures.l0020_duhamel_h1 import (
    full_mask_rho_star,
    lemma_l0020,
    omega_duhamel_H1,
)
from ns_exploration.conjectures.l0024_spectral_defect import (
    N_defect_bound,
    dealias_shell_stats,
)
from ns_exploration.conjectures.l0027_low_slab_cubic import lemma_l0027
from ns_exploration.conjectures.l0030_triad_N_bound import dealias_nonzero_modes
from ns_exploration.spectral.dealias import dealias_mask
from ns_exploration.spectral.fourier_conventions import wave_number_grids
from ns_exploration.validation.l0003_certificate import full_dealias_exact_stats


def dealias_radii(n: int = 24) -> list[int]:
    from collections import Counter
    from ns_exploration.spectral.fourier_conventions import k_squared

    mask = dealias_mask(n)
    k2 = k_squared(n)
    ctr = Counter(np.rint(k2[mask & (k2 > 0)]).astype(int).tolist())
    return sorted(ctr.keys())


def S_bounds(E: float, Om: float, radii: list[int]) -> tuple[float, float]:
    """Feasible range for S=Σ r² e_r given E,Ω (Cauchy + two-point extremes)."""
    E = max(float(E), 1e-15)
    lam = float(Om) / E
    r_min, r_max = radii[0], radii[-1]
    lam = min(max(lam, float(r_min)), float(r_max))
    # S ≥ Ω²/E
    S_lo = (Om * Om) / E
    # S ≤ two-point on (r_min, r_max)
    E_hi = (Om - r_min * E) / (r_max - r_min)
    E_lo = E - E_hi
    S_hi = r_min * r_min * max(E_lo, 0.0) + r_max * r_max * max(E_hi, 0.0)
    return float(S_lo), float(max(S_hi, S_lo))


def compatible_rhs(E: float, Om: float, nu: float, K2: int, stats: dict, radii: list[int]) -> float:
    """max_S (-2ν S + √(2S) N_defect) over feasible S."""
    Nbd = N_defect_bound(E, Om, stats)
    S_lo, S_hi = S_bounds(E, Om, radii)
    # f(S)= -2νS + √(2S) N is concave; max at endpoints or critical point
    def f(S: float) -> float:
        S = max(S, 1e-30)
        return -2.0 * nu * S + math.sqrt(2.0 * S) * Nbd

    crit = (Nbd * Nbd) / (8.0 * nu * nu) if nu > 0 else S_hi
    cands = [S_lo, S_hi]
    if S_lo <= crit <= S_hi:
        cands.append(crit)
    return max(f(S) for S in cands)


def lemma_l0031_compatible_worst(
    n: int = 24,
    E0: float = 0.5,
    nu: float = 0.1,
    T: float = 0.02,
    n_grid: int = 41,
) -> tuple[float, float]:
    stats = dealias_shell_stats(n)
    K2 = stats["K2"]
    radii = dealias_radii(n)

    def ode(Om0: float) -> float:
        def rhs(_t: float, y: np.ndarray) -> list[float]:
            E = max(float(y[0]), 1e-15)
            Om = min(max(float(y[1]), 1e-15), float(K2) * E)
            return [-2.0 * nu * Om, compatible_rhs(E, Om, nu, K2, stats, radii)]

        sol = solve_ivp(
            rhs, (0.0, T), [E0, Om0], rtol=1e-8, atol=1e-11, max_step=1e-4
        )
        return float(sol.y[1, -1])

    xs = np.linspace(E0, K2 * E0, n_grid)
    vals = [ode(float(o)) for o in xs]
    i = int(np.argmax(vals))
    return float(vals[i]), float(xs[i])


def full_mask_Smax(n: int = 24) -> float:
    modes = dealias_nonzero_modes(n)
    mset = set(modes)
    s_max = 0.0
    for p in modes:
        s = 0.0
        px, py, pz = p
        for qx, qy, qz in modes:
            if (px + qx, py + qy, pz + qz) in mset:
                s += float(qx * qx + qy * qy + qz * qz)
        if s > s_max:
            s_max = s
    return float(s_max)


def absolute_cubic_frobenius(n: int = 24) -> tuple[float, int]:
    """‖c‖_F for c_{pqs}=|s|/|p| on dealias triads; returns (F, n_tri)."""
    modes = dealias_nonzero_modes(n)
    mset = set(modes)
    abs_k = {
        p: math.sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]) for p in modes
    }
    F2 = 0.0
    n_tri = 0
    for p in modes:
        pk = abs_k[p]
        if pk < 1e-15:
            continue
        px, py, pz = p
        for q in modes:
            s = (px + q[0], py + q[1], pz + q[2])
            if s not in mset:
                continue
            c = abs_k[s] / pk
            F2 += c * c
            n_tri += 1
    return math.sqrt(F2), n_tri


def zero_alpha_shells_from_l0021() -> list[int]:
    path = Path("conjectures/proved_restricted/L-0021.json")
    if not path.exists():
        return []
    d = json.loads(path.read_text(encoding="utf-8"))
    alphas = d.get("alphas") or {}
    return sorted(int(r) for r, a in alphas.items() if abs(float(a)) < 1e-12)


@dataclass
class GalerkinBoundL0031_33:
    lemma_id: str = "L-0031+L-0032+L-0033"
    route: str = "B"
    status: str = "proved_restricted"
    evidence_level: str = "N7"
    n: int = 24
    c0007_M: float = 0.0
    C_dagger: float = 0.0
    # L-0031
    compatible_worst_OmT: float = 0.0
    compatible_worst_Om0: float = 0.0
    L0024_majorant: float = 0.0
    compatible_beats_L0024: bool = False
    # L-0032
    S_max_full: float = 0.0
    N_from_Smax: float = 0.0
    N_young: float = 0.0
    Smax_beats_young: bool = False
    # L-0033
    cubic_F: float = 0.0
    n_triads: int = 0
    C_from_cubic_F: float = 0.0
    # C-R-0003 seeds
    zero_alpha_shells: list[int] | None = None
    zero_alpha_duhamel_OmT: float = 0.0
    closes_c0007: bool = False
    clay_implication: str = (
        "None. Techos on compatible-S / full S_max / absolute cubic; "
        "plus α=0 mono-radial subclass; not continuum; not Clay."
    )
    notes: str = ""

    def as_dict(self) -> dict:
        return asdict(self)


def lemma_l0031_33(
    n: int = 24,
    E0: float = 0.5,
    nu: float = 0.1,
    T: float = 0.02,
    n_grid: int = 41,
) -> GalerkinBoundL0031_33:
    b7 = lemma_l0027(n=n, empirical=False)
    from ns_exploration.conjectures.l0026_c0007_techo import lemma_l0026

    b26 = lemma_l0026(n=n, n_grid=41)
    omT, om0 = lemma_l0031_compatible_worst(n, E0, nu, T, n_grid=n_grid)
    Smax = full_mask_Smax(n)
    rho, _, _ = full_mask_rho_star(n)
    N_s = 2.0 * E0 * math.sqrt(Smax)
    N_y = 2.0 * E0 * math.sqrt(rho)
    F, n_tri = absolute_cubic_frobenius(n)
    C_F = 2.0 * math.sqrt(2.0) * F
    zshells = zero_alpha_shells_from_l0021()
    b20 = lemma_l0020(n=n)
    z_om = (
        omega_duhamel_H1(b20.S_T, E0, 0.0, b20.I_sigma) if zshells else float("nan")
    )
    notes = (
        f"L-0031: compatible-S worst Ω(T)={omT:.6f} at Ω0={om0:.4f} "
        f"(L-0024={b26.L0024_majorant:.6f}, M={b7.c0007_M:.6f}). "
        f"L-0032: S_max={Smax:.0f}, N≤{N_s:.4f} vs Young {N_y:.4f}. "
        f"L-0033: ‖c‖_F={F:.4f}, C≤{C_F:.4f}≫C_†={b7.C_dagger:.4f}. "
        f"α=0 shells {zshells}: Duhamel={z_om:.6f}≤M. C-0007 all-IC still open."
    )
    return GalerkinBoundL0031_33(
        n=n,
        c0007_M=b7.c0007_M,
        C_dagger=b7.C_dagger,
        compatible_worst_OmT=omT,
        compatible_worst_Om0=om0,
        L0024_majorant=b26.L0024_majorant,
        compatible_beats_L0024=omT < b26.L0024_majorant - 1e-9,
        S_max_full=Smax,
        N_from_Smax=N_s,
        N_young=N_y,
        Smax_beats_young=N_s < N_y - 1e-9,
        cubic_F=F,
        n_triads=n_tri,
        C_from_cubic_F=C_F,
        zero_alpha_shells=zshells,
        zero_alpha_duhamel_OmT=float(z_om),
        closes_c0007=False,
        notes=notes,
    )


def save_lemma_l0031_33(
    bound: GalerkinBoundL0031_33,
    path: str | Path = "conjectures/proved_restricted/L-0031_0032_0033.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(bound.as_dict(), indent=2), encoding="utf-8")
