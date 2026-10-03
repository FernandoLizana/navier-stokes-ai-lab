"""
L-0034: Shell-factorized convective bound on ‖N‖ (signed-geometry light).

Split u=Σ_r u_r by radial shells. The crude estimate
  ‖(u_s·∇)u_t‖₂ ≤ ‖u_s‖_∞ ‖∇u_t‖₂
with ‖u_s‖_∞ ≤ √(m_s) √(2 e_s) and ‖∇u_t‖₂ = √(2 r_t e_t) yields
  ‖N‖ ≤ Σ_{s,t} √(m_s)√(2 e_s)·√(2 r_t e_t)
       = 2 (Σ_s √(m_s e_s)) (Σ_t √(r_t e_t)).

Each factor Σ √(· e) is concave in the shell energies, hence maximized at
extreme (two-point / singleton) spectra under the linear constraints
Σ e_r=E, Σ r e_r=Ω. Therefore
  ‖N‖ ≤ 2 A_★(E,Ω) B_★(E,Ω),
where A_★,B_★ are those two-point maxima — a rigorous state-dependent bound.

Consequences on N=24 dealias:
- At (E,Ω)=(0.5,0.5): ‖N‖≤≈2.45 ≪ Young ≈78.7 (large gain at low enstrophy).
- Effective stretch C_eff=√(2K²Ω)‖N‖/Ω^{3/2} still ≳65–129 ≫ C_†≈9.56.
- Uniform CS envelope C≤2√(2K²)√(Σ m_r/r)√(#shells)≈3298 (useless for C_†).
- Low-slab comparison ODE with this N floors ≈87.5 > M≈41.28.

So shell factorization helps ‖N‖ at low Ω but does not reach C_† / C-0007.
Still missing: cancellation / signed triad structure inside each shell block.

FINITE dealias Galerkin only. Not continuum. Not Clay.
Evidence: N7 (shell estimates + two-point extrema) + N2 (ODE probe).
"""

from __future__ import annotations

import json
import math
from collections import Counter
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from ns_exploration.conjectures.l0020_duhamel_h1 import full_mask_rho_star
from ns_exploration.conjectures.l0027_low_slab_cubic import lemma_l0027
from ns_exploration.spectral.dealias import dealias_mask
from ns_exploration.spectral.fourier_conventions import k_squared
from ns_exploration.validation.l0003_certificate import full_dealias_exact_stats


def shell_multiplicities(n: int = 24) -> tuple[list[int], list[int]]:
    mask = dealias_mask(n)
    k2 = k_squared(n)
    ctr = Counter(np.rint(k2[mask & (k2 > 0)]).astype(int).tolist())
    radii = sorted(ctr.keys())
    return radii, [ctr[r] for r in radii]


def max_two_point_factor(
    E: float,
    Om: float,
    radii: list[int],
    masses: list[int],
    kind: str,
) -> float:
    """Max Σ √(m e) or Σ √(r e) over two-point spectra with given (E,Ω)."""
    best = 0.0
    E = float(E)
    Om = float(Om)
    for i, r1 in enumerate(radii):
        m1 = masses[i]
        for j in range(i, len(radii)):
            r2 = radii[j]
            m2 = masses[j]
            if r1 == r2:
                if abs(Om - r1 * E) <= 1e-9 * max(E, 1.0):
                    if kind == "m":
                        best = max(best, math.sqrt(m1 * E))
                    else:
                        best = max(best, math.sqrt(r1 * E))
                continue
            if Om < r1 * E - 1e-12 or Om > r2 * E + 1e-12:
                continue
            e2 = (Om - r1 * E) / (r2 - r1)
            e1 = E - e2
            if e1 < 0.0 or e2 < 0.0:
                continue
            if kind == "m":
                best = max(best, math.sqrt(m1 * e1) + math.sqrt(m2 * e2))
            else:
                best = max(best, math.sqrt(r1 * e1) + math.sqrt(r2 * e2))
    return float(best)


def N_shell_bound(E: float, Om: float, radii: list[int], masses: list[int], rho_star: float) -> float:
    a = max_two_point_factor(E, Om, radii, masses, "m")
    b = max_two_point_factor(E, Om, radii, masses, "r")
    return min(2.0 * a * b, 2.0 * float(E) * math.sqrt(float(rho_star)))


def uniform_C_shell_cs(radii: list[int], masses: list[int], K2: int) -> float:
    smr = sum(m / r for m, r in zip(masses, radii))
    return 2.0 * math.sqrt(2.0 * float(K2)) * math.sqrt(smr) * math.sqrt(len(radii))


def low_slab_shell_ode_worst(
    E0: float,
    nu: float,
    T: float,
    K2: int,
    Omega_hi: float,
    radii: list[int],
    masses: list[int],
    rho_star: float,
    n_grid: int = 11,
) -> float:
    def ode(Om0: float) -> float:
        def rhs(_t: float, y: np.ndarray) -> list[float]:
            E = max(float(y[0]), 1e-15)
            Om = min(max(float(y[1]), 1e-15), float(K2) * E)
            Nbd = N_shell_bound(E, Om, radii, masses, rho_star)
            return [-2.0 * nu * Om, -(2.0 * nu / E) * Om * Om + math.sqrt(2.0 * K2 * Om) * Nbd]

        sol = solve_ivp(
            rhs, (0.0, float(T)), [float(E0), float(Om0)], rtol=1e-7, atol=1e-10, max_step=5e-4
        )
        return float(sol.y[1, -1])

    xs = np.linspace(E0, Omega_hi, n_grid)
    return max(ode(float(o)) for o in xs)


@dataclass
class GalerkinBoundL0034:
    lemma_id: str = "L-0034"
    route: str = "B"
    status: str = "proved_restricted"
    evidence_level: str = "N7"
    n: int = 24
    c0007_M: float = 0.0
    C_dagger: float = 0.0
    n_shells: int = 0
    sum_m_over_r: float = 0.0
    N_at_equipartition: float = 0.0
    N_young: float = 0.0
    C_eff_at_equipartition: float = 0.0
    C_eff_at_Omega_star: float = 0.0
    uniform_C_cs: float = 0.0
    low_slab_ode_worst: float = 0.0
    closes_c0007: bool = False
    clay_implication: str = (
        "None. Shell-factorized ‖N‖ bound without triad cancellation; "
        "C_† not reached; not continuum; not Clay."
    )
    notes: str = ""

    def as_dict(self) -> dict:
        return asdict(self)


def lemma_l0034(
    n: int = 24,
    E0: float = 0.5,
    nu: float = 0.1,
    T: float = 0.02,
    n_grid: int = 11,
) -> GalerkinBoundL0034:
    b7 = lemma_l0027(n=n, empirical=False)
    radii, masses = shell_multiplicities(n)
    rho, _, _ = full_mask_rho_star(n)
    K2, _ = full_dealias_exact_stats(n)
    smr = sum(m / r for m, r in zip(masses, radii))
    N_eq = N_shell_bound(E0, E0, radii, masses, rho)
    N_y = 2.0 * E0 * math.sqrt(rho)
    N_star = N_shell_bound(E0, b7.Omega_star, radii, masses, rho)
    C_eq = math.sqrt(2.0 * K2 * E0) * N_eq / (E0**1.5)
    C_st = math.sqrt(2.0 * K2 * b7.Omega_star) * N_star / (b7.Omega_star**1.5)
    C_uni = uniform_C_shell_cs(radii, masses, int(K2))
    low = low_slab_shell_ode_worst(
        E0, nu, T, int(K2), b7.Omega_star, radii, masses, rho, n_grid=n_grid
    )
    notes = (
        f"Shell N≤2 A_★ B_★: at Ω=E={E0}, N≤{N_eq:.4g} vs Young {N_y:.4g}. "
        f"C_eff(E)={C_eq:.4g}, C_eff(Ω★)={C_st:.4g}, uniform CS C≤{C_uni:.4g} "
        f"(C_†={b7.C_dagger:.4g}). Low-slab ODE worst={low:.4g} (M={b7.c0007_M:.4g}). "
        f"Closes C-0007: False."
    )
    return GalerkinBoundL0034(
        n=n,
        c0007_M=b7.c0007_M,
        C_dagger=b7.C_dagger,
        n_shells=len(radii),
        sum_m_over_r=float(smr),
        N_at_equipartition=N_eq,
        N_young=N_y,
        C_eff_at_equipartition=C_eq,
        C_eff_at_Omega_star=C_st,
        uniform_C_cs=C_uni,
        low_slab_ode_worst=low,
        closes_c0007=False,
        notes=notes,
    )


def save_lemma_l0034(
    bound: GalerkinBoundL0034,
    path: str | Path = "conjectures/proved_restricted/L-0034.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(bound.as_dict(), indent=2), encoding="utf-8")
