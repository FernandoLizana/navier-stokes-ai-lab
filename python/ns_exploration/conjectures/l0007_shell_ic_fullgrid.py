"""
L-0007: Short-time a priori bound for *shell-supported ICs* evolving on the
*full* dealiased Galerkin grid (the C-S conjecture class).

Gap closed relative to L-0005/L-0006:
  L-0005/L-0006 assume the dynamics stay inside the shell (hard truncation).
  C-S conjectures allow the nonlinear term to excite |k| > k_shell.
  L-0007 uses the honest worst-case stretch constants of the FULL dealias
  support, but the *initial* enstrophy ceiling of the SHELL:
      Ω(0) ≤ K_IC² E0,
  where K_IC is the max Euclidean |k| among modes admitted in the IC class
  (ℓ^∞ shell |k|_∞ ≤ k_inf, or Euclidean |k| ≤ k_eucl).

Differential inequality (same as L-0006, with M = M_full):
  dΩ/dt ≤ 2 a Ω^{3/2},   a = √2 √(3 M_full),
hence while a √Ω0 t < 1,
  Ω(t) ≤ Ω0 / (1 - a √Ω0 t)².
Fallback: L-0001-style exponential with (K_full, M_full) but Ω0 = K_IC² E0,
and L-0003 uniform ceiling on the full grid (independent of the IC shell).

This can be strictly sharper than "full-grid L-0006 from t=0" because Ω0 uses
K_IC ≤ K_full. It does NOT in general reach the sharp conjectured M of C-S-0002.

FINITE Galerkin ONLY. Not continuum NS. Not Clay. Evidence: N7.
"""

from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from ns_exploration.conjectures.l0003_uniform_galerkin import explicit_uniform_bound
from ns_exploration.conjectures.l0006_improved_shell import algebraic_short_time_bound
from ns_exploration.spectral.dealias import dealias_mask
from ns_exploration.spectral.fourier_conventions import k_squared
from ns_exploration.validation.l0003_certificate import full_dealias_exact_stats


def linf_shell_K2(k_inf: int) -> int:
    """Max Euclidean |k|² for integer modes with |k|_∞ ≤ k_inf."""
    # Corner mode (k_inf, k_inf, k_inf)
    return 3 * int(k_inf) ** 2


def euclidean_shell_K2(n: int, k_eucl: int) -> int:
    """Max |k|² on the dealiased grid with |k| ≤ k_eucl."""
    mask = dealias_mask(n)
    k2 = k_squared(n)
    kr = np.sqrt(k2)
    shell = mask & (kr <= float(k_eucl) + 1e-12)
    if not np.any(shell):
        return 0
    return int(np.rint(k2[shell]).max())


@dataclass
class GalerkinBoundL0007:
    lemma_id: str = "L-0007"
    route: str = "B"
    status: str = "proved_restricted"
    evidence_level: str = "N7"
    n: int = 16
    ic_kind: str = "linf"  # "linf" | "euclidean"
    k_ic: int = 4
    K_ic_squared: int = 0
    K_full: float = 0.0
    K_full_squared: int = 0
    M_full: int = 0
    E0: float = 0.5
    nu: float = 0.1
    t: float = 0.02
    Omega0_cap: float = 0.0
    a_stretch: float = 0.0
    algebraic_closes: bool = False
    algebraic_cap: float | None = None
    t_star_algebraic: float | None = None
    L0001_cap: float = 0.0
    L0003_cap: float = 0.0
    best_cap: float = 0.0
    which_best: str = ""
    proves_cs0002: bool = False
    cs0002_M: float | None = None
    clay_implication: str = (
        "None. Full-grid finite Galerkin with shell IC; constants grow with M_full."
    )
    notes: str = ""

    def as_dict(self) -> dict:
        return asdict(self)


def lemma_l0007(
    n: int = 16,
    k_ic: int = 4,
    ic_kind: str = "linf",
    E0: float = 0.5,
    t: float = 0.02,
    nu: float = 0.1,
    cs0002_M: float | None = 48.15928260103266,
) -> GalerkinBoundL0007:
    if ic_kind == "linf":
        K_ic2 = linf_shell_K2(k_ic)
    elif ic_kind == "euclidean":
        K_ic2 = euclidean_shell_K2(n, k_ic)
    else:
        raise ValueError(f"unknown ic_kind: {ic_kind}")

    K2_full, M_full = full_dealias_exact_stats(n)
    K_full = math.sqrt(float(K2_full))
    K_ic = math.sqrt(float(K_ic2))
    Omega0 = float(K_ic2) * E0

    closes, alg, a, t_star = algebraic_short_time_bound(E0, t, K_ic, M_full)
    # L-0001 exponential uses full-grid stretch scale but shell Ω0:
    # explicit_enstrophy_bound uses Ω0=K²E0 with the K passed in; pass K_ic for
    # Ω0 and ... wait, the alpha also uses K. For honest full-grid stretch,
    # alpha should use K_full. So compute manually:
    # α = 2 K_full √(3 M_full) √(2 E0), Ω(t) ≤ Ω0 e^{α t}
    alpha = 2.0 * K_full * math.sqrt(3.0 * M_full) * math.sqrt(2.0 * E0)
    L0001_cap = Omega0 * math.exp(alpha * t)
    b3 = explicit_uniform_bound(E0, K_full, M_full, nu=nu)

    candidates: list[tuple[str, float]] = [
        ("L-0007-exp", L0001_cap),
        ("L-0003-full", b3.Omega_uniform_cap),
    ]
    if closes and alg is not None:
        candidates.append(("L-0007-algebraic", alg))

    which, best = min(candidates, key=lambda x: x[1])
    proves = bool(cs0002_M is not None and best <= float(cs0002_M))

    notes = (
        f"IC {ic_kind} k={k_ic} on N={n}: K_IC²={K_ic2}, Ω0≤{Omega0:.6g}; "
        f"full dealias K²={K2_full}, M={M_full}, a={a:.6g}, t*={t_star:.6g}. "
        f"Algebraic closes={closes}. Best={which} → {best:.6e}. "
        f"Proves C-S-0002 M={cs0002_M}: {proves}."
    )
    return GalerkinBoundL0007(
        n=n,
        ic_kind=ic_kind,
        k_ic=k_ic,
        K_ic_squared=K_ic2,
        K_full=K_full,
        K_full_squared=K2_full,
        M_full=M_full,
        E0=E0,
        nu=nu,
        t=t,
        Omega0_cap=Omega0,
        a_stretch=a,
        algebraic_closes=closes,
        algebraic_cap=alg,
        t_star_algebraic=t_star,
        L0001_cap=L0001_cap,
        L0003_cap=b3.Omega_uniform_cap,
        best_cap=best,
        which_best=which,
        proves_cs0002=proves,
        cs0002_M=cs0002_M,
        notes=notes,
    )


def save_lemma_l0007(bound: GalerkinBoundL0007, path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(bound.as_dict(), indent=2), encoding="utf-8")
