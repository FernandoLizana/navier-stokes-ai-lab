"""Spectral weights W(t) = A exp(-2 nu A (T-t)) for terminal functional."""

from __future__ import annotations

import math

from ns_exploration.conjectures.l0021_quartic_shell import shell_modes_by_r
from ns_exploration.terminal_weighted.constants import FROZEN


def terminal_weight(r: float, t: float, nu: float, T: float) -> float:
    """w_r(t) = r exp(-2 nu r (T - t)). At t=T returns r."""
    return float(r) * math.exp(-2.0 * float(nu) * float(r) * (float(T) - float(t)))


def terminal_weight_at_T_equals_r(r: float) -> float:
    """Regression check: w_r(T) = r."""
    return terminal_weight(r, FROZEN.T, FROZEN.nu, FROZEN.T)


def stokes_floor_hi(
    n: int = 24,
    E0: float = FROZEN.E0,
    nu: float = FROZEN.nu,
    T: float = FROZEN.T,
    *,
    margin: float = 1e-12,
) -> float:
    """
    Certified upper bound on Phi(0) for E <= E0:
      Phi(0) <= E0 * max_{r in dealias shells} r exp(-2 nu r T).
    """
    by = shell_modes_by_r(n)
    best = 0.0
    for r in by:
        w = terminal_weight(float(r), 0.0, nu, T)
        if w > best:
            best = w
    return float(E0) * best + margin


def phi_terminal(u_coeff_energy: float, r: float, t: float, nu: float, T: float) -> float:
    """Phi(t) contribution from a single mode with |u_k|^2 = 2 * mode_energy."""
    return 0.5 * terminal_weight(r, t, nu, T) * (2.0 * u_coeff_energy)
