"""Route D: perturbation from t=T with shell blocks."""

from __future__ import annotations

from gmpy2 import mpfr

from ns_exploration.terminal_weighted.constants import FROZEN
from ns_exploration.terminal_weighted.intervals import (
    add_up,
    div_down,
    div_up,
    exp_down,
    exp_up,
    i_star_interval,
    mul_up,
    sqrt_up,
    sub_down,
    sub_up,
    _ctx,
)
from ns_exploration.terminal_weighted.shell_decomposition import build_shell_manifest


def route_d_perturbation_certificate(
    n: int = 24,
    prec: int = 200,
    manifest: dict | None = None,
) -> dict:
    """
    int_0^T C(t) dt <= C(T)*T + sum_r B_r * (T - (1-e^{-2 nu r T})/(2 nu r))

    with C(T) bounded by full Frobenius at t=T and B_r per-shell Frobenius bounds.
    """
    _ctx(prec)
    manifest = manifest or build_shell_manifest(n, prec=prec)
    nu = mpfr(FROZEN.nu, precision=prec)
    T = mpfr(FROZEN.T, precision=prec)
    two = mpfr(2, precision=prec)

    C_T_hi = mpfr(
        manifest.get("full_best_hi", manifest["full_frobenius_hi"])["C_term_hi"],
        precision=prec,
    )
    integral_hi = mul_up(C_T_hi, T)

    for blk in manifest["blocks"]:
        r = mpfr(blk["shell"], precision=prec)
        B_r = mpfr(blk["C_term_hi"], precision=prec)
        expo_neg = mul_up(mul_up(mul_up(two, nu), r), T)
        e_up = exp_up(mul_up(mpfr(-1, precision=prec), expo_neg))
        one = mpfr(1, precision=prec)
        correction_lo = div_down(sub_down(one, e_up), mul_up(mul_up(two, nu), r))
        int_phi = sub_up(T, correction_lo)
        contrib = mul_up(B_r, int_phi)
        integral_hi = add_up(integral_hi, contrib)

    I_lo, I_hi = i_star_interval(prec=prec)
    two_sqrt2 = sqrt_up(mpfr(8, precision=prec))
    with _ctx(prec, up=True):
        floor_hi = mpfr(FROZEN.stokes_floor)
    omega_hi = add_up(floor_hi, div_up(integral_hi, two_sqrt2))

    return {
        "method": "route_D_terminal_perturbation",
        "integral_hi": str(integral_hi),
        "I_star_lo": str(I_lo),
        "I_star_hi": str(I_hi),
        "closes_strict": bool(integral_hi < I_lo),
        "omega_T_hi": str(omega_hi),
        "evidence_level": "N5",
        "honesty": "Adds C(T)*T base; typically looser than Route A.",
    }
