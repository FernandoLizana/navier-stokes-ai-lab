"""Recompute certified integral bounds from shell block data (verifier support)."""

from __future__ import annotations

from gmpy2 import mpfr

from ns_exploration.terminal_weighted.constants import FROZEN
from ns_exploration.terminal_weighted.intervals import (
    add_up,
    div_up,
    exp_down,
    mul_up,
    sqrt_up,
    _ctx,
)


def shell_integral_hi_from_blocks(
    blocks: list[dict],
    *,
    prec: int = 200,
) -> mpfr:
    """Route A: sum_r B_r * (1 - exp(-2 nu r T)) / (2 nu r), upward rounded."""
    _ctx(prec)
    nu = mpfr(FROZEN.nu, precision=prec)
    T = mpfr(FROZEN.T, precision=prec)
    two = mpfr(2, precision=prec)
    integral_hi = mpfr(0, precision=prec)
    for blk in blocks:
        r = mpfr(blk["shell"], precision=prec)
        B_r = mpfr(blk["C_term_hi"], precision=prec)
        expo = mul_up(mul_up(mul_up(mpfr(-2, precision=prec), nu), r), T)
        e_down = exp_down(expo)
        with _ctx(prec, up=True):
            factor_num_hi = mpfr(1) - e_down
        denom = mul_up(mul_up(two, nu), r)
        factor = div_up(factor_num_hi, denom)
        integral_hi = add_up(integral_hi, mul_up(B_r, factor))
    return integral_hi


def omega_T_hi_from_integral(integral_hi: mpfr, *, prec: int = 200) -> mpfr:
    two_sqrt2 = sqrt_up(mpfr(8, precision=prec))
    with _ctx(prec, up=True):
        floor_hi = mpfr(FROZEN.stokes_floor)
    return add_up(floor_hi, div_up(integral_hi, two_sqrt2))
