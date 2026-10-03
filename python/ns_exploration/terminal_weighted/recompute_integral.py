"""Recompute certified integral bounds from shell block data (verifier support)."""

from __future__ import annotations

from gmpy2 import mpfr

from ns_exploration.terminal_weighted.constants import EXACT
from ns_exploration.terminal_weighted.intervals import (
    add_up,
    div_up,
    exp_down,
    mpfr_const,
    mul_up,
    omega_T_hi_from_I_term,
    _ctx,
)


def shell_integral_hi_from_blocks(
    blocks: list[dict],
    *,
    prec: int = 200,
) -> mpfr:
    """Route A: sum_r B_r * (1 - exp(-2 nu r T)) / (2 nu r), upward rounded."""
    _ctx(prec)
    nu = mpfr_const(EXACT.NU, prec=prec)
    T = mpfr_const(EXACT.T, prec=prec)
    two = mpfr(2, precision=prec)
    integral_hi = mpfr(0, precision=prec)
    for blk in blocks:
        r = mpfr(blk["shell"], precision=prec)
        B_r = mpfr(blk["C_term_hi"], precision=prec)
        expo = mul_up(mul_up(mpfr(-2, precision=prec), nu), mul_up(r, T))
        e_down = exp_down(expo)
        with _ctx(prec, up=True):
            factor_num_hi = mpfr(1) - e_down
        denom = mul_up(mul_up(two, nu), r)
        factor = div_up(factor_num_hi, denom)
        integral_hi = add_up(integral_hi, mul_up(B_r, factor))
    return integral_hi


def cluster_integral_hi_from_blocks(blocks: list[dict], *, prec: int = 200) -> mpfr:
    total = mpfr(0, precision=prec)
    for blk in blocks:
        total = add_up(total, mpfr(blk["contrib_hi"], precision=prec))
    return total


def omega_T_hi_from_integral(integral_hi: mpfr, *, prec: int = 200) -> mpfr:
    return omega_T_hi_from_I_term(integral_hi, prec=prec)
