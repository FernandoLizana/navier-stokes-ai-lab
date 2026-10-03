"""Certified ||M||_op bounds via sqrt(||M||_1 ||M||_inf) (MPFR upward)."""

from __future__ import annotations

from gmpy2 import mpfr

from ns_exploration.conjectures.l0040_sym_shor import _sym_col_index
from ns_exploration.terminal_weighted.intervals import add_up, mul_up, sqrt_up, _ctx


def abs_up(delta: float, prec: int) -> mpfr:
    with _ctx(prec, up=True):
        d = mpfr(delta)
        return mpfr(abs(d))


def one_inf_opnorm_hi(row_sum: dict[int, mpfr], col_sum: dict[int, mpfr]) -> mpfr:
    if not row_sum or not col_sum:
        return mpfr(0)
    norm_inf = mpfr(0, precision=next(iter(row_sum.values())).precision)
    norm_1 = mpfr(0, precision=norm_inf.precision)
    for v in row_sum.values():
        if v > norm_inf:
            norm_inf = v
    for v in col_sum.values():
        if v > norm_1:
            norm_1 = v
    return sqrt_up(mul_up(norm_1, norm_inf))


def one_inf_opnorm_hi_rowlist(row_sum: list[mpfr], col_sum: dict[int, mpfr]) -> mpfr:
    if not row_sum or not col_sum:
        return mpfr(0)
    prec = row_sum[0].precision if row_sum else 128
    norm_inf = mpfr(0, precision=prec)
    norm_1 = mpfr(0, precision=prec)
    for v in row_sum:
        if v > norm_inf:
            norm_inf = v
    for v in col_sum.values():
        if v > norm_1:
            norm_1 = v
    return sqrt_up(mul_up(norm_1, norm_inf))


def merge_sum_dicts(
    target: dict[int, mpfr], source: dict[int, mpfr], *, prec: int
) -> dict[int, mpfr]:
    for k, v in source.items():
        ki = int(k) if isinstance(k, str) else k
        if ki in target:
            target[ki] = add_up(target[ki], mpfr(v, precision=prec))
        else:
            target[ki] = mpfr(v, precision=prec)
    return target
