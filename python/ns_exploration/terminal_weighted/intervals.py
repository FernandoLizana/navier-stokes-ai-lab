"""
Directed rounding utilities for certified bounds (gmpy2 MPFR).

Critical constants enter as decimal strings (see constants.EXACT).
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any

import gmpy2
from gmpy2 import RoundDown, RoundUp, mpfr, sqrt

from ns_exploration.terminal_weighted.constants import EXACT

_DEFAULT_PREC = 200
_STOKES_META: dict[str, Any] = {}


@dataclass(frozen=True)
class ArithmeticBackend:
    name: str = "gmpy2_mpfr"
    version: str = ""
    precision_bits: int = _DEFAULT_PREC
    rounding_up: str = "RoundUp"
    rounding_down: str = "RoundDown"
    platform: str = "python"

    def as_dict(self) -> dict:
        import platform

        return {
            "name": self.name,
            "version": gmpy2.version(),
            "precision_bits": self.precision_bits,
            "rounding_up": self.rounding_up,
            "rounding_down": self.rounding_down,
            "platform": platform.platform(),
        }


def _ctx(prec: int = _DEFAULT_PREC, *, up: bool = True):
    ctx = gmpy2.get_context()
    ctx.precision = prec
    ctx.round = RoundUp if up else RoundDown
    return ctx


def mpfr_const(s: str, *, prec: int = _DEFAULT_PREC) -> mpfr:
    return mpfr(s, precision=prec)


def add_up(a: mpfr, b: mpfr) -> mpfr:
    with _ctx(a.precision, up=True):
        return mpfr(a + b)


def add_down(a: mpfr, b: mpfr) -> mpfr:
    with _ctx(a.precision, up=False):
        return mpfr(a + b)


def mul_up(a: mpfr, b: mpfr) -> mpfr:
    with _ctx(a.precision, up=True):
        return mpfr(a * b)


def mul_down(a: mpfr, b: mpfr) -> mpfr:
    with _ctx(a.precision, up=False):
        return mpfr(a * b)


def sqrt_up(a: mpfr) -> mpfr:
    with _ctx(a.precision, up=True):
        return sqrt(a)


def sqrt_lo(a: mpfr) -> mpfr:
    with _ctx(a.precision, up=False):
        return sqrt(a)


def sub_down(a: mpfr, b: mpfr) -> mpfr:
    with _ctx(a.precision, up=False):
        return mpfr(a - b)


def sub_up(a: mpfr, b: mpfr) -> mpfr:
    with _ctx(a.precision, up=True):
        return mpfr(a - b)


def div_up(a: mpfr, b: mpfr) -> mpfr:
    with _ctx(a.precision, up=True):
        return mpfr(a / b)


def div_down(a: mpfr, b: mpfr) -> mpfr:
    with _ctx(a.precision, up=False):
        return mpfr(a / b)


def exp_up(x: mpfr) -> mpfr:
    with _ctx(x.precision, up=True):
        return gmpy2.exp(x)


def exp_down(x: mpfr) -> mpfr:
    with _ctx(x.precision, up=False):
        return gmpy2.exp(x)


def exp_lo(x: mpfr) -> mpfr:
    return exp_down(x)


def exp_hi(x: mpfr) -> mpfr:
    return exp_up(x)


def to_float_hi(x: mpfr) -> float:
    f = float(x)
    if math.isfinite(f):
        return math.nextafter(f, math.inf)
    return f


def to_float_lo(x: mpfr) -> float:
    f = float(x)
    if math.isfinite(f):
        return math.nextafter(f, -math.inf)
    return f


def sqrt8_interval(prec: int = _DEFAULT_PREC) -> tuple[mpfr, mpfr]:
    eight = mpfr(8, precision=prec)
    return sqrt_lo(eight), sqrt_up(eight)


def stokes_floor_interval(
    n: int = 24,
    *,
    prec: int = _DEFAULT_PREC,
) -> tuple[mpfr, mpfr]:
    """
    Phi(0) <= E0 * max_r r exp(-2 nu r T) over dealias shells.
    Stores argmax shell in stokes_floor_meta().
    """
    from ns_exploration.conjectures.l0021_quartic_shell import shell_modes_by_r

    nu = mpfr_const(EXACT.NU, prec=prec)
    E0 = mpfr_const(EXACT.E0, prec=prec)
    T = mpfr_const(EXACT.T, prec=prec)
    neg_two = mpfr(-2, precision=prec)

    by = shell_modes_by_r(n)
    shells = sorted(by.keys())
    best_lo = mpfr(0, precision=prec)
    best_hi = mpfr(0, precision=prec)
    argmax = shells[0] if shells else 0

    for r in shells:
        rf = mpfr(r, precision=prec)
        expo = mul_up(mul_up(neg_two, nu), mul_up(rf, T))
        w_lo = mul_down(rf, exp_lo(expo))
        w_hi = mul_up(rf, exp_hi(expo))
        if w_hi > best_hi:
            best_hi = w_hi
            argmax = r
        if w_lo > best_lo:
            best_lo = w_lo

    floor_lo = mul_down(E0, best_lo)
    floor_hi = mul_up(E0, best_hi)
    _STOKES_META.update({"argmax_shell": argmax, "shells": shells, "n": n})
    return floor_lo, floor_hi


def stokes_floor_meta() -> dict[str, Any]:
    return dict(_STOKES_META)


def i_star_interval(prec: int = _DEFAULT_PREC) -> tuple[mpfr, mpfr]:
    """
    I_* = 2√2 (M_target - Stokes_floor).

    I_star_lo = sqrt8_lo * (M_target_lo - floor_hi).
    """
    sqrt8_lo, sqrt8_hi = sqrt8_interval(prec)
    _, floor_hi = stokes_floor_interval(prec=prec)
    floor_lo, _ = stokes_floor_interval(prec=prec)
    M = mpfr_const(EXACT.M_TARGET, prec=prec)
    gap_lo = sub_down(M, floor_hi)
    gap_hi = sub_up(M, floor_lo)
    I_lo = mul_down(sqrt8_lo, gap_lo)
    I_hi = mul_up(sqrt8_hi, gap_hi)
    return I_lo, I_hi


def omega_T_hi_from_I_term(integral_hi: mpfr, *, prec: int = _DEFAULT_PREC) -> mpfr:
    """Omega_hi = floor_hi + I_term_hi / sqrt8_lo."""
    _, floor_hi = stokes_floor_interval(prec=prec)
    sqrt8_lo, _ = sqrt8_interval(prec)
    return add_up(floor_hi, div_up(integral_hi, sqrt8_lo))
