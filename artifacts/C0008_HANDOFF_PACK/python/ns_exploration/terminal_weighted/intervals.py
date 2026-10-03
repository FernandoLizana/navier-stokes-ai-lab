"""
Directed rounding utilities for N5 certificates (gmpy2 MPFR).

Uses context.round = RoundUp / RoundDown for certified bounds.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import gmpy2
from gmpy2 import RoundDown, RoundUp, mpfr, sqrt

_DEFAULT_PREC = 200


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


def add_up(a: mpfr, b: mpfr) -> mpfr:
    with _ctx(a.precision, up=True):
        return mpfr(a + b)


def mul_up(a: mpfr, b: mpfr) -> mpfr:
    with _ctx(a.precision, up=True):
        return mpfr(a * b)


def sqrt_up(a: mpfr) -> mpfr:
    with _ctx(a.precision, up=True):
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


def i_star_interval(prec: int = _DEFAULT_PREC) -> tuple[mpfr, mpfr]:
    from ns_exploration.terminal_weighted.constants import FROZEN
    from ns_exploration.terminal_weighted.weights import stokes_floor_hi

    with _ctx(prec, up=True):
        two_sqrt2 = sqrt(mpfr(8))
    with _ctx(prec, up=True):
        M = mpfr(FROZEN.c0008_target_M)
        floor = mpfr(stokes_floor_hi(24))
    with _ctx(prec, up=False):
        gap_lo = mpfr(M - floor)
    with _ctx(prec, up=True):
        gap_hi = mpfr(M - floor)
    with _ctx(prec, up=False):
        I_lo = mpfr(two_sqrt2 * gap_lo)
    with _ctx(prec, up=True):
        I_hi = mpfr(two_sqrt2 * gap_hi)
    return I_lo, I_hi
