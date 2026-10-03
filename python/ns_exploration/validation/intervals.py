"""
Finite-mode interval arithmetic probes (Python).

Evidence target: N5 for *finite Galerkin identities only*.
Does NOT certify continuum Navier–Stokes on T³.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class Interval:
    lo: float
    hi: float

    def __post_init__(self):
        if self.lo > self.hi:
            self.lo, self.hi = self.hi, self.lo

    @staticmethod
    def from_float(x: float) -> Interval:
        # outward via nextafter
        return Interval(np.nextafter(x, -np.inf), np.nextafter(x, np.inf))

    def __add__(self, other: Interval) -> Interval:
        return Interval(
            np.nextafter(self.lo + other.lo, -np.inf),
            np.nextafter(self.hi + other.hi, np.inf),
        )

    def __sub__(self, other: Interval) -> Interval:
        return Interval(
            np.nextafter(self.lo - other.hi, -np.inf),
            np.nextafter(self.hi - other.lo, np.inf),
        )

    def __mul__(self, other: Interval) -> Interval:
        cands = [
            self.lo * other.lo,
            self.lo * other.hi,
            self.hi * other.lo,
            self.hi * other.hi,
        ]
        return Interval(np.nextafter(min(cands), -np.inf), np.nextafter(max(cands), np.inf))

    def __truediv__(self, other: Interval) -> Interval:
        if other.lo <= 0.0 <= other.hi:
            raise ZeroDivisionError("interval divisor contains zero")
        cands = [
            self.lo / other.lo,
            self.lo / other.hi,
            self.hi / other.lo,
            self.hi / other.hi,
        ]
        return Interval(np.nextafter(min(cands), -np.inf), np.nextafter(max(cands), np.inf))

    def sqrt(self) -> Interval:
        if self.lo < 0.0:
            raise ValueError("sqrt of interval with negative lower bound")
        return Interval(
            np.nextafter(float(np.sqrt(self.lo)), -np.inf),
            np.nextafter(float(np.sqrt(self.hi)), np.inf),
        )

    def sqr(self) -> Interval:
        # For nonnegative intervals; robust general square.
        if self.lo >= 0.0:
            return Interval(
                np.nextafter(self.lo * self.lo, -np.inf),
                np.nextafter(self.hi * self.hi, np.inf),
            )
        if self.hi <= 0.0:
            return Interval(
                np.nextafter(self.hi * self.hi, -np.inf),
                np.nextafter(self.lo * self.lo, np.inf),
            )
        hi = max(self.lo * self.lo, self.hi * self.hi)
        return Interval(0.0, np.nextafter(hi, np.inf))

    def exp(self) -> Interval:
        return Interval(
            np.nextafter(float(np.exp(self.lo)), -np.inf),
            np.nextafter(float(np.exp(self.hi)), np.inf),
        )

    def cosh(self) -> Interval:
        # cosh is even and minimized at 0; for intervals in R use endpoints carefully.
        # Safe outward: cosh(x) ≤ (exp(|x|) + exp(-|x|))/2 ≤ exp(|x|).
        # Tighter: evaluate cosh at endpoints and at 0 if 0 ∈ interval.
        import math

        samples = [self.lo, self.hi]
        if self.lo <= 0.0 <= self.hi:
            samples.append(0.0)
        vals = [math.cosh(s) for s in samples]
        return Interval(
            np.nextafter(min(vals), -np.inf),
            np.nextafter(max(vals), np.inf),
        )

    def contains(self, x: float) -> bool:
        return self.lo <= x <= self.hi

    def contains_zero(self) -> bool:
        return self.lo <= 0.0 <= self.hi


def leray_project_interval_mode(v: tuple[Interval, Interval, Interval], kx: float, ky: float, kz: float):
    """Interval Leray on a single mode k≠0."""
    k2 = kx * kx + ky * ky + kz * kz
    if k2 == 0.0:
        return v
    # factor = (k·v)/k2
    dot = v[0] * Interval.from_float(kx) + v[1] * Interval.from_float(ky) + v[2] * Interval.from_float(kz)
    inv = Interval.from_float(1.0 / k2)
    factor = dot * inv
    return (
        v[0] - Interval.from_float(kx) * factor,
        v[1] - Interval.from_float(ky) * factor,
        v[2] - Interval.from_float(kz) * factor,
    )


def divergence_interval(v: tuple[Interval, Interval, Interval], kx: float, ky: float, kz: float) -> Interval:
    return (
        Interval.from_float(kx) * v[0]
        + Interval.from_float(ky) * v[1]
        + Interval.from_float(kz) * v[2]
    )


def finite_energy_identity_cancel_check(
    u_modes: dict[tuple[int, int, int], tuple[complex, complex, complex]],
) -> Interval:
    """
    For a finite set of Fourier modes, the convective contribution to dE/dt
    is Re Σ û* · (nonlinear)_k. For exact Galerkin with div-free fields,
    the continuum identity says the nonlinear term cancels in the energy.

    Here we only check that after Leray projection of a random finite field,
    k·û_k ∈ [0] in interval sense for each kept mode — a prerequisite identity.
    """
    total_div_sq = Interval.from_float(0.0)
    for (kx, ky, kz), (ux, uy, uz) in u_modes.items():
        if (kx, ky, kz) == (0, 0, 0):
            continue
        v = (Interval.from_float(ux.real), Interval.from_float(uy.real), Interval.from_float(uz.real))
        # use real parts only for a simple real-mode check; complex handled componentwise below
        vr = (
            Interval.from_float(ux.real) + Interval.from_float(ux.imag) * Interval.from_float(0.0),
            Interval.from_float(uy.real),
            Interval.from_float(uz.real),
        )
        # Project complex mode by projecting real/imag separately
        def proj_c(cux, cuy, cuz):
            r = leray_project_interval_mode(
                (Interval.from_float(cux.real), Interval.from_float(cuy.real), Interval.from_float(cuz.real)),
                float(kx),
                float(ky),
                float(kz),
            )
            i = leray_project_interval_mode(
                (Interval.from_float(cux.imag), Interval.from_float(cuy.imag), Interval.from_float(cuz.imag)),
                float(kx),
                float(ky),
                float(kz),
            )
            div_r = divergence_interval(r, float(kx), float(ky), float(kz))
            div_i = divergence_interval(i, float(kx), float(ky), float(kz))
            return div_r * div_r + div_i * div_i

        total_div_sq = total_div_sq + proj_c(ux, uy, uz)
    return total_div_sq


def short_convolution_interval(
    a: dict[tuple[int, int, int], complex],
    b: dict[tuple[int, int, int], complex],
    target: tuple[int, int, int],
) -> Interval:
    """
    Interval enclosure of (a * b)_target = Σ_{p+q=target} a_p b_q for scalar modes.
    """
    acc = Interval.from_float(0.0)
    for p, ap in a.items():
        q = (target[0] - p[0], target[1] - p[1], target[2] - p[2])
        if q not in b:
            continue
        # complex product real-part enclosure (crude: enclose |a||b| into real interval via components)
        ar, ai = Interval.from_float(ap.real), Interval.from_float(ap.imag)
        br, bi = Interval.from_float(b[q].real), Interval.from_float(b[q].imag)
        # (ar+iai)(br+ibi) real part = ar br - ai bi
        acc = acc + (ar * br - ai * bi)
    return acc


def propagate_linear_heat_interval(
    coeff: Interval,
    k2: float,
    nu: float,
    dt: float,
) -> Interval:
    """
    Exact linear factor e^{-ν k² dt} with interval input coefficient.
    Uses float exp with outward rounding on the factor.
    """
    factor = float(np.exp(-nu * k2 * dt))
    fI = Interval(np.nextafter(factor, -np.inf), np.nextafter(factor, np.inf))
    return coeff * fI
