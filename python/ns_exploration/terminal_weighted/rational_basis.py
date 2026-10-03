"""
Rational MPFR transversal polarization and Leray projection for integer k.
"""

from __future__ import annotations

import numpy as np
from gmpy2 import mpfr, sqrt

from ns_exploration.conjectures.l0021_quartic_shell import _leray_vec, _pol_basis
from ns_exploration.terminal_weighted.intervals import _ctx


class Mpfr3:
    __slots__ = ("x", "y", "z")

    def __init__(self, x: mpfr, y: mpfr, z: mpfr) -> None:
        self.x, self.y, self.z = x, y, z

    def cross(self, other: Mpfr3, *, up: bool) -> Mpfr3:
        with _ctx(self.x.precision, up=up):
            return Mpfr3(
                mpfr(self.y * other.z - self.z * other.y),
                mpfr(self.z * other.x - self.x * other.z),
                mpfr(self.x * other.y - self.y * other.x),
            )

    def dot(self, other: Mpfr3, *, up: bool) -> mpfr:
        with _ctx(self.x.precision, up=up):
            return mpfr(self.x * other.x + self.y * other.y + self.z * other.z)

    def scale(self, s: mpfr, *, up: bool) -> Mpfr3:
        with _ctx(self.x.precision, up=up):
            return Mpfr3(mpfr(self.x * s), mpfr(self.y * s), mpfr(self.z * s))

    def norm(self, *, up: bool) -> mpfr:
        with _ctx(self.x.precision, up=up):
            return sqrt(self.dot(self, up=up))

    def normalized(self, *, up: bool) -> Mpfr3:
        n = self.norm(up=up)
        if n == 0:
            z = mpfr(0, precision=self.x.precision)
            return Mpfr3(z, z, z)
        with _ctx(self.x.precision, up=up):
            inv = mpfr(1) / n
        return self.scale(inv, up=up)

    def to_float(self) -> np.ndarray:
        return np.array([float(self.x), float(self.y), float(self.z)], dtype=float)


def mpfr_from_int3(k: tuple[int, int, int], *, prec: int) -> Mpfr3:
    return Mpfr3(
        mpfr(k[0], precision=prec),
        mpfr(k[1], precision=prec),
        mpfr(k[2], precision=prec),
    )


def mpfr_pol_basis(k: tuple[int, int, int], *, prec: int) -> tuple[np.ndarray, np.ndarray]:
    kf = mpfr_from_int3(k, prec=prec)
    ex = Mpfr3(mpfr(1, precision=prec), mpfr(0, precision=prec), mpfr(0, precision=prec))
    ey = Mpfr3(mpfr(0, precision=prec), mpfr(1, precision=prec), mpfr(0, precision=prec))
    if abs(k[0]) < 0.9:
        a = kf.cross(ex, up=True).normalized(up=True)
    else:
        a = kf.cross(ey, up=True).normalized(up=True)
    b = kf.cross(a, up=True).normalized(up=True)
    return a.to_float(), b.to_float()


def _aligned_rel(pref: np.ndarray, pmp: np.ndarray) -> float:
    scale = max(1.0, np.linalg.norm(pref), np.linalg.norm(pmp))
    return min(
        np.linalg.norm(pref - pmp) / scale,
        np.linalg.norm(pref + pmp) / scale,
    )


def mpfr_leray(s: tuple[int, int, int], v: np.ndarray, *, prec: int) -> np.ndarray:
    sf = mpfr_from_int3(s, prec=prec)
    vf = Mpfr3(
        mpfr(float(np.real(v[0])), precision=prec),
        mpfr(float(np.real(v[1])), precision=prec),
        mpfr(float(np.real(v[2])), precision=prec),
    )
    # Leray on real part; imaginary part projected separately if present
    vi = Mpfr3(
        mpfr(float(np.imag(v[0])), precision=prec),
        mpfr(float(np.imag(v[1])), precision=prec),
        mpfr(float(np.imag(v[2])), precision=prec),
    )
    sn2 = sf.dot(sf, up=True)
    if sn2 == 0:
        return np.zeros(3, dtype=np.complex128)
    with _ctx(prec, up=True):
        coeff_r = mpfr(sf.dot(vf, up=True) / sn2)
        coeff_i = mpfr(sf.dot(vi, up=True) / sn2)
    proj_r = sf.scale(coeff_r, up=True)
    proj_i = sf.scale(coeff_i, up=True)
    out = np.array(
        [
            complex(float(vf.x - proj_r.x), float(vi.x - proj_i.x)),
            complex(float(vf.y - proj_r.y), float(vi.y - proj_i.y)),
            complex(float(vf.z - proj_r.z), float(vi.z - proj_i.z)),
        ],
        dtype=np.complex128,
    )
    return out


def mpfr_vdot_real(a: np.ndarray, b: np.ndarray, *, prec: int, up: bool) -> mpfr:
    """Real part of conj(a)·b with directed MPFR."""
    with _ctx(prec, up=up):
        acc = mpfr(0)
        for i in range(3):
            ar = float(np.real(a[i]))
            ai = float(np.imag(a[i]))
            br = float(np.real(b[i]))
            bi = float(np.imag(b[i]))
            acc += mpfr(ar * br + ai * bi)
        return acc


def mpfr_triad_coef(va: np.ndarray, kq: tuple[int, int, int], *, prec: int, up: bool) -> complex:
    """1j * va · kq with MPFR real/imag parts (kq exact integer)."""
    kf = mpfr_from_int3(kq, prec=prec)
    with _ctx(prec, up=up):
        re = mpfr(
            float(np.real(va[0])) * float(kf.x)
            + float(np.real(va[1])) * float(kf.y)
            + float(np.real(va[2])) * float(kf.z)
        )
        im = mpfr(
            float(np.imag(va[0])) * float(kf.x)
            + float(np.imag(va[1])) * float(kf.y)
            + float(np.imag(va[2])) * float(kf.z)
        )
    return 1j * complex(float(re), float(im))


def verify_pol_leray_against_float(
    modes: list[tuple[int, int, int]] | None = None,
    *,
    prec: int = 200,
) -> dict:
    if modes is None:
        modes = [
            (1, 0, 0),
            (1, 1, 0),
            (1, 1, 1),
            (0, 1, 1),
            (2, 0, 0),
            (2, 1, 0),
        ]
    worst_pol = 0.0
    worst_leray = 0.0
    for k in modes:
        a_mp, b_mp = mpfr_pol_basis(k, prec=prec)
        pols = _pol_basis(np.array(k, float))
        for pref in pols[:2]:
            worst_pol = max(
                worst_pol,
                min(_aligned_rel(pref, a_mp), _aligned_rel(pref, b_mp)),
            )
        s = (2, 1, 0)
        v = np.array([0.3, -0.2, 0.7])
        lf = _leray_vec(np.array(s, float), v)
        lm = mpfr_leray(s, v, prec=prec)
        scale = max(1.0, np.linalg.norm(lf), np.linalg.norm(lm))
        worst_leray = max(worst_leray, np.linalg.norm(lf - lm) / scale)
    return {
        "modes_tested": len(modes),
        "worst_pol_rel": worst_pol,
        "worst_leray_rel": worst_leray,
        "ok": worst_pol <= 1e-12 and worst_leray <= 1e-12,
        "prec": prec,
    }
