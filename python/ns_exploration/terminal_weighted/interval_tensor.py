"""
MPFR-directed terminal tensor intervals for small bands.

Phase 7: optional rational Leray + directed inner products (mpfr_leray, mpfr_vdot_real).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from gmpy2 import mpfr

from ns_exploration.conjectures.l0021_quartic_shell import _leray_vec, shell_modes_by_r
from ns_exploration.experiments.sprintB_physical_tensor import build_hermitian_basis
from ns_exploration.terminal_weighted.constants import EXACT
from ns_exploration.terminal_weighted.intervals import (
    add_down,
    add_up,
    div_down,
    div_up,
    exp_lo,
    exp_up,
    mpfr_const,
    mul_down,
    mul_up,
    sub_down,
    sub_up,
    _ctx,
)
from ns_exploration.terminal_weighted.rational_basis import mpfr_leray, mpfr_vdot_real
from ns_exploration.terminal_weighted.tensor import build_terminal_G


def _mpfr_vec(v: np.ndarray, *, prec: int, up: bool) -> tuple[mpfr, mpfr, mpfr]:
    with _ctx(prec, up=up):
        return mpfr(v[0]), mpfr(v[1]), mpfr(v[2])


def _dot3(a: tuple[mpfr, mpfr, mpfr], b: tuple[mpfr, mpfr, mpfr], *, up: bool) -> mpfr:
    with _ctx(a[0].precision, up=up):
        return mpfr(a[0] * b[0] + a[1] * b[1] + a[2] * b[2])


def _terminal_weight_lo_hi(r: int, t: mpfr, *, nu: mpfr, T: mpfr, prec: int) -> tuple[mpfr, mpfr]:
    rf = mpfr(r, precision=prec)
    neg2 = mpfr(-2, precision=prec)
    expo_lo = mul_down(mul_down(mul_down(neg2, nu), rf), sub_down(T, t))
    expo_hi = mul_up(mul_up(mul_up(neg2, nu), rf), sub_up(T, t))
    return mul_down(rf, exp_lo(expo_lo)), mul_up(rf, exp_up(expo_hi))


def _build_G_directed(
    n: int,
    radii: tuple[int, ...],
    t: float,
    *,
    prec: int,
    up: bool,
    rational: bool = False,
) -> np.ndarray:
    """Single-directed MPFR rebuild of build_terminal_G."""
    b = build_hermitian_basis(n, radii)
    by = shell_modes_by_r(n)
    out_set: set[tuple[int, int, int]] = set()
    for r in radii:
        if r in by:
            out_set.update(by[r])
    nu = mpfr_const(EXACT.NU, prec=prec)
    T = mpfr_const(EXACT.T, prec=prec)
    tm = mpfr(t, precision=prec)
    D = b.D
    G = np.zeros((D, D, D), dtype=float)
    for a in range(D):
        for bb in range(D):
            Nloc: dict[tuple[int, int, int], np.ndarray] = {}
            for kp, va in b.psis[a]:
                for kq, vb in b.psis[bb]:
                    s = (kp[0] + kq[0], kp[1] + kq[1], kp[2] + kq[2])
                    if s not in out_set:
                        continue
                    coef = 1j * np.dot(va, np.array(kq, float))
                    if s not in Nloc:
                        Nloc[s] = np.zeros(3, dtype=np.complex128)
                    Nloc[s] += coef * vb
            if not Nloc:
                continue
            for s, vec in list(Nloc.items()):
                if rational:
                    Nloc[s] = -mpfr_leray(s, vec, prec=prec)
                else:
                    Nloc[s] = -_leray_vec(np.array(s, float), vec)
            for m in range(D):
                acc = mpfr(0, precision=prec)
                for km, vm in b.psis[m]:
                    if km not in Nloc:
                        continue
                    r2 = km[0] * km[0] + km[1] * km[1] + km[2] * km[2]
                    w_lo, w_hi = _terminal_weight_lo_hi(r2, tm, nu=nu, T=T, prec=prec)
                    w = w_hi if up else w_lo
                    if rational:
                        re = mpfr_vdot_real(vm, Nloc[km], prec=prec, up=up)
                    else:
                        re = mpfr(float(np.real(np.vdot(vm, Nloc[km]))), precision=prec)
                    term = mul_up(w, re) if up else mul_down(w, re)
                    acc = add_up(acc, term) if up else add_down(acc, term)
                G[m, a, bb] = float(acc)
    return G


def build_terminal_G_interval_band(
    n: int,
    radii: tuple[int, ...],
    t: float,
    *,
    prec: int = 200,
    rational: bool = False,
) -> tuple[np.ndarray, np.ndarray]:
    """Return (G_lo, G_hi) arrays enclosing the true terminal tensor."""
    G_lo = _build_G_directed(n, radii, t, prec=prec, up=False, rational=rational)
    G_hi = _build_G_directed(n, radii, t, prec=prec, up=True, rational=rational)
    return G_lo, G_hi


def verify_float_G_enclosed(
    n: int,
    radii: tuple[int, ...],
    t: float,
    *,
    prec: int = 200,
    rational: bool = False,
) -> dict:
    Gf, b = build_terminal_G(n, radii, t)
    G_lo, G_hi = build_terminal_G_interval_band(n, radii, t, prec=prec, rational=rational)
    worst = 0.0
    n_bad = 0
    for m in range(b.D):
        for a in range(b.D):
            for bb in range(b.D):
                fv = Gf[m, a, bb]
                if fv == 0.0 and G_lo[m, a, bb] == 0.0 and G_hi[m, a, bb] == 0.0:
                    continue
                if fv < G_lo[m, a, bb] - 1e-12 or fv > G_hi[m, a, bb] + 1e-12:
                    n_bad += 1
                    worst = max(worst, max(G_lo[m, a, bb] - fv, fv - G_hi[m, a, bb], 0.0))
    return {
        "n": n,
        "radii": list(radii),
        "D": b.D,
        "n_bad": n_bad,
        "worst_violation": worst,
        "enclosed": n_bad == 0,
        "max_width": float(np.max(G_hi - G_lo)),
        "rational": rational,
    }


def verify_rational_vs_float_interval(
    n: int,
    radii: tuple[int, ...],
    t: float,
    *,
    prec: int = 200,
) -> dict:
    """Rational MPFR rebuild should agree with float rebuild within interval width."""
    float_iv = verify_float_G_enclosed(n, radii, t, prec=prec, rational=False)
    rat_iv = verify_float_G_enclosed(n, radii, t, prec=prec, rational=True)
    Gf_lo, Gf_hi = build_terminal_G_interval_band(n, radii, t, prec=prec, rational=False)
    Gr_lo, Gr_hi = build_terminal_G_interval_band(n, radii, t, prec=prec, rational=True)
    max_mid_gap = float(np.max(np.abs((Gr_lo + Gr_hi) / 2 - (Gf_lo + Gf_hi) / 2)))
    return {
        "float_enclosed": float_iv["enclosed"],
        "rational_enclosed": rat_iv["enclosed"],
        "max_midpoint_gap": max_mid_gap,
        "ok": float_iv["enclosed"] and rat_iv["enclosed"] and max_mid_gap <= 1e-9,
    }
