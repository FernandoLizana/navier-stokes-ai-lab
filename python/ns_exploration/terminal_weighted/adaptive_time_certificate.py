"""Adaptive temporal cover (Route B) with Lipschitz shell bounds."""

from __future__ import annotations

import json
from pathlib import Path

from gmpy2 import mpfr

from ns_exploration.terminal_weighted.constants import FROZEN
from ns_exploration.terminal_weighted.intervals import (
    add_up,
    div_up,
    exp_down,
    exp_up,
    i_star_interval,
    mul_up,
    sqrt_up,
    _ctx,
)
from ns_exploration.terminal_weighted.shell_decomposition import build_shell_manifest


def cell_integral_hi(
    U_j: mpfr,
    L_j: mpfr,
    h_j: mpfr,
    two_sqrt2: mpfr,
) -> mpfr:
    """∫_{cell} 2√2 (U_j + |t-m| L_j) dt <= 2√2 (h U + h²/4 L)."""
    term1 = mul_up(h_j, U_j)
    h2 = mul_up(h_j, h_j)
    h2 = div_up(h2, mpfr(4, precision=h_j.precision))
    term2 = mul_up(h2, L_j)
    return mul_up(two_sqrt2, add_up(term1, term2))


def adaptive_lipschitz_certificate(
    n: int = 24,
    prec: int = 200,
    initial_cells: int = 4,
    max_cells: int = 64,
    min_width: float = 0.001,
    manifest: dict | None = None,
) -> dict:
    _ctx(prec)
    manifest = manifest or build_shell_manifest(n, prec=prec)
    T = mpfr(FROZEN.T, precision=prec)
    nu = mpfr(FROZEN.nu, precision=prec)
    two_sqrt2 = sqrt_up(mpfr(8, precision=prec))
    min_w = mpfr(min_width, precision=prec)
    I_lo, _ = i_star_interval(prec=prec)

    blocks = []
    for blk in manifest["blocks"]:
        r = mpfr(blk["shell"], precision=prec)
        B_r = mpfr(blk["C_term_hi"], precision=prec)
        L_r = mul_up(mul_up(mul_up(mpfr(2, precision=prec), nu), r), r)
        L_r = mul_up(L_r, B_r)
        blocks.append((r, B_r, L_r))

    queue: list[tuple[mpfr, mpfr]] = []
    h0 = div_up(T, mpfr(initial_cells, precision=prec))
    for i in range(initial_cells):
        a = mul_up(h0, mpfr(i, precision=prec))
        b = mul_up(h0, mpfr(i + 1, precision=prec))
        queue.append((a, b))

    total_hi = mpfr(0, precision=prec)
    cells_out: list[dict] = []
    n_cells = 0
    coverage_complete = False

    from ns_exploration.terminal_weighted.intervals import sub_down

    while queue and n_cells < max_cells:
        a, b = queue.pop(0)
        with _ctx(prec, up=True):
            h = mpfr(b - a)
        m = div_up(add_up(a, b), mpfr(2, precision=prec))

        U_j = mpfr(0, precision=prec)
        L_j = mpfr(0, precision=prec)
        for r, B_r, L_r in blocks:
            expo = mul_up(mul_up(mul_up(mpfr(-2, precision=prec), nu), r), sub_down(T, m))
            phi = exp_up(expo)
            U_j = add_up(U_j, mul_up(phi, B_r))
            L_j = add_up(L_j, mul_up(L_r, phi))

        if h > min_w and n_cells + 2 <= max_cells and L_j > mul_up(U_j, mpfr(4, precision=prec)):
            queue.insert(0, (m, b))
            queue.insert(0, (a, m))
            continue

        local = cell_integral_hi(U_j, L_j, h, two_sqrt2)
        total_hi = add_up(total_hi, local)
        cells_out.append(
            {"a": str(a), "b": str(b), "U_j_hi": str(U_j), "L_j_hi": str(L_j), "contrib_hi": str(local)}
        )
        n_cells += 1

    coverage_complete = len(queue) == 0 and n_cells > 0
    if queue:
        coverage_complete = False

    with _ctx(prec, up=True):
        floor_hi = mpfr(FROZEN.stokes_floor)
    from ns_exploration.terminal_weighted.intervals import omega_T_hi_from_I_term

    omega_hi = omega_T_hi_from_I_term(total_hi, prec=prec)

    return {
        "method": "adaptive_lipschitz_route_B",
        "n_cells": n_cells,
        "integral_hi": str(total_hi),
        "I_star_lo": str(I_lo),
        "closes_strict": bool(total_hi < I_lo),
        "omega_T_hi": str(omega_hi),
        "cells": cells_out,
        "evidence_level": "N4",
        "max_cells_hit": n_cells >= max_cells,
        "coverage_complete": coverage_complete,
        "min_width": min_width,
        "note": "Partial subtotal never reported as full integral unless queue empty.",
    }
