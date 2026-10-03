"""
Route A: certified shell-integrated upper bound for int_0^T C_term(t) dt.

  int_0^T C_term(t) dt <= sum_r B_r * (1 - exp(-2 nu r T)) / (2 nu r)

with B_r = certified ||M_r||_op upper bound (here Frobenius/MPFR).
"""

from __future__ import annotations

import json
from pathlib import Path

from gmpy2 import mpfr

from ns_exploration.conjectures.l0048_matrixfree_fullsym import all_dealias_radii
from ns_exploration.terminal_weighted.constants import FROZEN
from ns_exploration.terminal_weighted.intervals import (
    add_up,
    div_up,
    exp_down,
    i_star_interval,
    mul_up,
    sqrt_up,
    sub_down,
    _ctx,
)
from ns_exploration.terminal_weighted.shell_decomposition import build_shell_manifest


def shell_integral_certificate(
    n: int = 24,
    prec: int = 200,
    *,
    manifest: dict | None = None,
) -> dict:
    _ctx(prec)
    manifest = manifest or build_shell_manifest(n, prec=prec)
    nu = mpfr(FROZEN.nu, precision=prec)
    T = mpfr(FROZEN.T, precision=prec)
    two = mpfr(2, precision=prec)

    integral_hi = mpfr(0, precision=prec)
    per_shell = []

    for blk in manifest["blocks"]:
        r = mpfr(blk["shell"], precision=prec)
        B_r = mpfr(blk["C_term_hi"], precision=prec)
        expo = mul_up(mul_up(mul_up(mpfr(-2, precision=prec), nu), r), T)
        e_down = exp_down(expo)
        with _ctx(prec, up=True):
            factor_num_hi = mpfr(1) - e_down
        denom = mul_up(mul_up(two, nu), r)
        factor = div_up(factor_num_hi, denom)
        contrib = mul_up(B_r, factor)
        integral_hi = add_up(integral_hi, contrib)
        per_shell.append(
            {
                "shell": int(blk["shell"]),
                "B_r_hi": str(B_r),
                "factor_hi": str(factor),
                "contrib_hi": str(contrib),
            }
        )

    I_lo, I_hi = i_star_interval(prec=prec)
    closes = integral_hi < I_lo
    full_dealias = (
        manifest.get("D_full") == FROZEN.D_full
        and manifest.get("n_shells") == len(all_dealias_radii(n))
    )

    two_sqrt2 = sqrt_up(mpfr(8, precision=prec))
    with _ctx(prec, up=True):
        floor_hi = mpfr(FROZEN.stokes_floor)
    omega_hi = add_up(floor_hi, div_up(integral_hi, two_sqrt2))
    with _ctx(prec, up=False):
        M_target_lo = mpfr(FROZEN.c0008_target_M)

    return {
        "method": "shell_integral_route_A",
        "n": n,
        "evidence_level": "N5",
        "integral_hi": str(integral_hi),
        "I_star_lo": str(I_lo),
        "I_star_hi": str(I_hi),
        "closes_strict": bool(closes),
        "omega_T_hi": str(omega_hi),
        "M_target_lo": str(M_target_lo),
        "closes_c0008": bool(closes and omega_hi < M_target_lo and full_dealias),
        "full_dealias": full_dealias,
        "per_shell": per_shell,
        "manifest_sha256": manifest.get("manifest_sha256"),
        "honesty": (
            "Triangle inequality on shell op-norms; Frobenius not tight. "
            "Full dealias all shells included."
        ),
    }


def save_shell_integral_certificate(cert: dict, path: str | Path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(cert, indent=2), encoding="utf-8")
    return path
