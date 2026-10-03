"""Band cluster residual (L-0073R) prototype with combined B_C one_inf pass."""

from __future__ import annotations

from gmpy2 import mpfr

from ns_exploration.terminal_weighted.cluster_bounds import (
    _shell_factor_hi,
    cluster_residual_vs_legacy_audit,
)
from ns_exploration.terminal_weighted.constants import EXACT
from ns_exploration.terminal_weighted.intervals import add_up, mpfr_const
from ns_exploration.terminal_weighted.opnorm_certified import (
    certified_cluster_one_inf_hi,
    certified_shell_blocks,
)


def run_cluster_residual_audit(
    *,
    n: int = 24,
    radii: tuple[int, ...] = (1, 2, 3),
    prec: int = 128,
) -> dict:
    blocks = certified_shell_blocks(
        n, radii, prec=prec, workers=1, bound_method="one_inf_only"
    )["blocks"]
    B_r = {int(b["shell"]): mpfr(b["C_term_hi"], precision=prec) for b in blocks}

    cluster_bound = certified_cluster_one_inf_hi(n, radii, radii, prec=prec)
    B_C_combined = mpfr(cluster_bound.C_term_hi, precision=prec)

    B_C_sum = mpfr(0, precision=prec)
    for r in radii:
        B_C_sum = add_up(B_C_sum, B_r[r])

    nu = mpfr_const(EXACT.NU, prec=prec)
    T = mpfr_const(EXACT.T, prec=prec)
    phi_r = {
        r: _shell_factor_hi(mpfr(r, precision=prec), nu=nu, T=T, prec=prec) for r in radii
    }

    audit_combined = cluster_residual_vs_legacy_audit(
        radii, B_r, B_C_combined, phi_r, prec=prec
    )
    audit_sum = cluster_residual_vs_legacy_audit(radii, B_r, B_C_sum, phi_r, prec=prec)

    return {
        "n": n,
        "radii": list(radii),
        "bound_method": "one_inf_only",
        "B_C_combined_hi": str(B_C_combined),
        "B_C_triangle_sum_hi": str(B_C_sum),
        "combined_le_sum": float(B_C_combined) <= float(B_C_sum) + 1e-9,
        "residual_combined": audit_combined,
        "residual_triangle_sum": audit_sum,
    }


def run_band123_cluster_residual(*, n: int = 24, prec: int = 128) -> dict:
    return run_cluster_residual_audit(n=n, radii=(1, 2, 3), prec=prec)


def run_band6_cluster_residual(*, n: int = 24, prec: int = 128) -> dict:
    return run_cluster_residual_audit(n=n, radii=(1, 2, 3, 4, 5, 6), prec=prec)


def main() -> dict:
    rep123 = run_band123_cluster_residual()
    print("band123:", rep123)
    rep6 = run_band6_cluster_residual()
    print("band6:", rep6)
    return {"band123": rep123, "band6": rep6}


if __name__ == "__main__":
    main()
