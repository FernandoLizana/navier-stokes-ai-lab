"""Tests for L-0073R in cluster manifest (one_inf_only)."""

from __future__ import annotations

from gmpy2 import mpfr

from ns_exploration.terminal_weighted.cluster_bounds import build_cluster_manifest
from ns_exploration.terminal_weighted.shell_decomposition import build_shell_manifest


def test_cluster_manifest_one_inf_uses_l0073r():
    band_radii = tuple(range(1, 7))
    per_shell = build_shell_manifest(
        24, prec=128, radii=band_radii, bound_method="one_inf_only"
    )
    cluster = build_cluster_manifest(
        24,
        radii=band_radii,
        n_clusters=3,
        per_shell_manifest=per_shell,
        bound_method="one_inf_only",
        prec=128,
        workers=1,
    )
    for blk in cluster["blocks"]:
        assert "cluster_residual_L0073R" in blk
        assert "L-0073R" in blk["method"]
        legacy = mpfr(blk["cluster_integral_hi"])
        contrib = mpfr(blk["contrib_hi"])
        assert contrib <= legacy * mpfr("1.0000001")
