"""Tests for band {1,2,3} exact prototype gate."""

from __future__ import annotations

from ns_exploration.terminal_weighted.exact_prototype import run_band123_prototype


def test_band123_four_path_gate():
    rep = run_band123_prototype()
    assert rep.four_path_ok
    assert rep.max_rel_G_vs_stream <= 1e-12
    assert rep.worst_rel["C_D"] <= 1e-9
    assert rep.interval_enclosed


def test_band123_mpfr_interval_enclosure():
    from ns_exploration.terminal_weighted.interval_tensor import verify_float_G_enclosed
    from ns_exploration.terminal_weighted.constants import FROZEN

    iv = verify_float_G_enclosed(24, (1, 2, 3), FROZEN.T, prec=200)
    assert iv["enclosed"]
    assert iv["n_bad"] == 0


def test_band123_cluster_residual_audit():
    from ns_exploration.terminal_weighted.cluster_residual_prototype import (
        run_band123_cluster_residual,
    )

    rep = run_band123_cluster_residual(prec=128)
    assert rep["combined_le_sum"]
    assert "B_C_combined_hi" in rep
    assert rep["bound_method"] == "one_inf_only"


def test_rational_pol_leray():
    from ns_exploration.terminal_weighted.rational_basis import verify_pol_leray_against_float

    rep = verify_pol_leray_against_float(prec=200)
    assert rep["ok"]


def test_cluster_combined_one_inf_le_triangle():
    from ns_exploration.terminal_weighted.cluster_residual_prototype import (
        run_band123_cluster_residual,
    )

    rep = run_band123_cluster_residual(prec=128)
    assert float(rep["B_C_combined_hi"]) <= float(rep["B_C_triangle_sum_hi"]) + 1e-6


def test_band6_gate_smoke():
    from ns_exploration.terminal_weighted.band_gate import run_band6_gate

    rep = run_band6_gate(n_probe=10)
    assert rep["four_path_ok"]
    assert rep["interval_enclosed"]
    assert rep["rational_interval_ok"]
    assert rep["rational_pol_ok"]


def test_rational_interval_tensor_band123():
    from ns_exploration.terminal_weighted.constants import FROZEN
    from ns_exploration.terminal_weighted.interval_tensor import verify_rational_vs_float_interval

    rep = verify_rational_vs_float_interval(24, (1, 2, 3), FROZEN.T, prec=200)
    assert rep["ok"]


def test_band_majorant_ladder_monotone():
    from ns_exploration.terminal_weighted.n_coverage import band_majorant_ladder

    rep = band_majorant_ladder((1, 2), prec=128)
    assert rep["monotone_non_decreasing"]
    assert len(rep["ladder"]) == 4


def test_frobenius_entry_agg_band123():
    from ns_exploration.terminal_weighted.opnorm_certified import certified_shell_blocks
    from ns_exploration.terminal_weighted.tensor import streaming_terminal_fullsym_M
    from ns_exploration.experiments.sprintA_stretch_tensor import full_sym
    from ns_exploration.terminal_weighted.constants import FROZEN
    from gmpy2 import mpfr

    radii = (1, 2, 3)
    blk = certified_shell_blocks(24, radii, prec=128, bound_method="frobenius")
    b1 = next(b for b in blk["blocks"] if b["shell"] == 1)
    M1, _ = streaming_terminal_fullsym_M(24, radii, FROZEN.T, output_shell=1)
    fro_stream = float((M1 * M1).sum())
    fro_cert = float(mpfr(b1["frobenius_sq_hi"]))
    assert fro_cert >= fro_stream * (1 - 1e-6)
    assert float(mpfr(b1["C_term_one_inf_hi"])) <= float(mpfr(b1["C_term_hi"])) + 1e-6
