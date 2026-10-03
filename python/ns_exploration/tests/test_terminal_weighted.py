"""Tests for C-0008 terminal-weighted first delivery (spec §21 subset)."""

from __future__ import annotations

import math

import numpy as np

from ns_exploration.conjectures.l0045_streaming_fullsym import streaming_C_fullsym
from ns_exploration.conjectures.l0048_fullsym_techo import C_FULLSYM_FULL_DEALIAS
from ns_exploration.terminal_weighted.constants import FROZEN
from ns_exploration.terminal_weighted.tensor import (
    terminal_fullsym_bound,
    terminal_sparse_at_T,
    validate_direct_vs_tensor,
)
from ns_exploration.terminal_weighted.weights import (
    stokes_floor_hi,
    terminal_weight,
    terminal_weight_at_T_equals_r,
)


BAND_18 = tuple(range(1, 9))


def test_terminal_weight_at_T_matches_A():
    assert terminal_weight(147, FROZEN.T, FROZEN.nu, FROZEN.T) == 147.0
    assert terminal_weight_at_T_equals_r(134.0) == 134.0


def test_terminal_phi_at_0_matches_stokes_formula():
    hi = stokes_floor_hi(24)
    assert abs(hi - FROZEN.stokes_floor) < 1e-6


def test_terminal_tensor_matches_direct():
    val = validate_direct_vs_tensor(24, tuple(range(1, 7)), FROZEN.T, n_probe=20, rtol=1e-11)
    assert val["ok"], val


def test_l0048_regression_at_T_band():
    tw = terminal_fullsym_bound(24, BAND_18, FROZEN.T)
    l45 = streaming_C_fullsym(24, BAND_18)
    rel = abs(tw["C_term_ub"] - l45["C_fullsym"]) / l45["C_fullsym"]
    assert rel < 1e-10


def test_terminal_weight_monotone_in_t_for_fixed_r():
    rs = [1.0, 50.0, 147.0]
    for r in rs:
        ws = [terminal_weight(r, t, FROZEN.nu, FROZEN.T) for t in np.linspace(0, FROZEN.T, 20)]
        assert ws[-1] == r
        assert all(ws[i] <= ws[i + 1] + 1e-15 for i in range(len(ws) - 1))


def test_integral_threshold_constants():
    assert abs(FROZEN.delta_M - (FROZEN.c0008_target_M - FROZEN.stokes_floor)) < 1e-12
    expected_I = 2.0 * math.sqrt(2.0) * FROZEN.delta_M
    assert abs(expected_I - FROZEN.integral_threshold) < 1e-9
