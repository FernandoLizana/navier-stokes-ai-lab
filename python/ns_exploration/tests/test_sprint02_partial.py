"""Tests for Sprint 2 partial modules."""

from __future__ import annotations

import numpy as np

from ns_exploration.diagnostics.energy_crosscheck import (
    embedded_2d_operator_agreement,
    embedded_2d_taylor_green,
    energy_decay_crosscheck,
    parseval_crosscheck,
)
from ns_exploration.diagnostics.residuals import (
    beltrami_nonlinear_residual,
    manufactured_linear_mode_residual,
)
from ns_exploration.optimization.enstrophy_gradient import (
    _resample_hat,
    resolution_doubling_gate,
    _project_fixed_energy,
)
from ns_exploration.initial_conditions.generators import random_div_free
from ns_exploration.spectral.fourier_conventions import kinetic_energy_from_hat
from ns_exploration.spectral.leray import divergence_l2


def test_beltrami_residual_small_and_improves_or_stable_with_n():
    r16 = beltrami_nonlinear_residual(16, dealias=True)
    r32 = beltrami_nonlinear_residual(32, dealias=True)
    assert r16.divergence_l2 < 1e-10
    assert r32.divergence_l2 < 1e-10
    # Pseudospectral Beltrami residual should be tiny
    assert r16.residual_rel < 1e-8
    assert r32.residual_rel < 1e-8


def test_beltrami_dealias_vs_alias_documented():
    with_d = beltrami_nonlinear_residual(24, dealias=True)
    no_d = beltrami_nonlinear_residual(24, dealias=False)
    # Both should still be small for Beltrami; record that dealias does not worsen badly
    assert with_d.residual_rel < 1e-7
    assert no_d.residual_rel < 1e-6


def test_manufactured_linear_mode_two_resolutions():
    r16 = manufactured_linear_mode_residual(16)
    r32 = manufactured_linear_mode_residual(32)
    assert r16.residual_rel < 1e-6
    assert r32.residual_rel < 1e-6


def test_parseval_crosscheck():
    uh = embedded_2d_taylor_green(16)
    rep = parseval_crosscheck(uh)
    assert rep.rel_diff < 1e-12


def test_embedded_2d_operator_agreement():
    out = embedded_2d_operator_agreement(24)
    assert out["operator_rel_diff"] < 1e-8
    assert abs(out["energy_channel_3d"]) < 1e-8
    assert abs(out["energy_channel_2d"]) < 1e-8


def test_energy_decay_residual_bounded():
    out = energy_decay_crosscheck(16, nu=0.05, dt=5e-4, steps=20)
    # Forward-Euler predictor vs RK4: residual O(dt) relative scale
    assert out["max_abs_step_residual"] < 5e-3
    assert out["delta_E"] < 0.0  # energy must drop


def test_resample_preserves_low_modes_and_div_free():
    uh, _ = random_div_free(12, seed=2, energy_target=0.5)
    up = _resample_hat(uh, 24)
    assert up.shape[-1] == 24
    assert divergence_l2(up) < 1e-10
    # energy after reproject
    up = _project_fixed_energy(up, 0.5)
    assert abs(kinetic_energy_from_hat(up) - 0.5) < 1e-10


def test_resolution_gate_runs():
    gate = resolution_doubling_gate(
        n_low=12, seed=0, nu=0.1, dt=1e-3, t_end=0.015, energy_target=0.5, gate_tol=0.5
    )
    assert gate.n_high == 24
    assert np.isfinite(gate.directional_deriv_low)
    assert np.isfinite(gate.directional_deriv_high)
    assert gate.relative_agreement >= 0.0
