"""Tests for L-0001 and stress helpers."""

from __future__ import annotations

import math

from ns_exploration.conjectures.l0001_galerkin_bound import (
    bound_for_resolution,
    explicit_enstrophy_bound,
    retained_mode_stats,
)
from ns_exploration.initial_conditions.generators import random_div_free
from ns_exploration.optimization.adjoint_ascent import enstrophy_integrator
from ns_exploration.optimization.enstrophy_gradient import _project_fixed_energy


def test_retained_modes_dealias_smaller_than_full():
    K_d, M_d = retained_mode_stats(16, dealias=True)
    K_f, M_f = retained_mode_stats(16, dealias=False)
    assert M_d < M_f
    assert K_d <= K_f
    assert K_d > 0


def test_explicit_bound_monotonic_in_t():
    b1 = explicit_enstrophy_bound(0.5, 0.01, K=5.0, M=100)
    b2 = explicit_enstrophy_bound(0.5, 0.02, K=5.0, M=100)
    assert b2.Omega_t_cap > b1.Omega_t_cap
    assert b1.Omega0_cap == 5.0**2 * 0.5


def test_numerical_enstrophy_below_L0001():
    """N2 check: a random trajectory stays under the crude analytical ceiling."""
    n = 12
    b = bound_for_resolution(n)
    uh, _ = random_div_free(n, seed=0, energy_target=0.5)
    uh = _project_fixed_energy(uh, 0.5)
    j = enstrophy_integrator(uh, 0.1, 1e-3, 0.02, "etd_rk2")
    assert j < b.Omega_t_cap
    assert b.Omega_t_cap > 1.0  # crude bound should be loose
    assert math.isfinite(b.Omega_t_cap)


def test_bound_for_n16_finite():
    b = bound_for_resolution(16)
    assert b.lemma_id == "L-0001"
    assert b.evidence_level == "N7"
    assert "Clay" in b.clay_implication or "None" in b.clay_implication
