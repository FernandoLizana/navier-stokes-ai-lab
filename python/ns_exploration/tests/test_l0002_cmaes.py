"""Tests for L-0002 and CMA-lite enstrophy ES."""

from __future__ import annotations

from ns_exploration.conjectures.l0002_viscous_galerkin import (
    bound_for_resolution_l0002,
    explicit_bound_viscous,
    stretch_constant,
)
from ns_exploration.optimization.cmaes_enstrophy import (
    evolve_enstrophy_es,
    theta_to_uhat,
    _pack_divfree_low_modes,
)
from ns_exploration.spectral.fourier_conventions import kinetic_energy_from_hat
from ns_exploration.spectral.leray import divergence_l2


def test_stretch_constant_scales_sqrt_M():
    assert abs(stretch_constant(100) / stretch_constant(25) - 2.0) < 1e-12


def test_l0002_small_M_closes():
    """Tiny mode count: estimate should remain finite on short time."""
    b = explicit_bound_viscous(E0=0.5, t=0.02, K=2.0, M=2, nu=0.1)
    assert not b.finite_time_blowup_of_estimate
    assert b.Omega_t_cap is not None
    assert b.Omega_t_cap >= b.Omega0_cap * 0.1


def test_l0002_large_M_may_fail_honestly():
    b = bound_for_resolution_l0002(16)
    # Either finite cap or honest failure — both OK; must not claim singularity
    assert "Clay" in b.clay_implication or "None" in b.clay_implication
    if b.finite_time_blowup_of_estimate:
        assert b.Omega_t_cap is None
        assert b.t_star_estimate is not None


def test_theta_to_uhat_energy_and_div():
    n = 12
    dofs = _pack_divfree_low_modes(n, k_max=2)
    import numpy as np

    theta = np.zeros(2 * len(dofs))
    theta[0] = 0.3
    theta[1] = -0.1
    uh = theta_to_uhat(theta, n, dofs, 0.5)
    assert abs(kinetic_energy_from_hat(uh) - 0.5) < 1e-10
    assert divergence_l2(uh) < 1e-10


def test_cmaes_short_run_improves_or_stable():
    _u, res = evolve_enstrophy_es(
        n=12,
        k_max=2,
        generations=4,
        population=6,
        mu=2,
        seed=0,
        M_target=1e9,
    )
    assert res.n_evals == 4 * 6
    assert res.best_J > 0
    assert res.history_best[-1] >= res.history_best[0] - 1e-9
