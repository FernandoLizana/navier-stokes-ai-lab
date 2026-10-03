"""Tests for L-0003 and candidate polish."""

from __future__ import annotations

from ns_exploration.conjectures.l0001_galerkin_bound import bound_for_resolution
from ns_exploration.conjectures.l0003_uniform_galerkin import (
    bound_for_resolution_l0003,
    explicit_uniform_bound,
)
from ns_exploration.initial_conditions.generators import random_div_free
from ns_exploration.optimization.polish_candidate import polish_candidate


def test_l0003_formula():
    b = explicit_uniform_bound(E0=0.5, K=2.0, M=10, nu=0.1)
    # Ω_eq = 6*10*(0.25)/0.01 = 1500
    assert abs(b.Omega_eq - 1500.0) < 1e-6
    assert abs(b.Omega_uniform_cap - max(4.0 * 0.5, 1500.0)) < 1e-6
    assert b.closes


def test_l0003_beats_l0001_at_n16_long_horizon_conceptually():
    """Uniform L-0003 finite; L-0001 grows in t — at fixed t=0.02 compare available."""
    b3 = bound_for_resolution_l0003(16)
    b1 = bound_for_resolution(16)
    assert b3.Omega_uniform_cap < b1.Omega_t_cap
    assert b3.lemma_id == "L-0003"


def test_l0003_scales_with_M():
    b_small = explicit_uniform_bound(0.5, K=1.0, M=4, nu=0.1)
    b_big = explicit_uniform_bound(0.5, K=1.0, M=16, nu=0.1)
    assert abs(b_big.Omega_eq / b_small.Omega_eq - 4.0) < 1e-12


def test_polish_runs_and_preserves_energy():
    uh, _ = random_div_free(12, seed=3, energy_target=0.5)
    u1, res = polish_candidate(uh, M_target=1e9, n_steps=2, step_size=0.04)
    assert abs(res.energy - 0.5) < 1e-8
    assert res.div_l2 < 1e-8
    assert res.J_etd > 0
    assert u1.shape[-1] == 12


def test_c0002_successor_fields():
    from ns_exploration.conjectures.c0002_enstrophy_bound import ConjectureC0002

    c = ConjectureC0002(proposed_bound_M=20.32, numerical_support_max=13.55)
    c.statement = (
        f"For all divergence-free Fourier fields on the N^3 grid with N<={c.n_max}, "
        f"enstrophy E_ens(T) <= M with M={c.proposed_bound_M}. FINITE-DIMENSIONAL only."
    )
    assert c.id == "C-0002"
    assert c.parent_refuted == "C-0001"
    assert "None" in c.clay_implication
