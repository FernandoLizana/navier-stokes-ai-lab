"""Tests for adjoint ascent and C-0001 adversarial update logic."""

from __future__ import annotations

from ns_exploration.conjectures.c0001_adversarial import update_c0001_after_campaign
from ns_exploration.initial_conditions.generators import random_div_free
from ns_exploration.optimization.adjoint_ascent import adjoint_ascent


def test_adjoint_ascent_nondecreasing():
    u0, _ = random_div_free(12, seed=1, energy_target=0.5)
    _u, res = adjoint_ascent(
        u0, nu=0.1, dt=1e-3, t_end=0.015, energy_target=0.5, n_steps=3, step_size=0.05
    )
    assert res.J_final >= res.J_history[0] - 1e-9
    assert res.energy > 0
    assert res.div_l2 < 1e-8


def test_update_tightens_when_not_refuted():
    campaign = {
        "M": 20.0,
        "n": 12,
        "best_J_etd": 6.0,
        "refuted_with_gates": False,
        "refutation_value": None,
    }
    c = update_c0001_after_campaign(campaign, tighten_factor=1.5)
    assert not c.refuted
    assert abs(c.proposed_bound_M - 9.0) < 1e-9
    assert "compact" in c.notes.lower() or "Compactness" in c.notes or "LEMMA" in c.notes


def test_update_refuted_branch():
    campaign = {
        "M": 5.0,
        "n": 12,
        "best_J_etd": 8.0,
        "refuted_with_gates": True,
        "refutation_value": 8.0,
    }
    c = update_c0001_after_campaign(campaign)
    assert c.refuted
    assert c.status == "refuted"
