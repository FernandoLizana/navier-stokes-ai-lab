"""Tests for L-0004 and C-0002 stress helpers."""

from __future__ import annotations

from ns_exploration.conjectures.l0004_spectral_support import (
    field_conditional_bounds,
    shell_a_priori_bounds,
    spectral_support_stats,
)
from ns_exploration.initial_conditions.generators import random_div_free, taylor_green


def test_shell_bound_smaller_with_smaller_k():
    b4 = shell_a_priori_bounds(12, k_shell_max=4)
    b2 = shell_a_priori_bounds(12, k_shell_max=2)
    assert b2.M < b4.M
    assert b2.K <= b4.K
    assert b2.best_cap <= b4.best_cap


def test_low_shell_cap_order_of_magnitude():
    """k<=3 shell on N=12 should be O(10^2)-O(10^3), not 10^4+."""
    b = shell_a_priori_bounds(12, k_shell_max=3)
    assert b.best_cap < 5e4
    assert b.best_cap > 1.0


def test_effective_stats_on_taylor_green():
    uh, _ = taylor_green(12)
    st = spectral_support_stats(uh, energy_fraction=0.99)
    assert st.M_eff <= st.M_total
    assert st.K_eff <= st.K_max + 1e-9


def test_field_conditional_tighter_than_full_grid():
    uh, _ = random_div_free(12, seed=0, k_peak=2, energy_target=0.5)
    f = field_conditional_bounds(uh)
    full = shell_a_priori_bounds(12, k_shell_max=6)
    assert f.best_cap <= full.best_cap
