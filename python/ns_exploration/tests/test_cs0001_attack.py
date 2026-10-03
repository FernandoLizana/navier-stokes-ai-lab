"""Tests for L-0006 and C-S-0001 attack helpers."""

from __future__ import annotations

from ns_exploration.conjectures.cs0001_attack import _shell_ic
from ns_exploration.conjectures.l0005_shell_galerkin import lemma_l0005
from ns_exploration.conjectures.l0006_improved_shell import (
    algebraic_short_time_bound,
    lemma_l0006,
)
from ns_exploration.spectral.leray import divergence_l2


def test_algebraic_closes_for_tiny_shell():
    closes, cap, a, t_star = algebraic_short_time_bound(0.5, 0.01, K=1.0, M=7)
    assert closes
    assert cap is not None and cap > 0
    assert t_star > 0.01


def test_l0006_k2_not_worse_than_l0005_style():
    b6 = lemma_l0006(n=12, k_shell=2, t=0.02)
    b5 = lemma_l0005(n=12, k_shell=2, t=0.02)
    assert b6.best_cap <= b5.Omega_bound + 1e-9


def test_l0006_records_which_bound():
    b = lemma_l0006(12, 3)
    assert b.which_best in ("L-0001", "L-0003", "L-0006-algebraic")
    assert b.evidence_level == "N7"


def test_shell_ic_div_free():
    uh = _shell_ic(12, 4, 0.5, "taylor_green")
    assert divergence_l2(uh) < 1e-10
