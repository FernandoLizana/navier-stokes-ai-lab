"""Tests for C-S-0001 / L-0005."""

from __future__ import annotations

import numpy as np

from ns_exploration.conjectures.cs0001_shell import ConjectureCS0001
from ns_exploration.conjectures.l0005_shell_galerkin import (
    lemma_l0005,
    numerical_check_below_l0005,
    project_to_shell,
    shell_mask,
)
from ns_exploration.initial_conditions.generators import random_div_free, taylor_green
from ns_exploration.optimization.enstrophy_gradient import _project_fixed_energy
from ns_exploration.spectral.leray import divergence_l2


def test_project_to_shell_kills_high_modes():
    uh, _ = random_div_free(16, seed=0, k_peak=6, energy_target=0.5)
    us = project_to_shell(uh, k_shell=3)
    mask = shell_mask(16, 3)
    assert float(np.sqrt(np.sum(np.abs(us[:, ~mask]) ** 2))) < 1e-14
    assert divergence_l2(us) < 1e-10


def test_l0005_bound_finite_and_ordered():
    b3 = lemma_l0005(n=12, k_shell=3)
    b4 = lemma_l0005(n=12, k_shell=4)
    assert b3.Omega_bound < b4.Omega_bound
    assert b3.Omega_bound < 100  # k=3 shell ~45
    assert b4.evidence_level == "N7"


def test_truncated_evolution_stays_in_shell_and_under_bound():
    uh, _ = taylor_green(12)
    uh = _project_fixed_energy(project_to_shell(uh, 4), 0.5)
    check = numerical_check_below_l0005(uh, k_shell=4)
    assert check["energy_outside_shell"] < 1e-12
    assert check["below_bound"]


def test_cs0001_statement_mentions_shell():
    c = ConjectureCS0001(proposed_bound_M=18.0, k_shell=4)
    c.ensure_statement()
    assert "|k|" in c.statement
    assert "None" in c.clay_implication
