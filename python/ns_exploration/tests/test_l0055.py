"""Tests for L-0055 greedy SOS feasibility."""

from ns_exploration.conjectures.l0055_greedy_sos_feasibility import lemma_l0055


def test_l0055_greedy_D_and_verdict():
    b = lemma_l0055()
    assert b.greedy_D > 1700
    assert b.extrap_closes_cr0012 is True
    assert b.laptop_feasible is False
    assert b.verdict in ("laptop_infeasible", "cluster_candidate", "marginal_laptop")
    assert b.est_solve_hours > 100


def test_l0055_prefix_monotone():
    b = lemma_l0055()
    assert b.prefix_rows
    Ds = [r.D for r in b.prefix_rows]
    assert Ds == sorted(Ds)
    assert b.prefix_rows[-1].D == b.greedy_D
