"""Tests for L-0030 triad ‖N‖ bound."""

from __future__ import annotations

from ns_exploration.conjectures.l0030_triad_N_bound import (
    N_triad_bound,
    N_young_bound,
    lemma_l0030,
)
from ns_exploration.validation.l0030_certificate import (
    build_l0030_certificate,
    verify_l0030_certificate,
)


def test_l0030_triad_beats_young_near_equipartition():
    b = lemma_l0030(n=24, n_grid=11)
    assert b.R_star == 3148
    assert b.R_star < b.M_modes
    assert b.N_triad_at_equipartition < b.N_young_at_E0
    assert b.C_from_Rstar > b.C_dagger
    assert b.low_slab_hybrid_worst > b.c0007_M
    assert not b.closes_c0007
    assert N_triad_bound(0.5, 0.5, b.R_star) == b.N_triad_at_equipartition
    assert N_young_bound(0.5, b.rho_star) == b.N_young_at_E0


def test_l0030_certificate():
    cert = build_l0030_certificate(n=24)
    ok, checks = verify_l0030_certificate(cert.as_dict())
    assert ok and all(checks.values())
