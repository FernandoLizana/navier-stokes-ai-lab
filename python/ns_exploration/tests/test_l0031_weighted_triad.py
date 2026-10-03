"""Tests for L-0031 weighted triad."""

from ns_exploration.conjectures.l0031_weighted_triad import lemma_l0031, weighted_triad_Rstar


def test_weighted_triad_rp_rq_smaller_than_unweighted():
    R_w, R_u, _ = weighted_triad_Rstar(24, weight="inv_sqrt_rp_rq")
    assert R_w < R_u
    b = lemma_l0031(weight="inv_sqrt_rp_rq")
    assert b.C_from_weighted_R < b.C_from_Rstar
    assert b.C_from_weighted_R > b.C_dagger  # still techo vs C_dagger
