"""Tests for L-0051 structured ladder."""

from ns_exploration.conjectures.l0051_structured_ladder import lemma_l0051


def test_l0051_brute_sos_is_techo():
    b = lemma_l0051()
    assert b.ratio_median_large > 0.5
    assert b.C_ub_sos_extrapolated_full > b.C_dagger
    assert b.brute_sos_closes_all_ic is False


def test_l0051_has_structured_routes():
    b = lemma_l0051()
    assert b.structured_routes
    ids = {r.id for r in b.structured_routes}
    assert "onepol-cr0008" in ids
    assert "sos-brute-scale" in ids
    assert "greedy-41" in ids
    assert "greedy-sos-scan" in ids


def test_l0051_onepol_sos_validated():
    b = lemma_l0051()
    by_id = {r.id: r for r in b.structured_routes or []}
    r = by_id["sos-onepol-d92"]
    assert r.status == "validated"
    assert r.C_bound_est is not None
    assert r.C_bound_est < b.C_dagger
    assert r.C_bound_est < 8.5
