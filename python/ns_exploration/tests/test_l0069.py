"""Tests for L-0069 laptop closure synthesis."""

from ns_exploration.conjectures.l0069_laptop_onepol_closure import lemma_l0069


def test_l0069_closure():
    b = lemma_l0069()
    assert b.lemma_id == "L-0069"
    assert b.n_proved_restricted >= 14
    assert b.best_onepol_cr == "C-R-0014"
    assert b.techo_shell == 25
    assert b.techo_C_shor > b.C_dagger
    assert not b.all_ic_open
    assert b.c0008_sharp_open
