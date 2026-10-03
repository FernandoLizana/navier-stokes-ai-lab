"""Tests for L-0054 structured SOS synthesis."""

from ns_exploration.conjectures.l0054_structured_sos_synthesis import lemma_l0054


def test_l0054_has_n5_rows():
    b = lemma_l0054()
    assert b.all_ic_open
    assert b.best_sos_n5_hi > 0
    assert b.best_sos_n5_hi < b.C_dagger
    ids = {r.id for r in b.rows or []}
    assert "band-123-n5" in ids
    assert "onepol-cr0008-n5" in ids
    assert "full-dealias-shor" in ids


def test_l0054_best_is_band123():
    b = lemma_l0054()
    assert "1,2,3" in b.best_sos_subclass
