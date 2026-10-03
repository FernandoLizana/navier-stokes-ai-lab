"""Tests for L-0059 Lean SOS cert bridge."""

from ns_exploration.conjectures.l0059_lean_sos_cert_bridge import CERTS, lemma_l0059


def test_l0059_pins():
    l59 = lemma_l0059(run_lean=False)
    assert l59.lemma_id == "L-0059"
    assert len(l59.pins or []) >= 3
    for p in l59.pins or []:
        assert p.C_ub_hi_num / p.C_ub_hi_den < p.C_dagger_num / p.C_dagger_den


def test_l0059_lean_build():
    l59 = lemma_l0059(run_lean=True)
    assert l59.lean_build_ok
