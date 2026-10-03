"""Tests for L-0065 one-pol greedy-second SOS."""

from pathlib import Path

import pytest

from ns_exploration.conjectures.l0065_onepol_greedy_second_sos import lemma_l0065
from ns_exploration.validation.l0065_certificate import verify_l0065_certificate

TSSOS = Path("tools/sos_julia/data/onepol_greedy_second/cosmo_tssos_result.json")


@pytest.mark.skipif(not TSSOS.is_file(), reason="run COSMO TSSOS on onepol_greedy_second first")
def test_l0065_interval():
    b = lemma_l0065()
    assert b.lemma_id == "L-0065"
    assert b.D == 154
    assert b.closes_vs_Cdagger
    assert b.beats_shor
    assert b.C_ub_hi < b.C_dagger - 1e-9


@pytest.mark.skipif(not TSSOS.is_file(), reason="run COSMO TSSOS on onepol_greedy_second first")
def test_l0065_certificate():
    from ns_exploration.validation.l0065_certificate import build_l0065_certificate

    cert = build_l0065_certificate()
    ok, checks = verify_l0065_certificate(cert.as_dict())
    assert ok, checks
