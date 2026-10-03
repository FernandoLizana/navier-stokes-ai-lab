"""Tests for L-0061 band {1..8} N5."""

import pytest

from ns_exploration.conjectures.l0061_band12345678_sos_n5 import TSSOS_JSON, lemma_l0061
from ns_exploration.validation.l0061_certificate import build_l0061_certificate, verify_l0061_certificate


@pytest.mark.skipif(not TSSOS_JSON.is_file(), reason="missing cosmo_tssos_result")
def test_l0061_interval():
    b = lemma_l0061(require_gram=False)
    assert b.closes_vs_Cdagger
    assert b.beats_fullsym
    assert b.status == "validated_constant"


@pytest.mark.skipif(not TSSOS_JSON.is_file(), reason="missing cosmo_tssos_result")
def test_l0061_certificate():
    cert = build_l0061_certificate()
    ok, checks = verify_l0061_certificate(cert.as_dict())
    assert ok, checks
