"""Tests for L-0060 band {1..6} N5."""

import pytest

from ns_exploration.conjectures.l0060_band123456_sos_gram_n5 import (
    GRAM_JSON,
    TSSOS_JSON,
    lemma_l0060,
)
from ns_exploration.validation.l0060_certificate import build_l0060_certificate, verify_l0060_certificate


@pytest.mark.skipif(not TSSOS_JSON.is_file(), reason="missing band_123456 tssos_result")
def test_l0060_interval_or_gram():
    b = lemma_l0060(require_gram=GRAM_JSON.is_file())
    assert b.closes_vs_Cdagger
    assert b.beats_fullsym
    assert b.status in ("validated", "validated_constant")


@pytest.mark.skipif(not TSSOS_JSON.is_file(), reason="missing band_123456 tssos_result")
def test_l0060_certificate_verify():
    cert = build_l0060_certificate()
    ok, checks = verify_l0060_certificate(cert.as_dict())
    assert ok, checks
