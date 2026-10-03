"""Tests for L-0058 band {1..5} N5 Gram."""

from pathlib import Path

import pytest

from ns_exploration.conjectures.l0058_band12345_sos_gram_n5 import GRAM_JSON, lemma_l0058
from ns_exploration.validation.l0058_certificate import build_l0058_certificate, verify_l0058_certificate


@pytest.mark.skipif(not GRAM_JSON.is_file(), reason="run export_tssos_gram.jl on band_12345 first")
def test_l0058_gram_psd():
    b = lemma_l0058()
    assert b.closes_vs_Cdagger
    assert b.beats_fullsym
    assert b.n_gram_blocks > 0


@pytest.mark.skipif(not GRAM_JSON.is_file(), reason="run export_tssos_gram.jl on band_12345 first")
def test_l0058_certificate_verify():
    cert = build_l0058_certificate()
    ok, checks = verify_l0058_certificate(cert.as_dict())
    assert ok, checks
