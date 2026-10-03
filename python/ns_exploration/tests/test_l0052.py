"""Tests for L-0052 N5 Gram rationalization."""

from pathlib import Path

import pytest

from ns_exploration.conjectures.l0052_sos_gram_n5 import GRAM_JSON, lemma_l0052, verify_gram_blocks
from ns_exploration.validation.l0052_certificate import verify_l0052_certificate


@pytest.mark.skipif(not GRAM_JSON.is_file(), reason="run export_tssos_gram.jl first")
def test_l0052_gram_psd():
    b = lemma_l0052()
    assert b.closes_vs_Cdagger
    assert b.beats_fullsym
    assert b.n_gram_blocks > 0
    assert b.min_psd_eig >= -5e-8


@pytest.mark.skipif(not GRAM_JSON.is_file(), reason="run export_tssos_gram.jl first")
def test_l0052_certificate_verify():
    from ns_exploration.validation.l0052_certificate import build_l0052_certificate

    cert = build_l0052_certificate()
    ok, checks = verify_l0052_certificate(cert.as_dict())
    assert ok, checks
