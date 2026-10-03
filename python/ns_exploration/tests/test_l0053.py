"""Tests for L-0053 one-pol N5 Gram."""

import pytest

from ns_exploration.conjectures.l0053_onepol_sos_gram_n5 import GRAM_JSON, lemma_l0053


@pytest.mark.skipif(not GRAM_JSON.is_file(), reason="run export_tssos_gram.jl on onepol_cr0008 first")
def test_l0053_gram_psd():
    b = lemma_l0053()
    assert b.closes_vs_Cdagger
    assert b.beats_fullsym
    assert b.n_gram_blocks > 0
    assert b.min_psd_eig >= -5e-8


@pytest.mark.skipif(not GRAM_JSON.is_file(), reason="run export_tssos_gram.jl on onepol_cr0008 first")
def test_l0053_certificate_verify():
    from ns_exploration.validation.l0053_certificate import build_l0053_certificate, verify_l0053_certificate

    cert = build_l0053_certificate()
    ok, checks = verify_l0053_certificate(cert.as_dict())
    assert ok, checks
