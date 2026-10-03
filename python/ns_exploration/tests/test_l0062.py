"""Tests for L-0062 COSMO ladder."""

from ns_exploration.conjectures.l0062_cosmo_ladder_n5 import COSMO_BANDS, lemma_l0062
from ns_exploration.validation.l0062_certificate import (
    build_l0062_certificates,
    verify_l0062_certificate,
)


def test_l0062_three_bands():
    l62 = lemma_l0062()
    assert len(l62.rows or []) == len(COSMO_BANDS)
    for r in l62.rows or []:
        assert r.C_ub_hi < r.C_dagger
        assert r.C_ub_hi < r.C_fullsym


def test_l0062_certs():
    certs = build_l0062_certificates()
    assert len(certs) == 3
    for c in certs:
        ok, checks = verify_l0062_certificate(c.as_dict())
        assert ok, checks
