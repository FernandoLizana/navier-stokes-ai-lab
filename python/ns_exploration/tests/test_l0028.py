"""Tests for L-0028 embedding techo."""

from __future__ import annotations

from ns_exploration.conjectures.l0028_embedding_techo import lemma_l0028
from ns_exploration.validation.l0028_certificate import (
    build_l0028_certificate,
    verify_l0028_certificate,
)


def test_l0028_embedding_cannot_meet_cdagger():
    b = lemma_l0028(n=24, empirical=False)
    assert b.M_star < 6
    assert b.M_shell_r1 == 6
    assert b.C_shell_r1 > b.C_dagger
    assert b.C_hard_R2 > b.C_dagger
    assert not b.embedding_can_meet_Cdagger
    assert not b.proved_closes_c0007


def test_l0028_certificate():
    cert = build_l0028_certificate(n=24)
    ok, checks = verify_l0028_certificate(cert.as_dict())
    assert ok and all(checks.values())
