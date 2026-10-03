"""Tests for L-0029 Frobenius stretch embedding."""

from __future__ import annotations

import math

from ns_exploration.conjectures.l0029_frobenius_stretch import (
    lemma_l0029,
    stretch_constant_frobenius,
)
from ns_exploration.validation.l0029_certificate import (
    build_l0029_certificate,
    verify_l0029_certificate,
)


def test_l0029_frobenius_improves_but_techo():
    b = lemma_l0029(n=24)
    assert abs(b.C_L0002_full / b.C_F_full - math.sqrt(3.0)) < 1e-9
    assert b.C_F_full > b.C_dagger
    assert b.M_star_F < 12
    assert b.shell_r1_meets_Cdagger
    assert not b.shell_r1_reaches_low_slab_edge
    assert not b.frobenius_closes_c0007
    assert abs(stretch_constant_frobenius(6) - 2 * math.sqrt(12)) < 1e-12


def test_l0029_certificate():
    cert = build_l0029_certificate(n=24)
    ok, checks = verify_l0029_certificate(cert.as_dict())
    assert ok and all(checks.values())
