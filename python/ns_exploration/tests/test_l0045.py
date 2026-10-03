"""Tests for L-0045 streaming physical fullsym and C-R-0010.

Slow (~15–25 min): recomputes streamed SVD on D≈1028 and D≈1172.
"""

from __future__ import annotations

from ns_exploration.conjectures.l0045_physical_band import (
    SUPPORT_CR0010,
    cr0010_record,
    lemma_l0045,
)
from ns_exploration.validation.l0045_certificate import (
    build_l0045_certificate,
    verify_l0045_certificate,
)


def test_l0045_cr0010_and_certificate():
    b = lemma_l0045(n=24)
    assert b.support_cr0010 == list(SUPPORT_CR0010)
    assert b.C_fullsym_rmax25 <= b.C_dagger
    assert b.C_fullsym_rmax26 > b.C_dagger
    assert b.closes_cr0010
    assert not b.closes_c0007_all_ic
    cr = cr0010_record(b)
    assert cr["id"] == "C-R-0010"
    assert cr["C_fullsym"] <= cr["C_dagger"]
    cert = build_l0045_certificate(n=24, bound=b)
    ok, checks = verify_l0045_certificate(cert.as_dict())
    assert ok and all(checks.values())
