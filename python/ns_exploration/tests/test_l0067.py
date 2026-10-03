"""Tests for L-0067 / L-0068 / C-R-0014."""

from pathlib import Path

import pytest

from ns_exploration.conjectures.l0065_onepol_greedy_second_sos import SUPPORT_GREEDY_SECOND
from ns_exploration.conjectures.l0067_onepol_greedy_r24_sos import SUPPORT_CR0014, lemma_l0067
from ns_exploration.conjectures.l0068_cr0014_onepol_greedy import lemma_l0068

COSMO = Path("tools/sos_julia/data/onepol_greedy_second_r24/cosmo_tssos_result.json")


@pytest.mark.skipif(not COSMO.is_file(), reason="COSMO on onepol_greedy_second_r24 required")
def test_l0067_closes():
    b = lemma_l0067()
    assert b.closes_vs_Cdagger
    assert b.D == 178
    assert b.support_radii == SUPPORT_CR0014


@pytest.mark.skipif(not COSMO.is_file(), reason="COSMO required")
def test_l0068_extends_cr0013():
    b = lemma_l0068()
    assert b.closes_cr0014
    assert set(SUPPORT_GREEDY_SECOND).issubset(set(b.support_cr0014 or []))
    assert 24 in (b.support_cr0014 or [])


def test_cr0014_file_after_sprint():
    p = Path("conjectures/proved_restricted/C-R-0014.json")
    if p.is_file():
        import json

        assert json.loads(p.read_text())["id"] == "C-R-0014"
