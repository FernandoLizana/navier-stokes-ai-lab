"""Tests for L-0066 / C-R-0013."""

from pathlib import Path

from ns_exploration.conjectures.l0043_onepol_shellblock import SUPPORT_CR0008
from ns_exploration.conjectures.l0065_onepol_greedy_second_sos import SUPPORT_GREEDY_SECOND
from ns_exploration.conjectures.l0066_cr0013_onepol_greedy import lemma_l0066


def test_l0066_closes():
    b = lemma_l0066()
    assert b.closes_cr0013
    assert b.C_shor_sym_cr0013 <= b.C_dagger + 1e-9
    assert len(b.support_cr0013 or []) == len(SUPPORT_GREEDY_SECOND)


def test_l0066_extends_cr0008():
    b = lemma_l0066()
    assert set(SUPPORT_CR0008).issubset(set(b.support_cr0013 or []))
    assert b.extends_cr0008


def test_cr0013_json_exists():
    p = Path("conjectures/proved_restricted/C-R-0013.json")
    if p.is_file():
        import json

        d = json.loads(p.read_text(encoding="utf-8"))
        assert d["id"] == "C-R-0013"
        assert d["status"] == "proved_restricted"
