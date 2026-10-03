"""Tests for L-0064 one-pol scanner."""

import pytest

from ns_exploration.conjectures.l0043_onepol_shellblock import SUPPORT_CR0008
from ns_exploration.conjectures.l0064_onepol_scanner import lemma_l0064


@pytest.fixture(scope="module")
def scan():
    return lemma_l0064(
        n=24,
        max_consecutive=12,
        max_extend_shell=15,
        max_greedy_shell=15,
    )


def test_l0064_scan_runs(scan):
    assert scan.lemma_id == "L-0064"
    assert scan.n_candidates > 0
    assert scan.n_closing >= 1
    assert scan.rows is not None


def test_l0064_cr0008_closes(scan):
    cr_rows = [
        r
        for r in scan.rows or []
        if tuple(r.radii) == tuple(SUPPORT_CR0008) and r.branch == "first"
    ]
    assert cr_rows
    assert cr_rows[0].closes_C_dagger


def test_l0064_consecutive_9_fails(scan):
    r9 = next(r for r in scan.rows or [] if r.scan_id == "consec-1to9-first")
    assert not r9.closes_C_dagger


def test_l0064_greedy_at_least_cr0008(scan):
    assert scan.greedy_growth_first is not None
    assert len(scan.greedy_growth_first) >= len(SUPPORT_CR0008)
