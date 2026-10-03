"""Tests for L-0048 sparse fullsym techo (fast: uses frozen witness + small validate)."""

from __future__ import annotations

from ns_exploration.conjectures.l0048_fullsym_techo import lemma_l0048, save_lemma_l0048
from ns_exploration.conjectures.l0048_matrixfree_fullsym import validate_sparse_vs_dense


def test_l0048_validate_and_techo_frozen():
    val = validate_sparse_vs_dense()
    assert val["ok"]
    b = lemma_l0048(n=24, recompute=False)
    assert b.fullsym_shor_all_ic_techo
    assert not b.closes_c0007_all_ic
    assert b.C_fullsym_full > b.C_dagger
    assert b.D_full == 6748
    save_lemma_l0048(b, path="conjectures/proved_restricted/L-0048.json")
