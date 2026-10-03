"""Tests for L-0056 closure map."""

from __future__ import annotations

from ns_exploration.conjectures.l0056_closure_map import (
    lemma_l0056,
    render_closure_map_md,
    write_closure_map_report,
)


def test_l0056_lemma():
    l56 = lemma_l0056()
    assert l56.lemma_id == "L-0056"
    assert l56.n_proved_subclasses >= 2
    assert l56.n_open_subclasses >= 1
    assert l56.best_n5_sos_hi is not None
    assert l56.best_n5_sos_hi < l56.C_dagger
    assert l56.best_structured_sos_hi is not None
    assert l56.best_structured_sos_hi > l56.best_n5_sos_hi


def test_l0056_report(tmp_path):
    l56 = lemma_l0056()
    md = render_closure_map_md(l56)
    assert "C-0007" in md
    assert "C-R-0012" in md
    assert "L-0055" in md
    p = write_closure_map_report(l56, tmp_path / "closure_map.md")
    assert p.exists()
