"""Tests for L-0063 greedy export status."""

from ns_exploration.conjectures.l0063_greedy_cluster_export import lemma_l0063


def test_l0063_status():
    l63 = lemma_l0063()
    assert l63.lemma_id == "L-0063"
    assert l63.n_shells == 41
    assert l63.greedy_D == 1772
    assert l63.status in ("cluster_export_ready", "cluster_export_pending")
