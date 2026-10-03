"""L-0076 N<=24 inclusion audit tests."""

from __future__ import annotations

from gmpy2 import mpfr

from ns_exploration.terminal_weighted.n_coverage import (
    full_dealias_inclusion_audit,
    load_repair_manifest,
    n_domination_analysis,
    repair_manifest_integral_hi,
    wavevector_shells_subset,
)


def test_wavevector_subsets_embed_in_n24():
    for n in (12, 16, 20):
        rep = wavevector_shells_subset(n, 24)
        assert rep["subset"] is True


def test_embed_basis_wavevectors_in_n24():
    from ns_exploration.terminal_weighted.n_coverage import embed_basis_index_check

    for n in (12, 16, 20):
        rep = embed_basis_index_check(n, 24)
        assert rep["all_wavevectors_embedded"] is True


def test_n24_repair_manifest_integral_matches_cert_order():
    manifest = load_repair_manifest(24)
    if manifest is None:
        return
    I_hi = repair_manifest_integral_hi(manifest, prec=128)
    assert float(I_hi) > 100.0
    assert len(manifest["blocks"]) == 87


def test_full_dealias_inclusion_audit_structure():
    rep = full_dealias_inclusion_audit(prec=128)
    assert rep["lemma"] == "L-0076"
    assert rep["wavevector_domination"]["all_subset"] is True
    r24 = rep.get("repair_full_dealias_n24", {})
    if r24:
        assert r24["repair_complete"] is True
        assert mpfr(r24["I_hi_route_A"]) > mpfr("170")
