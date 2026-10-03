"""L-0072 terminal route refutation (expected verifier FAIL)."""

from __future__ import annotations

from pathlib import Path

import pytest

from ns_exploration.validation.l0072_certificate import verify_l0072_route_refuted

ROOT = Path(__file__).resolve().parents[3]
CERT = ROOT / "certificates/CERT-L0072-C0008-terminal-full-dealias.json"
MANIFEST = ROOT / "experiments/terminal_weighted/shell_manifest_best_full.json"
REGISTRY = ROOT / "conjectures/active/L-0072.json"


@pytest.mark.skipif(not CERT.is_file(), reason="full-dealias certificate missing")
@pytest.mark.skipif(not MANIFEST.is_file(), reason="Phase E manifest missing")
def test_l0072_route_refuted_verifier_fail():
    refuted, msgs = verify_l0072_route_refuted(CERT, MANIFEST)
    assert refuted, msgs
    assert any("ROUTE_REFUTED" in m for m in msgs)
    assert any("I_term_hi" in m for m in msgs)


@pytest.mark.skipif(not REGISTRY.is_file(), reason="L-0072 registry missing")
def test_l0072_registry_status():
    import json

    reg = json.loads(REGISTRY.read_text(encoding="utf-8"))
    assert reg["status"] == "route_refuted"
    assert reg["closes_c0008"] is False
    assert reg["refutes_conjecture"] is False
    assert reg["verifier"] == "FAIL"
