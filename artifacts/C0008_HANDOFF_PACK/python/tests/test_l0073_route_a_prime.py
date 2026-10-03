"""L-0073 Route A' cluster certificate (expected closure FAIL, integrity checks)."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[3]
CERT = ROOT / "certificates/CERT-L0073-C0008-route-a-prime-best.json"
MANIFEST = ROOT / "experiments/terminal_weighted/shell_manifest_cluster_best_full.json"


@pytest.mark.skipif(not CERT.is_file(), reason="L-0073 certificate missing")
@pytest.mark.skipif(not MANIFEST.is_file(), reason="cluster best manifest missing")
def test_l0073_certificate_integrity_and_expected_fail():
    import sys

    sys.path.insert(0, str(ROOT))
    from verify_c0008_cluster_certificate import verify_certificate

    ok, msgs = verify_certificate(CERT, manifest_path=MANIFEST)
    assert not ok
    assert any("I_term_hi" in m for m in msgs)
    assert any("omega_T_hi" in m for m in msgs)
    assert not any("hash mismatch" in m for m in msgs)


@pytest.mark.skipif(not (ROOT / "conjectures/active/L-0073.json").is_file(), reason="registry missing")
def test_l0073_registry():
    reg = json.loads((ROOT / "conjectures/active/L-0073.json").read_text(encoding="utf-8"))
    assert reg["status"] == "best_known_bound"
    assert reg["closes_c0008"] is False
    assert reg["I_term_hi"] < reg["phase_e_route_d_I_hi"]
