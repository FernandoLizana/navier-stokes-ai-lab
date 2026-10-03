"""Detect incompatible Clay claims in literature CSV (heuristic)."""

from __future__ import annotations

import csv
from pathlib import Path


def find_claim_conflicts(csv_path: str | Path) -> list[str]:
    """
    Flag if both 'global regularity claimed' and 'blow-up claimed' appear
    under ns_mrl_class == claimed_solutions without acceptance.
    """
    path = Path(csv_path)
    if not path.exists():
        return [f"missing:{path}"]
    claims_reg = []
    claims_blow = []
    with path.open(newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            if row.get("ns_mrl_class") != "claimed_solutions":
                continue
            text = (row.get("conclusion", "") + " " + row.get("title", "")).lower()
            if "regular" in text or "smooth" in text:
                claims_reg.append(row.get("id", "?"))
            if "blow" in text or "singul" in text:
                claims_blow.append(row.get("id", "?"))
    out = []
    if claims_reg and claims_blow:
        out.append(
            f"INCOMPATIBLE_CLAIM_CLUSTER: regularity-like {claims_reg} vs blowup-like {claims_blow}"
        )
    return out
