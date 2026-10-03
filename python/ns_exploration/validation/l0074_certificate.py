"""L-0074 Route A ladder partition refutation — integrity checks."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
LADDER = ROOT / "experiments/terminal_weighted/route_a_ladder_full.json"
REGISTRY = ROOT / "conjectures/active/L-0074.json"


def verify_l0074_ladder(*, ladder_path: Path | None = None) -> tuple[bool, list[str]]:
    ladder_path = ladder_path or LADDER
    msgs: list[str] = []
    if not ladder_path.is_file():
        return False, [f"missing {ladder_path}"]

    ladder = json.loads(ladder_path.read_text(encoding="utf-8"))
    if not ladder.get("complete"):
        msgs.append("FAIL: ladder not complete")
        return False, msgs

    ok_rows = [r for r in ladder["rows"] if "error" not in r]
    if len(ok_rows) != 4:
        msgs.append(f"FAIL: expected 4 rows, got {len(ok_rows)}")
        return False, msgs

    best = ladder["best"]
    if best["label"] != "equal_12":
        msgs.append(f"FAIL: unexpected best label {best['label']}")
        return False, msgs

    I_best = float(best["best_I_hi"])
    for r in ok_rows:
        if r["label"] == "equal_12":
            continue
        if float(r["best_I_hi"]) < I_best:
            msgs.append(f"FAIL: {r['label']} beats equal_12")
            return False, msgs

    if float(best["omega_T_hi"]) >= float(ladder["M_target"]):
        msgs.append("omega_T_hi >= M_target (expected FAIL for C-0008)")

    msgs.append("L-0074 partition ladder consistent: equal_12 optimal in family")
    return True, msgs
