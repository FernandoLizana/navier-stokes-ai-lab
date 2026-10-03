"""Per-shell majorants and shell-sum temporal upper bounds (spec §13)."""

from __future__ import annotations

import json
import math
from pathlib import Path

from ns_exploration.conjectures.l0048_matrixfree_fullsym import all_dealias_radii
from ns_exploration.terminal_weighted.constants import FROZEN
from ns_exploration.terminal_weighted.tensor import terminal_fullsym_bound
from ns_exploration.terminal_weighted.weights import terminal_weight


def compute_shell_bounds(
    n: int = 24,
    t: float = FROZEN.T,
    *,
    radii: tuple[int, ...] | None = None,
) -> dict:
    """
    For each shell r, B_r = 2√2 ‖T_r‖ with output restricted to shell r at t=T
    (weight w_r(T)=r on that shell). Sum bound: C_ub(t) <= sum_r e^{-2νr(T-t)} B_r.
    """
    radii = radii or all_dealias_radii(n)
    shell_radii = sorted(set(radii))
    rows = []
    for r in shell_radii:
        row = terminal_fullsym_bound(
            n, tuple(radii), t, output_shell=r, method="streaming"
        )
        row["shell"] = r
        row["B_r"] = row["C_term_ub"]
        rows.append(row)
    return {
        "n": n,
        "t_reference": t,
        "shells": rows,
        "sum_B_r": sum(r["B_r"] for r in rows),
        "evidence_level": "N2",
        "note": "sum B_r >= C_term(T) by triangle inequality on shell blocks",
    }


def shell_sum_bound_at_t(
    shell_bounds: list[dict],
    t: float,
    nu: float = FROZEN.nu,
    T: float = FROZEN.T,
) -> float:
    acc = 0.0
    for row in shell_bounds:
        r = float(row["shell"])
        B_r = float(row["B_r"])
        acc += math.exp(-2.0 * nu * r * (T - t)) * B_r
    return acc


def save_shell_bounds(data: dict, path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2), encoding="utf-8")
