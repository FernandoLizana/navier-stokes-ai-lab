"""Exploratory temporal scan (N2) for terminal-weighted fullsym bounds."""

from __future__ import annotations

import json
from pathlib import Path

from ns_exploration.conjectures.l0045_streaming_fullsym import streaming_C_fullsym
from ns_exploration.terminal_weighted.constants import FROZEN
from ns_exploration.terminal_weighted.tensor import terminal_fullsym_bound


def _time_grid(T: float = FROZEN.T, step: float = 0.00025) -> list[float]:
    n = int(round(T / step))
    return [round(i * step, 10) for i in range(n + 1)]


def _integrate_trapz(ts: list[float], ys: list[float]) -> float:
    if len(ts) < 2:
        return 0.0
    acc = 0.0
    for i in range(len(ts) - 1):
        h = ts[i + 1] - ts[i]
        acc += 0.5 * h * (ys[i] + ys[i + 1])
    return acc


def _integrate_left(ts: list[float], ys: list[float]) -> float:
    acc = 0.0
    for i in range(len(ts) - 1):
        h = ts[i + 1] - ts[i]
        acc += h * ys[i]
    return acc


def _integrate_right(ts: list[float], ys: list[float]) -> float:
    acc = 0.0
    for i in range(len(ts) - 1):
        h = ts[i + 1] - ts[i]
        acc += h * ys[i + 1]
    return acc


def run_exploratory_scan(
    n: int = 24,
    *,
    step: float = 0.00025,
    band_radii: tuple[int, ...] = tuple(range(1, 9)),
) -> dict:
    """
    N2 exploratory scan on band {1..8} with exact streaming fullsym at each t.

    Full dealias all-IC envelope requires shell-sum or sparse per-t (Phase D).
    """
    ts = _time_grid(FROZEN.T, step)
    grid = []
    bounds = []
    dominant = {"t": 0.0, "C_ub": 0.0, "shell_note": "band only"}

    for t in ts:
        band = terminal_fullsym_bound(n, band_radii, t, method="streaming")
        row = {
            "t": t,
            "C_ub_band_exact": band["C_term_ub"],
            "L_op": band["L_op"],
            "band_radii": list(band_radii),
            "method": "streaming_fullsym_band",
            "evidence_level": "N2",
        }
        grid.append(row)
        bounds.append(band["C_term_ub"])
        if band["C_term_ub"] > dominant["C_ub"]:
            dominant = {
                "t": t,
                "C_ub": band["C_term_ub"],
                "shell_note": f"band {list(band_radii)}",
            }

    max_bound = max(bounds)
    t_max = ts[bounds.index(max_bound)]
    trap = _integrate_trapz(ts, bounds)
    left = _integrate_left(ts, bounds)
    right = _integrate_right(ts, bounds)
    margin = FROZEN.integral_threshold - trap

    # Regression: band at t=T vs L-0045 streaming (unweighted path)
    l45 = streaming_C_fullsym(n, band_radii)
    band_at_T = bounds[-1]

    return {
        "n": n,
        "nu": FROZEN.nu,
        "T": FROZEN.T,
        "step": step,
        "n_times": len(ts),
        "band_radii": list(band_radii),
        "integral_threshold": FROZEN.integral_threshold,
        "uniform_threshold": FROZEN.uniform_threshold,
        "max_grid_bound": max_bound,
        "t_at_max": t_max,
        "dominant": dominant,
        "trapezoid_estimate": trap,
        "left_riemann_estimate": left,
        "right_riemann_estimate": right,
        "margin_vs_threshold": margin,
        "closes_integral_exploratory": trap <= FROZEN.integral_threshold,
        "closes_uniform_exploratory": max_bound <= FROZEN.uniform_threshold,
        "band_at_T": band_at_T,
        "l0045_band_at_T": l45["C_fullsym"],
        "band_T_regression_rel_err": abs(band_at_T - l45["C_fullsym"])
        / max(l45["C_fullsym"], 1e-30),
        "grid": grid,
        "evidence_level": "N2",
        "honesty": (
            "Band {1..8} scan only; full dealias all-IC needs shell-sum or "
            "interval cover (Phase D). Not a C-0008 certificate."
        ),
    }


def save_scan(data: dict, path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2), encoding="utf-8")
