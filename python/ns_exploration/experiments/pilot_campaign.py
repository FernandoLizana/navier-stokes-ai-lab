"""Pilot campaign: infrastructure validation only. NOT singularity search.

Budget caps (Sprint 1):
- ≤ 100 ICs
- 3 resolutions
- 3 integrators
- 5 seeds per random family
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

from ns_exploration.initial_conditions.generators import (
    abc_flow,
    random_div_free,
    taylor_green,
    vortex_tubes_periodic,
)
from ns_exploration.spectral.fourier_conventions import kinetic_energy_from_hat
from ns_exploration.spectral.leray import divergence_l2
from ns_exploration.spectral.solver import NavierStokesSolver, SolverConfig


def run_pilot(
    out_dir: str | Path = "experiments/exploratory/pilot_sprint01",
    max_ics: int = 100,
) -> Path:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    resolutions = [16, 24, 32]
    integrators = ["rk4", "etd_rk2", "semi_implicit"]
    seeds = [0, 1, 2, 3, 4]
    nu = 0.05
    dt = 5e-4
    t_end = 0.02

    rows = []
    count = 0

    def add_run(name, u_hat, n, integ, seed=None):
        nonlocal count
        if count >= max_ics:
            return
        # resample / regenerate at resolution n if needed
        cfg = SolverConfig(n=n, nu=nu, dt=dt, t_end=t_end, integrator=integ, dealias=True)
        if u_hat.shape[-1] != n:
            raise ValueError("IC resolution mismatch")
        solver = NavierStokesSolver(cfg)
        state = solver.run(u_hat)
        final = state.history[-1]
        rows.append(
            {
                "name": name,
                "n": n,
                "integrator": integ,
                "seed": seed if seed is not None else "",
                "t_final": state.t,
                "steps": state.step,
                "energy0": state.history[0]["energy"],
                "energy_final": final["energy"],
                "enstrophy_final": final["enstrophy"],
                "div_final": final["div_l2"],
                "stopped": state.stopped_reason or "",
                "route": "B",
                "evidence": "N2",
            }
        )
        count += 1

    # Deterministic families × resolutions × integrators
    for n in resolutions:
        for integ in integrators:
            uh, _ = taylor_green(n)
            add_run("taylor_green", uh, n, integ)
            uh, _ = abc_flow(n)
            add_run("abc", uh, n, integ)
            uh, _ = vortex_tubes_periodic(n)
            add_run("vortex_tubes", uh, n, integ)
            for s in seeds:
                uh, _ = random_div_free(n, seed=s, energy_target=0.5)
                add_run("random_div_free", uh, n, integ, seed=s)

    csv_path = out / "pilot_results.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    summary = {
        "n_runs": len(rows),
        "max_div": max(r["div_final"] for r in rows),
        "n_stopped": sum(1 for r in rows if r["stopped"]),
        "route": "B",
        "purpose": "infrastructure_validation",
        "evidence": "N2",
        "not_a_singularity_search": True,
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return csv_path


if __name__ == "__main__":
    run_pilot()
