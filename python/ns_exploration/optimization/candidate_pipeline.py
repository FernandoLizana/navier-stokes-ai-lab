"""
Multi-seed resolution gate battery + single candidate extraction.

Evidence ≤ N2. Route B. Not a singularity search at scale.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from ns_exploration.diagnostics.metrics import compute_diagnostics
from ns_exploration.initial_conditions.generators import random_div_free, save_ic
from ns_exploration.optimization.adjoint import verify_adjoint_vs_fd
from ns_exploration.optimization.enstrophy_gradient import (
    _project_fixed_energy,
    _resample_hat,
    evolve_enstrophy,
    one_gradient_ascent_step,
    resolution_doubling_gate,
)
from ns_exploration.spectral.fourier_conventions import kinetic_energy_from_hat
from ns_exploration.spectral.integrators import step_etd_rk2, step_rk4
from ns_exploration.spectral.leray import divergence_l2, leray_project_hat


@dataclass
class CandidateRecord:
    seed: int
    n: int
    n_verify: int
    J0: float
    J_after_ascent: float
    J_etd_n: float
    J_rk4_n: float
    J_etd_2n: float
    integrator_rel_diff_n: float
    resolution_rel_diff_etd: float
    div_l2: float
    energy: float
    accepted: bool
    reject_reason: str
    evidence_level: str = "N2"
    route: str = "B"

    def as_dict(self) -> dict:
        return asdict(self)


def run_gate_battery(
    seeds: range | list[int] | None = None,
    n_low: int = 12,
    gate_tol: float = 0.35,
    nu: float = 0.1,
    dt: float = 1e-3,
    t_end: float = 0.02,
    energy_target: float = 0.5,
) -> list[dict]:
    if seeds is None:
        seeds = range(20)
    rows = []
    for s in seeds:
        g = resolution_doubling_gate(
            n_low=n_low,
            seed=s,
            nu=nu,
            dt=dt,
            t_end=t_end,
            energy_target=energy_target,
            gate_tol=gate_tol,
        )
        rows.append(g.as_dict())
    return rows


def evolve_enstrophy_integrator(
    u_hat: np.ndarray,
    nu: float,
    dt: float,
    t_end: float,
    integrator: str,
) -> float:
    u = u_hat.copy()
    step = step_etd_rk2 if integrator == "etd_rk2" else step_rk4
    nsteps = int(round(t_end / dt))
    for _ in range(nsteps):
        u = step(u, dt, nu, dealias=True)
    return compute_diagnostics(u, nu).enstrophy


def extract_one_candidate(
    gate_rows: list[dict],
    nu: float = 0.1,
    dt: float = 1e-3,
    t_end: float = 0.02,
    energy_target: float = 0.5,
    integrator_tol: float = 0.05,
    resolution_tol: float = 0.08,
    out_dir: str | Path = "datasets/candidates",
) -> CandidateRecord | None:
    """
    Pick the first gate-passing seed, do one ascent step, cross-check RK4 vs ETD
    and N vs 2N. Accept only if both gates pass.
    """
    passed = [r for r in gate_rows if r["gate_passed"]]
    if not passed:
        return None
    # Prefer smallest relative agreement
    passed.sort(key=lambda r: r["relative_agreement"])
    seed = int(passed[0]["seed"])
    n = int(passed[0]["n_low"])

    u0, meta = random_div_free(n, seed=seed, energy_target=energy_target)
    u1, j0, j1 = one_gradient_ascent_step(
        u0, nu=nu, dt=dt, t_end=t_end, energy_target=energy_target, seed=seed + 99
    )
    u1 = _project_fixed_energy(u1, energy_target)

    j_etd = evolve_enstrophy_integrator(u1, nu, dt, t_end, "etd_rk2")
    j_rk4 = evolve_enstrophy_integrator(u1, nu, dt, t_end, "rk4")
    integ_rel = abs(j_etd - j_rk4) / (abs(j_etd) + abs(j_rk4) + 1e-30)

    u2n = _project_fixed_energy(_resample_hat(u1, 2 * n), energy_target)
    j_etd_2n = evolve_enstrophy_integrator(u2n, nu, dt, t_end, "etd_rk2")
    res_rel = abs(j_etd - j_etd_2n) / (abs(j_etd) + abs(j_etd_2n) + 1e-30)

    reason = ""
    accepted = True
    if integ_rel > integrator_tol:
        accepted = False
        reason = f"integrator_disagreement:{integ_rel}"
    elif res_rel > resolution_tol:
        accepted = False
        reason = f"resolution_disagreement:{res_rel}"

    rec = CandidateRecord(
        seed=seed,
        n=n,
        n_verify=2 * n,
        J0=j0,
        J_after_ascent=j1,
        J_etd_n=j_etd,
        J_rk4_n=j_rk4,
        J_etd_2n=j_etd_2n,
        integrator_rel_diff_n=integ_rel,
        resolution_rel_diff_etd=res_rel,
        div_l2=divergence_l2(u1),
        energy=kinetic_energy_from_hat(u1),
        accepted=accepted,
        reject_reason=reason,
    )

    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    save_ic(str(out / f"candidate_seed{seed}_n{n}.npz"), u1, meta)
    (out / f"candidate_seed{seed}_n{n}.json").write_text(
        json.dumps(rec.as_dict(), indent=2), encoding="utf-8"
    )
    return rec
