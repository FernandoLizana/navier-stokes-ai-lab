"""Minimal 3D periodic pseudospectral Navier–Stokes solver (exploratory ≤ N2)."""

from __future__ import annotations

import json
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path

import numpy as np

from ns_exploration.diagnostics.metrics import Diagnostics, compute_diagnostics
from ns_exploration.spectral.dealias import dealias_mask
from ns_exploration.spectral.fourier_conventions import enforce_reality, kinetic_energy_from_hat
from ns_exploration.spectral.integrators import INTEGRATORS
from ns_exploration.spectral.leray import divergence_l2, leray_project_hat


@dataclass
class SolverConfig:
    n: int = 32
    nu: float = 0.01
    dt: float = 0.001
    t_end: float = 0.1
    integrator: str = "rk4"
    dealias: bool = True
    force_hat: np.ndarray | None = None
    checkpoint_every: int = 0
    checkpoint_dir: str = "datasets/trajectories"
    # Safety stops (exploratory heuristics — not proofs)
    max_div: float = 1e-6
    max_energy_drift_rel: float = 0.5  # relative; with f=0 expect decay
    min_spectral_tail_ratio: float = 1e-16  # stop if unused
    cfl_limit: float = 0.5
    route: str = "B"
    evidence_level: str = "N2"


@dataclass
class SolverState:
    t: float
    u_hat: np.ndarray
    step: int
    history: list[dict] = field(default_factory=list)
    stopped_reason: str | None = None


class NavierStokesSolver:
    """Pseudospectral NS on T³. Does NOT prove regularity or blow-up."""

    def __init__(self, config: SolverConfig):
        if config.integrator not in INTEGRATORS:
            raise ValueError(f"Unknown integrator {config.integrator}")
        self.config = config
        self.mask = dealias_mask(config.n) if config.dealias else None
        self._step_fn = INTEGRATORS[config.integrator]

    def run(self, u_hat0: np.ndarray) -> SolverState:
        cfg = self.config
        u_hat = leray_project_hat(enforce_reality(u_hat0.copy()))
        if u_hat.shape != (3, cfg.n, cfg.n, cfg.n):
            raise ValueError(f"Expected shape (3,{cfg.n},{cfg.n},{cfg.n}), got {u_hat.shape}")

        state = SolverState(t=0.0, u_hat=u_hat, step=0)
        e0 = kinetic_energy_from_hat(u_hat)
        nsteps = int(np.ceil(cfg.t_end / cfg.dt))

        for i in range(nsteps):
            diag = compute_diagnostics(u_hat, cfg.nu)
            state.history.append({"t": state.t, **diag.as_dict()})

            reason = self._check_stop(u_hat, diag, e0)
            if reason:
                state.stopped_reason = reason
                break

            u_hat = self._step_fn(
                u_hat,
                cfg.dt,
                cfg.nu,
                f_hat=cfg.force_hat,
                dealias=cfg.dealias,
            )
            state.t += cfg.dt
            state.step = i + 1
            state.u_hat = u_hat

            if cfg.checkpoint_every and (i + 1) % cfg.checkpoint_every == 0:
                self._checkpoint(state)

        # final diagnostics
        diag = compute_diagnostics(u_hat, cfg.nu)
        state.history.append({"t": state.t, **diag.as_dict()})
        state.u_hat = u_hat
        return state

    def _check_stop(self, u_hat: np.ndarray, diag: Diagnostics, e0: float) -> str | None:
        cfg = self.config
        div = divergence_l2(u_hat)
        if div > cfg.max_div:
            return f"divergence_violation:{div}"
        if e0 > 0 and diag.energy > e0 * (1.0 + cfg.max_energy_drift_rel) and cfg.force_hat is None:
            return f"energy_growth_unexpected:{diag.energy}"
        # crude CFL using max |u|
        if diag.u_inf * cfg.dt * (cfg.n / (2 * np.pi)) > cfg.cfl_limit:
            return f"cfl_exceeded:u_inf={diag.u_inf}"
        return None

    def _checkpoint(self, state: SolverState) -> None:
        path = Path(self.config.checkpoint_dir)
        path.mkdir(parents=True, exist_ok=True)
        fname = path / f"ckpt_n{self.config.n}_{self.config.integrator}_t{state.t:.6f}.npz"
        np.savez_compressed(
            fname,
            t=state.t,
            u_hat_real=state.u_hat.real,
            u_hat_imag=state.u_hat.imag,
            config=json.dumps({k: v for k, v in asdict(self.config).items() if k != "force_hat"}),
        )
