"""
Finite-dimensional refutable conjecture (NOT a Clay claim).

C-0001 (Route B context, Galerkin truncation only):
  For the semi-discrete pseudospectral NS on T³ with dealiasing, viscosity ν=0.1,
  energy E(u0)=1/2, resolution N≤16, and time T=0.02, the enstrophy satisfies
      E_ens(T) ≤ M
  for an explicit numerical M proposed from sampling — to be REFUTED or
  tightened, never promoted to continuum without proof.

Status workflow: propose → attempt refutation by adversarial FD ascent →
  reformulate if broken.

Evidence starts at N6 (formal statement) with numerical probes ≤ N2.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from ns_exploration.initial_conditions.generators import (
    abc_flow,
    random_div_free,
    taylor_green,
    vortex_tubes_periodic,
)
from ns_exploration.optimization.enstrophy_gradient import (
    evolve_enstrophy,
    one_gradient_ascent_step,
    _project_fixed_energy,
)


CONJECTURE_ID = "C-0001"


@dataclass
class ConjectureC0001:
    id: str = CONJECTURE_ID
    route: str = "B"
    status: str = "proposed"
    evidence_level: str = "N6"
    domain: str = "Galerkin/pseudospectral T^3, N<=16 (NOT continuum PDE)"
    nu: float = 0.1
    energy: float = 0.5
    t_end: float = 0.02
    n_max: int = 16
    proposed_bound_M: float = 20.0
    statement: str = (
        "For all divergence-free Fourier fields on the N^3 grid with N<=16, "
        "dealiased convective term, nu=0.1, kinetic energy = 0.5, evolved with "
        "ETD-RK2 (dt=1e-3) to T=0.02, the enstrophy E_ens(T) <= M "
        "with M=proposed_bound_M. This is a FINITE-DIMENSIONAL conjecture only."
    )
    clay_implication: str = (
        "None. Even if true for all N, passing N->inf requires uniform estimates "
        "not provided here. Does not address Clay Routes A-D."
    )
    numerical_support_max: float | None = None
    refuted: bool = False
    refutation_value: float | None = None
    notes: str = "Prioritize refutation."

    def as_dict(self) -> dict:
        return asdict(self)


def probe_enstrophy_samples(
    n: int = 12,
    nu: float = 0.1,
    dt: float = 1e-3,
    t_end: float = 0.02,
    energy: float = 0.5,
    n_random: int = 30,
) -> dict:
    values = []
    # Structured families
    for gen in (taylor_green, abc_flow, vortex_tubes_periodic):
        uh, _ = gen(n)
        uh = _project_fixed_energy(uh, energy)
        values.append(evolve_enstrophy(uh, nu, dt, t_end))
    for s in range(n_random):
        uh, _ = random_div_free(n, seed=s, energy_target=energy)
        values.append(evolve_enstrophy(uh, nu, dt, t_end))
    # Mild adversarial bumps on a few seeds
    for s in range(5):
        uh, _ = random_div_free(n, seed=100 + s, energy_target=energy)
        uh2, _, j1 = one_gradient_ascent_step(
            uh, nu, dt, t_end, energy, step_size=0.08, seed=200 + s
        )
        values.append(j1)
    arr = np.array(values, dtype=float)
    return {
        "n": n,
        "count": int(arr.size),
        "max": float(arr.max()),
        "mean": float(arr.mean()),
        "min": float(arr.min()),
        "values_hash_max": float(arr.max()),
    }


def evaluate_conjecture(
    M: float = 20.0,
    n: int = 12,
) -> ConjectureC0001:
    probe = probe_enstrophy_samples(n=n)
    c = ConjectureC0001(proposed_bound_M=M, numerical_support_max=probe["max"])
    if probe["max"] > M:
        c.status = "refuted"
        c.refuted = True
        c.refutation_value = probe["max"]
        c.evidence_level = "N2"  # numerical refutation evidence
    else:
        c.status = "exploring"
        c.notes = (
            f"No refutation in {probe['count']} probes; max={probe['max']:.4f} < M={M}. "
            "Still unproved; continue adversarial search."
        )
    return c


def save_conjecture(c: ConjectureC0001, path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(c.as_dict(), indent=2), encoding="utf-8")
