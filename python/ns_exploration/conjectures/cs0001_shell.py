"""
C-S-0001 (C-0002-shell): enstrophy bound for low-mode initial data only.

Hypothesis (finite-dimensional, NOT Clay):
  - u0 divergence-free on N^3 grid, N ≤ 16
  - Fourier support of u0 satisfies |k|_∞ ≤ k_shell (default 4)
  - kinetic energy E(u0) = 0.5
  - evolved with dealiased ETD-RK2, ν=0.1, dt=1e-3, to T=0.02
    on the FULL dealias grid (nonlinear term may excite |k|>k_shell)

Conclusion:
  Ω(T) ≤ M_shell

This is STRICTLY weaker than C-0002 (all ICs). High modes can appear during
evolution, so L-0004 shell a priori bounds do NOT close the proof by themselves.

Companion proved statement: L-0005 (Galerkin truncation stays in shell).
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path


@dataclass
class ConjectureCS0001:
    id: str = "C-S-0001"
    route: str = "B"
    status: str = "proposed"
    evidence_level: str = "N6"
    domain: str = (
        "Pseudospectral T^3, N<=16, IC with |k|_inf <= k_shell only "
        "(evolution may leave the shell). NOT continuum PDE."
    )
    k_shell: int = 4
    nu: float = 0.1
    energy: float = 0.5
    t_end: float = 0.02
    n_max: int = 16
    proposed_bound_M: float = 18.0
    parent: str = "C-0002"
    statement: str = ""
    clay_implication: str = (
        "None. Shell-restricted IC class; does not address Clay Routes A-D."
    )
    numerical_support_max: float | None = None
    refuted: bool = False
    refutation_value: float | None = None
    notes: str = ""

    def as_dict(self) -> dict:
        return asdict(self)

    def ensure_statement(self) -> None:
        self.statement = (
            f"For all divergence-free Fourier ICs on N^3 (N≤{self.n_max}) with "
            f"|k|_∞ ≤ {self.k_shell}, kinetic energy = {self.energy}, evolved with "
            f"dealiased ETD-RK2 (ν={self.nu}, dt=1e-3) to T={self.t_end}, "
            f"the enstrophy Ω(T) ≤ M with M={self.proposed_bound_M}. "
            f"FINITE-DIMENSIONAL; IC shell only (flow may excite higher modes)."
        )


def save_cs0001(c: ConjectureCS0001, path: str | Path) -> None:
    c.ensure_statement()
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(c.as_dict(), indent=2), encoding="utf-8")
