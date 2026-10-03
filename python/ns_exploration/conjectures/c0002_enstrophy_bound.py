"""
C-0002: Reformulated finite-dimensional enstrophy bound after C-0001 refutation.

Same setting as C-0001 but with a larger explicit M informed by gated adversarial
search (CMA + adjoint polish). Still Galerkin/pseudospectral only — not Clay.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path


@dataclass
class ConjectureC0002:
    id: str = "C-0002"
    route: str = "B"
    status: str = "proposed"
    evidence_level: str = "N6"
    domain: str = "Galerkin/pseudospectral T^3, N<=16 (NOT continuum PDE)"
    nu: float = 0.1
    energy: float = 0.5
    t_end: float = 0.02
    n_max: int = 16
    proposed_bound_M: float = 20.0
    parent_refuted: str = "C-0001"
    statement: str = ""
    clay_implication: str = (
        "None. Finite-dimensional only; no continuum regularity or blow-up claim."
    )
    numerical_support_max: float | None = None
    refuted: bool = False
    refutation_value: float | None = None
    notes: str = ""

    def as_dict(self) -> dict:
        return asdict(self)

    def __post_init__(self):
        if not self.statement:
            self.statement = (
                f"For all divergence-free Fourier fields on the N^3 grid with N<={self.n_max}, "
                f"dealiased convective term, nu={self.nu}, kinetic energy = {self.energy}, evolved with "
                f"ETD-RK2 (dt=1e-3) to T={self.t_end}, the enstrophy E_ens(T) <= M "
                f"with M={self.proposed_bound_M}. FINITE-DIMENSIONAL only. "
                f"Replaces refuted C-0001 (M≈10.99)."
            )


def save_conjecture_c0002(c: ConjectureC0002, path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if not c.statement:
        c.__post_init__()
    path.write_text(json.dumps(c.as_dict(), indent=2), encoding="utf-8")
