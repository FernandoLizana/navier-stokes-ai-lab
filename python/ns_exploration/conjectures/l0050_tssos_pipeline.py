"""
L-0050: TSSOS+Clarabel pipeline for physical cubic (beats Shor on small bands).

On Hermitian support {1,2,3} (D=52), order-2 TSSOS with sphere equality and
Clarabel yields

  C_ub ≈ 0.844 ≪ C_fullsym ≈ 1.370 ≪ C_† ≈ 9.562.

So the SOS hierarchy can strictly improve on physical fullsym Shor (L-0044/48).
Evidence is float SDP (N2) until Gram rationalization. Full dealias still open.

FINITE Galerkin only. Not continuum. Not Clay.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path

from ns_exploration.conjectures.l0027_low_slab_cubic import lemma_l0027
from ns_exploration.conjectures.l0048_fullsym_techo import C_FULLSYM_FULL_DEALIAS

# Frozen from tools/sos_julia/data/band_123/tssos_result.json (TSSOS 1.5.3 + Clarabel)
C_UB_BAND_123 = 0.843822489307523
C_FULLSYM_BAND_123 = 1.36986977843755
D_BAND_123 = 52


@dataclass
class GalerkinBoundL0050:
    lemma_id: str = "L-0050"
    route: str = "B"
    status: str = "exploring"
    evidence_level: str = "N2"
    n: int = 24
    C_dagger: float = 0.0
    band_radii: list[int] | None = None
    D: int = D_BAND_123
    C_ub_tssos: float = C_UB_BAND_123
    C_fullsym_band: float = C_FULLSYM_BAND_123
    C_fullsym_full_dealias: float = C_FULLSYM_FULL_DEALIAS
    beats_fullsym_on_band: bool = True
    closes_c0007_all_ic: bool = False
    clay_implication: str = (
        "None. Float TSSOS on a tiny band; not all-IC; not continuum; not Clay."
    )
    notes: str = ""

    def as_dict(self) -> dict:
        return asdict(self)


def lemma_l0050(n: int = 24) -> GalerkinBoundL0050:
    Cd = lemma_l0027(n=n, empirical=False).C_dagger
    assert C_UB_BAND_123 < C_FULLSYM_BAND_123
    assert C_UB_BAND_123 < Cd
    notes = (
        f"TSSOS+Clarabel order-2 on {{1,2,3}}: C_ub≈{C_UB_BAND_123:.4f} "
        f"< C_fullsym≈{C_FULLSYM_BAND_123:.4f} < C_†≈{Cd:.4f}. "
        f"Pipeline ready; full-dealias Shor techo C≈{C_FULLSYM_FULL_DEALIAS:.2f} "
        "still needs scaled SOS / structure."
    )
    return GalerkinBoundL0050(
        n=n,
        C_dagger=Cd,
        band_radii=[1, 2, 3],
        notes=notes,
    )


def save_lemma_l0050(
    bound: GalerkinBoundL0050,
    path: str | Path = "conjectures/active/L-0050.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(bound.as_dict(), indent=2), encoding="utf-8")
