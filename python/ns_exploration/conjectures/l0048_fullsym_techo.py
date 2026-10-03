"""
L-0048: Sparse physical fullsym on full dealias — Shor techo for all-IC.

Streaming dense M cannot fit full Hermitian dealias (D=6748, ~0.6 TB).
M is sparse (nnz≈3.78e7). Assembling CSR + Gram σmax yields

  C_fullsym ≈ 25.926 ≫ C_† ≈ 9.562.

Therefore the physical fullsym Shor majorant on the full dealias mask
cannot certify the low-slab hypothesis of L-0027 for all-IC C-0007.

This does NOT refute C-0007: C_true ≤ C_fullsym, so C_fullsym > C_† is
only a certificate-method techo (cf. L-0036 for ambient/non-Hermitian Shor).
Subclasses C-R-0009…0012 remain valid (restricted supports with C≤C_†).

FINITE dealias Galerkin only. Not continuum. Not Clay.
Evidence: N7 (sparse CSR Gram of exact triad tensor; validated vs dense on {1..6}).
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path

from ns_exploration.conjectures.l0026_c0007_techo import lemma_l0026
from ns_exploration.conjectures.l0027_low_slab_cubic import lemma_l0027
from ns_exploration.conjectures.l0048_matrixfree_fullsym import (
    all_dealias_radii,
    sparse_C_fullsym,
    validate_sparse_vs_dense,
)

# Frozen witness from reports/l0048_full_dealias_sparse.json
C_FULLSYM_FULL_DEALIAS = 25.925922025377965
NNZ_FULL_DEALIAS = 37818119
D_FULL_DEALIAS = 6748


@dataclass
class GalerkinBoundL0048:
    lemma_id: str = "L-0048"
    route: str = "B"
    status: str = "techo"
    evidence_level: str = "N7"
    n: int = 24
    c0007_M: float = 0.0
    C_dagger: float = 0.0
    Omega_star: float = 0.0
    D_full: int = D_FULL_DEALIAS
    nnz_full: int = NNZ_FULL_DEALIAS
    C_fullsym_full: float = C_FULLSYM_FULL_DEALIAS
    validate_rel_err: float = 0.0
    closes_c0007_all_ic: bool = False
    fullsym_shor_all_ic_techo: bool = True
    clay_implication: str = (
        "None. Full-dealias physical fullsym Shor exceeds C_† — method techo "
        "for all-IC via this majorant; not a refutation of C-0007; not Clay."
    )
    notes: str = ""

    def as_dict(self) -> dict:
        return asdict(self)


def lemma_l0048(n: int = 24, recompute: bool = False) -> GalerkinBoundL0048:
    b7 = lemma_l0027(n=n, empirical=False)
    b6 = lemma_l0026(n=n, n_grid=21)
    Cd = b7.C_dagger
    val = validate_sparse_vs_dense(n=n)
    assert val["ok"]
    if recompute:
        full = sparse_C_fullsym(n, all_dealias_radii(n), progress_every=500_000)
        C_full = float(full["C_fullsym"])
        D = int(full["D"])
        nnz = int(full["nnz"])
    else:
        C_full = C_FULLSYM_FULL_DEALIAS
        D = D_FULL_DEALIAS
        nnz = NNZ_FULL_DEALIAS
    assert C_full > Cd
    assert D == D_FULL_DEALIAS
    notes = (
        f"Sparse CSR Gram full dealias: C={C_full:.4f}>>C_dagger={Cd:.4f} "
        f"(D={D}, nnz={nnz}). Fullsym Shor cannot close all-IC. "
        f"Validate {{1..6}} rel={val['rel_err']:.3e}."
    )
    return GalerkinBoundL0048(
        n=n,
        c0007_M=b7.c0007_M,
        C_dagger=Cd,
        Omega_star=b6.Omega_star,
        D_full=D,
        nnz_full=nnz,
        C_fullsym_full=C_full,
        validate_rel_err=float(val["rel_err"]),
        closes_c0007_all_ic=False,
        fullsym_shor_all_ic_techo=True,
        notes=notes,
    )


def save_lemma_l0048(
    bound: GalerkinBoundL0048,
    path: str | Path = "conjectures/proved_restricted/L-0048.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(bound.as_dict(), indent=2), encoding="utf-8")
