"""
L-0039: Sparse-SVD shell-support extension + C-R-0005.

L-0038 closed C-0007 on Fourier support in {1,2,3} via dense SVD Shor.
For larger supports, build the stretch matricization L as a sparse D×(#pairs)
matrix and compute ‖L‖_op by sparse SVD (SciPy svds). Then
  C ≤ C_Shor := 2√2 ‖L‖_op
as in L-0036/38.

Greedy growth (ascending shell order) from L-0038 seeds yields supports with
C_Shor ≤ C_† ≈ 9.562, e.g.
  A = {1,2,3,10,12,19,27,35,43,58,69,73,123,147}  → C_Shor ≈ 9.519
  B = {1,3,4,6,8,17,32,43,58,76,85,102,108,134}   → C_Shor ≈ 9.548
Also the 4-shell set {1,3,4,6} has C_Shor ≈ 8.485 ≤ C_† (dense-checked).

C-R-0005: C-0007 on support A via L-0039 + L-0026 + L-0027.
All-IC still open (consecutive {1..4} and many additions fail C_†).

FINITE dealias Galerkin only. Not continuum. Not Clay.
Evidence: N7 (sparse SVD on finite matrix) + L-0026/L-0027.
"""

from __future__ import annotations

import json
import math
from collections import defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import svds

from ns_exploration.conjectures.l0021_quartic_shell import (
    _leray_vec,
    _pol_basis,
    shell_modes_by_r,
)
from ns_exploration.conjectures.l0026_c0007_techo import lemma_l0026
from ns_exploration.conjectures.l0027_low_slab_cubic import lemma_l0027
from ns_exploration.conjectures.l0038_band_svd_shor import exact_band_C_shor
from ns_exploration.spectral.dealias import dealias_mask
from ns_exploration.spectral.fourier_conventions import wave_number_grids

# Greedy-maximal supports (ascending shell order; frozen witnesses).
SUPPORT_A = (1, 2, 3, 10, 12, 19, 27, 35, 43, 58, 69, 73, 123, 147)
SUPPORT_B = (1, 3, 4, 6, 8, 17, 32, 43, 58, 76, 85, 102, 108, 134)
QUAD_1346 = (1, 3, 4, 6)


def sparse_band_C_shor(n: int, radii: tuple[int, ...] | list[int]) -> dict:
    """C_Shor via sparse SVD of the Ω-weighted stretch matricization."""
    by = shell_modes_by_r(n)
    radii_t = tuple(int(r) for r in radii)
    modes: list[tuple[int, int, int]] = []
    for r in radii_t:
        if r not in by:
            raise KeyError(f"shell {r} missing on n={n}")
        modes.extend(list(by[r]))
    basis: list[tuple[np.ndarray, np.ndarray, float]] = []
    for k in modes:
        lam = float(k[0] * k[0] + k[1] * k[1] + k[2] * k[2])
        for p in _pol_basis(np.array(k, float)):
            basis.append((np.array(k, int), np.asarray(p, float), lam))
    D = len(basis)
    mask = dealias_mask(n)
    kx, ky, kz = wave_number_grids(n)
    out: dict[tuple[int, int, int], int] = {}
    for ix, iy, iz in np.argwhere(mask & ~((kx == 0) & (ky == 0) & (kz == 0))):
        s = (int(kx[ix, iy, iz]), int(ky[ix, iy, iz]), int(kz[ix, iy, iz]))
        out[s] = len(out)
    by_s: dict[int, list[tuple[int, np.ndarray, float]]] = defaultdict(list)
    for a, (k, p, lam) in enumerate(basis):
        by_s[out[tuple(int(x) for x in k)]].append((a, p, lam))

    col_id: dict[tuple[int, int], int] = {}
    rows: list[int] = []
    cols: list[int] = []
    vals: list[float] = []
    n_pairs = 0
    for i, (ki, ai, li) in enumerate(basis):
        for j, (kj, aj, lj) in enumerate(basis):
            s = (int(ki[0] + kj[0]), int(ki[1] + kj[1]), int(ki[2] + kj[2]))
            if s not in out:
                continue
            n_pairs += 1
            key = (i, j)
            if key not in col_id:
                col_id[key] = len(col_id)
            c = col_id[key]
            coef = float(np.dot(ai, kj.astype(float))) / math.sqrt(li * lj)
            Pv = _leray_vec(np.array(s, float), aj)
            vec = coef * Pv
            for a, p, lam in by_s[out[s]]:
                rows.append(a)
                cols.append(c)
                vals.append(math.sqrt(lam) * float(np.dot(p, vec)))
    n_col = len(col_id)
    if D == 0 or n_col == 0:
        return {
            "n": n,
            "radii": list(radii_t),
            "D": D,
            "n_pairs": 0,
            "n_cols": 0,
            "L_op": 0.0,
            "C_shor": 0.0,
        }
    L = lil_matrix((D, n_col), dtype=np.float64)
    for r, c, v in zip(rows, cols, vals):
        L[r, c] += v
    Lc = L.tocsr()
    if min(D, n_col) <= 3:
        Lop = float(np.linalg.norm(Lc.toarray(), ord=2))
    else:
        svals = svds(Lc, k=1, which="LM", return_singular_vectors=False)
        Lop = float(svals[-1])
    C = 2.0 * math.sqrt(2.0) * Lop
    return {
        "n": n,
        "radii": list(radii_t),
        "D": D,
        "n_pairs": n_pairs,
        "n_cols": n_col,
        "L_op": Lop,
        "C_shor": C,
    }


@dataclass
class GalerkinBoundL0039:
    lemma_id: str = "L-0039"
    route: str = "B"
    status: str = "proved_restricted"
    evidence_level: str = "N7"
    n: int = 24
    c0007_M: float = 0.0
    C_dagger: float = 0.0
    quad_1346_C: float = 0.0
    support_A: list | None = None
    support_A_C: float = 0.0
    support_A_D: int = 0
    support_B: list | None = None
    support_B_C: float = 0.0
    support_B_D: int = 0
    closes_cr0005: bool = True
    closes_c0007_all_ic: bool = False
    clay_implication: str = (
        "None. Sparse-SVD Shor on enlarged shell supports; C-R-0005 only; "
        "not all-IC; not continuum; not Clay."
    )
    notes: str = ""

    def as_dict(self) -> dict:
        return asdict(self)


def lemma_l0039(n: int = 24) -> GalerkinBoundL0039:
    b7 = lemma_l0027(n=n, empirical=False)
    Cd = b7.C_dagger
    q = exact_band_C_shor(n, QUAD_1346)  # small: dense
    a = sparse_band_C_shor(n, SUPPORT_A)
    b = sparse_band_C_shor(n, SUPPORT_B)
    assert q["C_shor"] <= Cd + 1e-8
    assert a["C_shor"] <= Cd + 1e-6
    assert b["C_shor"] <= Cd + 1e-6
    notes = (
        f"Sparse SVD: support A ({len(SUPPORT_A)} shells) C={a['C_shor']:.4g}; "
        f"B ({len(SUPPORT_B)}) C={b['C_shor']:.4g}; quad {{1,3,4,6}} C={q['C_shor']:.4g} "
        f"(all ≤ C_†={Cd:.4g}). C-R-0005 on A. All-IC: False."
    )
    return GalerkinBoundL0039(
        n=n,
        c0007_M=b7.c0007_M,
        C_dagger=Cd,
        quad_1346_C=float(q["C_shor"]),
        support_A=list(SUPPORT_A),
        support_A_C=float(a["C_shor"]),
        support_A_D=int(a["D"]),
        support_B=list(SUPPORT_B),
        support_B_C=float(b["C_shor"]),
        support_B_D=int(b["D"]),
        closes_cr0005=True,
        closes_c0007_all_ic=False,
        notes=notes,
    )


def build_cr0005(n: int = 24) -> dict:
    b7 = lemma_l0027(n=n, empirical=False)
    b26 = lemma_l0026(n=n, n_grid=21)
    d = sparse_band_C_shor(n, SUPPORT_A)
    assert d["C_shor"] <= b7.C_dagger + 1e-6
    return {
        "id": "C-R-0005",
        "route": "B",
        "status": "proved_restricted",
        "evidence_level": "N7",
        "lemma": "L-0039+L-0026+L-0027",
        "domain": (
            "Pseudospectral T^3, N<=24, 2/3 dealias, Fourier support in the "
            f"{len(SUPPORT_A)}-shell set A={list(SUPPORT_A)}. NOT continuum."
        ),
        "nu": 0.1,
        "energy": 0.5,
        "t_end": 0.02,
        "n_max": n,
        "proposed_bound_M": b7.c0007_M,
        "shells": list(SUPPORT_A),
        "C_shor": d["C_shor"],
        "C_dagger": b7.C_dagger,
        "D": d["D"],
        "Omega_star": b26.Omega_star,
        "parent_open": "C-0007",
        "alt_support_B": list(SUPPORT_B),
        "statement": (
            f"For divergence-free fields on N≤{n} with Fourier support in shells "
            f"{list(SUPPORT_A)}, dealiased Galerkin NS (ν=0.1, E≤0.5) satisfies "
            f"Ω(0.02)≤{b7.c0007_M:.8f}. Proof: sparse SVD Shor gives "
            f"C≤{d['C_shor']:.8f}≤C_† ⇒ L-0027 low slab; L-0026 high slab. FINITE only."
        ),
        "clay_implication": (
            "None. Multi-shell support subclass only; not all-IC C-0007; not Clay."
        ),
        "notes": (
            "Greedy enlargement of C-R-0004={1,2,3}. Alternate support B also ≤C_†."
        ),
    }


def save_lemma_l0039(
    bound: GalerkinBoundL0039,
    path: str | Path = "conjectures/proved_restricted/L-0039.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(bound.as_dict(), indent=2), encoding="utf-8")


def save_cr0005(
    d: dict,
    path: str | Path = "conjectures/proved_restricted/C-R-0005.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(d, indent=2), encoding="utf-8")
