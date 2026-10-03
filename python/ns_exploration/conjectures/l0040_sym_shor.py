"""
L-0040: Sym-restricted stretch matricization + C-R-0006.

L-0036/38/39 bound |stretch| via ‖L‖_{F→2} on all matrices Z. But
stretch = z · L(zzᵀ) only samples symmetric Z=zzᵀ. Therefore
  |stretch| ≤ ‖L|_{Sym}‖_{F→2} ‖z‖³,
and C ≤ C_Sym := 2√2 ‖L|_{Sym}‖_op, which is ≤ the unrestricted Shor constant.

On N=24 dealias (exact dense SVD on an orthonormal Sym basis):
  {1..5}  → C_Sym ≈ 9.015 ≤ C_†   (consecutive K=5)
  {1..6}  → C_Sym ≈ 11.344 > C_†  (techo)
Unrestricted Shor on {1..4} was already > C_† (~10.07); Sym closes it (~6.11).

C-R-0006: Fourier support in shells {1,2,3,4,5} ⇒ C-0007 via L-0040+L-0026+L-0027.
A greedy Sym enlargement (9 shells) is recorded as an alternate witness.

FINITE dealias Galerkin only. Not continuum. Not Clay.
Evidence: N7 (Sym SVD arithmetic) + L-0026/L-0027.
"""

from __future__ import annotations

import json
import math
from collections import defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from ns_exploration.conjectures.l0021_quartic_shell import (
    _leray_vec,
    _pol_basis,
    shell_modes_by_r,
)
from ns_exploration.conjectures.l0026_c0007_techo import lemma_l0026
from ns_exploration.conjectures.l0027_low_slab_cubic import lemma_l0027
from ns_exploration.spectral.dealias import dealias_mask
from ns_exploration.spectral.fourier_conventions import wave_number_grids

# Consecutive low shells closed by Sym Shor.
SUPPORT_CR0006 = (1, 2, 3, 4, 5)
# Greedy ascending Sym enlargement (frozen witness).
SUPPORT_GREEDY9 = (1, 2, 3, 4, 5, 11, 16, 24, 25)


def _sym_col_index(i: int, j: int, D: int) -> int:
    if i == j:
        return i
    if i > j:
        i, j = j, i
    return D + i * (2 * D - i - 1) // 2 + (j - i - 1)


def sym_band_C_shor(n: int, radii: tuple[int, ...] | list[int]) -> dict:
    """C_Sym = 2√2 ‖L|_{Sym}‖_op via dense SVD."""
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

    n_sym = D * (D + 1) // 2
    Lmat = np.zeros((D, n_sym), dtype=np.float64)
    sq2 = math.sqrt(0.5)
    for i, (ki, ai, li) in enumerate(basis):
        for j, (kj, aj, lj) in enumerate(basis):
            s = (int(ki[0] + kj[0]), int(ki[1] + kj[1]), int(ki[2] + kj[2]))
            if s not in out:
                continue
            coef = float(np.dot(ai, kj.astype(float))) / math.sqrt(li * lj)
            v = coef * _leray_vec(np.array(s, float), aj)
            vec = np.zeros(D, dtype=np.float64)
            for a, p, lam in by_s[out[s]]:
                vec[a] += math.sqrt(lam) * float(np.dot(p, v))
            if i == j:
                Lmat[:, i] += vec
            else:
                Lmat[:, _sym_col_index(i, j, D)] += sq2 * vec
    Lop = float(np.linalg.svd(Lmat, compute_uv=False)[0]) if D else 0.0
    C = 2.0 * math.sqrt(2.0) * Lop
    return {
        "n": n,
        "radii": list(radii_t),
        "D": D,
        "n_sym": n_sym,
        "L_op_sym": Lop,
        "C_shor_sym": C,
    }


@dataclass
class GalerkinBoundL0040:
    lemma_id: str = "L-0040"
    route: str = "B"
    status: str = "proved_restricted"
    evidence_level: str = "N7"
    n: int = 24
    c0007_M: float = 0.0
    C_dagger: float = 0.0
    band_12345_C: float = 0.0
    band_123456_C: float = 0.0
    band_1234_C_sym: float = 0.0
    greedy9_C: float = 0.0
    greedy9_shells: list | None = None
    closes_cr0006: bool = True
    closes_c0007_all_ic: bool = False
    clay_implication: str = (
        "None. Sym-restricted Shor on low-shell bands; C-R-0006 only; "
        "not all-IC; not continuum; not Clay."
    )
    notes: str = ""

    def as_dict(self) -> dict:
        return asdict(self)


def lemma_l0040(n: int = 24) -> GalerkinBoundL0040:
    b7 = lemma_l0027(n=n, empirical=False)
    Cd = b7.C_dagger
    d5 = sym_band_C_shor(n, SUPPORT_CR0006)
    d6 = sym_band_C_shor(n, (1, 2, 3, 4, 5, 6))
    d4 = sym_band_C_shor(n, (1, 2, 3, 4))
    dg = sym_band_C_shor(n, SUPPORT_GREEDY9)
    assert d5["C_shor_sym"] <= Cd + 1e-9
    assert d6["C_shor_sym"] > Cd
    assert dg["C_shor_sym"] <= Cd + 1e-8
    notes = (
        f"Sym Shor: {{1..5}} C={d5['C_shor_sym']:.4g}≤C_†; {{1..6}} C={d6['C_shor_sym']:.4g}>C_†; "
        f"{{1..4}} C={d4['C_shor_sym']:.4g}; greedy9 C={dg['C_shor_sym']:.4g}. "
        f"C-R-0006 on {{1..5}}. All-IC: False."
    )
    return GalerkinBoundL0040(
        n=n,
        c0007_M=b7.c0007_M,
        C_dagger=Cd,
        band_12345_C=float(d5["C_shor_sym"]),
        band_123456_C=float(d6["C_shor_sym"]),
        band_1234_C_sym=float(d4["C_shor_sym"]),
        greedy9_C=float(dg["C_shor_sym"]),
        greedy9_shells=list(SUPPORT_GREEDY9),
        closes_cr0006=True,
        closes_c0007_all_ic=False,
        notes=notes,
    )


def build_cr0006(n: int = 24) -> dict:
    b7 = lemma_l0027(n=n, empirical=False)
    b26 = lemma_l0026(n=n, n_grid=21)
    d = sym_band_C_shor(n, SUPPORT_CR0006)
    assert d["C_shor_sym"] <= b7.C_dagger + 1e-9
    return {
        "id": "C-R-0006",
        "route": "B",
        "status": "proved_restricted",
        "evidence_level": "N7",
        "lemma": "L-0040+L-0026+L-0027",
        "domain": (
            "Pseudospectral T^3, N<=24, 2/3 dealias, Fourier support in shells "
            "|k|^2 in {1,2,3,4,5}. NOT continuum."
        ),
        "nu": 0.1,
        "energy": 0.5,
        "t_end": 0.02,
        "n_max": n,
        "proposed_bound_M": b7.c0007_M,
        "shells": list(SUPPORT_CR0006),
        "C_shor_sym": d["C_shor_sym"],
        "C_dagger": b7.C_dagger,
        "D": d["D"],
        "Omega_star": b26.Omega_star,
        "parent_open": "C-0007",
        "statement": (
            f"For divergence-free fields on N≤{n} with Fourier support only on "
            f"shells |k|²∈{{1,2,3,4,5}}, dealiased Galerkin NS (ν=0.1, E≤0.5) satisfies "
            f"Ω(0.02)≤{b7.c0007_M:.8f}. Proof: Sym-restricted SVD Shor gives "
            f"C≤{d['C_shor_sym']:.8f}≤C_† ⇒ L-0027; high slab via L-0026. FINITE only."
        ),
        "clay_implication": (
            "None. Consecutive low-shell subclass only; not all-IC C-0007; not Clay."
        ),
        "notes": "Unrestricted Shor on {1..4} already exceeded C_†; Sym restriction closes {1..5}.",
    }


def save_lemma_l0040(
    bound: GalerkinBoundL0040,
    path: str | Path = "conjectures/proved_restricted/L-0040.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(bound.as_dict(), indent=2), encoding="utf-8")


def save_cr0006(
    d: dict,
    path: str | Path = "conjectures/proved_restricted/C-R-0006.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(d, indent=2), encoding="utf-8")
