"""
L-0042: Strengthenings beyond Sym Shor on consecutive {1..6} (techo).

L-0040/41: C_Sym({1..5})≤C_† but C_Sym({1..6})≈11.344>C_†≈9.562.
Rank-1 / low-rank PSD probes on {1..6} give C≲1.3≪C_†, so the obstruction
is the Sym relaxation gap, not large physical stretch.

Attempted strengthenings (all fail to prove C≤C_† on {1..6}):
  1. Column-triangle majorant from Sym column norms → C_tri≈80 ≫ C_Sym.
  2. A/B block hybrid with A={1..5}, B={6} (triangle on Z_AA/Z_BB/Z_AB)
     → C_hyb≈16.21 > C_Sym.
  3. Rank-k PSD Rayleigh lower bounds remain ≪C_† (do not certify an upper bound).

So Sym-restricted Shor remains the best *proved* majorant here, and it does
not close consecutive {1..6}. Need a true rank-constrained / SOS upper bound
(or a different identity) to absorb shell 6 into the consecutive low band.

FINITE dealias Galerkin only. Not continuum. Not Clay.
Evidence: N7 (block/triangle SVD arithmetic) + N2 (rank-k probes).
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
from ns_exploration.conjectures.l0027_low_slab_cubic import lemma_l0027
from ns_exploration.conjectures.l0040_sym_shor import _sym_col_index, sym_band_C_shor
from ns_exploration.spectral.dealias import dealias_mask
from ns_exploration.spectral.fourier_conventions import wave_number_grids

BAND_123456 = (1, 2, 3, 4, 5, 6)
BAND_12345 = (1, 2, 3, 4, 5)


def _build_sym_Lmat(n: int, radii: tuple[int, ...]) -> tuple[np.ndarray, np.ndarray]:
    by = shell_modes_by_r(n)
    modes: list[tuple[int, int, int]] = []
    shell_of: list[int] = []
    for r in radii:
        for k in by[r]:
            modes.append(k)
            shell_of.append(r)
    basis: list[tuple[np.ndarray, np.ndarray, float]] = []
    basis_shell: list[int] = []
    for k, r in zip(modes, shell_of):
        lam = float(k[0] * k[0] + k[1] * k[1] + k[2] * k[2])
        for p in _pol_basis(np.array(k, float)):
            basis.append((np.array(k, int), np.asarray(p, float), lam))
            basis_shell.append(r)
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
    return Lmat, np.asarray(basis_shell, dtype=int)


def column_triangle_C(n: int, radii: tuple[int, ...]) -> dict:
    Lmat, _ = _build_sym_Lmat(n, radii)
    D = Lmat.shape[0]
    alpha = np.linalg.norm(Lmat, axis=0)
    A = np.zeros((D, D), dtype=np.float64)
    for i in range(D):
        A[i, i] = alpha[i]
    c = D
    for i in range(D):
        for j in range(i + 1, D):
            A[i, j] = A[j, i] = math.sqrt(0.5) * alpha[c]
            c += 1
    lam = float(np.linalg.eigvalsh(A)[-1])
    C = 2.0 * math.sqrt(2.0) * lam
    Lop = float(np.linalg.svd(Lmat, compute_uv=False)[0])
    return {
        "D": D,
        "lam_A": lam,
        "C_triangle": C,
        "C_sym": 2.0 * math.sqrt(2.0) * Lop,
    }


def ab_block_hybrid_C(n: int = 24) -> dict:
    Lmat, shells = _build_sym_Lmat(n, BAND_123456)
    D = Lmat.shape[0]
    setA, setB = {1, 2, 3, 4, 5}, {6}
    idxA = np.array([i for i in range(D) if int(shells[i]) in setA], dtype=int)
    idxB = np.array([i for i in range(D) if int(shells[i]) in setB], dtype=int)

    def col_for(i: int, j: int) -> int:
        return i if i == j else _sym_col_index(i, j, D)

    def cols_block(I: np.ndarray, J: np.ndarray) -> list[int]:
        cols: set[int] = set()
        for i in I:
            for j in J:
                ii, jj = int(i), int(j)
                cols.add(col_for(ii, jj) if ii <= jj else col_for(jj, ii))
        return sorted(cols)

    def op(rows: np.ndarray, cols: list[int]) -> float:
        if len(rows) == 0 or not cols:
            return 0.0
        M = Lmat[np.ix_(rows, np.asarray(cols, dtype=int))]
        return float(np.linalg.svd(M, compute_uv=False)[0])

    cols_AA = cols_block(idxA, idxA)
    cols_BB = cols_block(idxB, idxB)
    cols_AB = cols_block(idxA, idxB)
    ops = {
        "L_A_from_AA": op(idxA, cols_AA),
        "L_B_from_BB": op(idxB, cols_BB),
        "L_A_from_BB": op(idxA, cols_BB),
        "L_B_from_AA": op(idxB, cols_AA),
        "L_A_from_AB": op(idxA, cols_AB),
        "L_B_from_AB": op(idxB, cols_AB),
    }
    laa, lbb = ops["L_A_from_AA"], ops["L_B_from_BB"]
    lab_b, lba_a = ops["L_A_from_BB"], ops["L_B_from_AA"]
    la_ab, lb_ab = ops["L_A_from_AB"], ops["L_B_from_AB"]

    def stretch_ub(u: float) -> float:
        v = 1.0 - u
        za, zb = math.sqrt(u), math.sqrt(v)
        zab = math.sqrt(2.0) * za * zb
        return za * (laa * u + lab_b * v + la_ab * zab) + zb * (
            lbb * v + lba_a * u + lb_ab * zab
        )

    us = np.linspace(0.0, 1.0, 201)
    vals = [stretch_ub(float(u)) for u in us]
    i = int(np.argmax(vals))
    max_s = float(vals[i])
    return {
        **ops,
        "max_stretch_unit": max_s,
        "C_hybrid": 2.0 * math.sqrt(2.0) * max_s,
        "u_star": float(us[i]),
    }


def rank1_C_lower(n: int, radii: tuple[int, ...], trials: int = 40) -> float:
    """N2 lower bound on true cubic C via tensor power iteration."""
    Lmat, _ = _build_sym_Lmat(n, radii)
    D = Lmat.shape[0]
    # Reconstruct apply on zz^T via Sym coordinates is awkward; use pair apply.
    by = shell_modes_by_r(n)
    modes = []
    for r in radii:
        modes.extend(list(by[r]))
    basis = []
    for k in modes:
        lam = float(k[0] * k[0] + k[1] * k[1] + k[2] * k[2])
        for p in _pol_basis(np.array(k, float)):
            basis.append((np.array(k, int), np.asarray(p, float), lam))
    mask = dealias_mask(n)
    kx, ky, kz = wave_number_grids(n)
    out: dict[tuple[int, int, int], int] = {}
    for ix, iy, iz in np.argwhere(mask & ~((kx == 0) & (ky == 0) & (kz == 0))):
        s = (int(kx[ix, iy, iz]), int(ky[ix, iy, iz]), int(kz[ix, iy, iz]))
        out[s] = len(out)
    n_out = len(out)
    I: list[int] = []
    J: list[int] = []
    SI: list[int] = []
    VX: list[float] = []
    VY: list[float] = []
    VZ: list[float] = []
    for i, (ki, ai, li) in enumerate(basis):
        for j, (kj, aj, lj) in enumerate(basis):
            s = (int(ki[0] + kj[0]), int(ki[1] + kj[1]), int(ki[2] + kj[2]))
            if s not in out:
                continue
            coef = float(np.dot(ai, kj.astype(float))) / math.sqrt(li * lj)
            Pv = _leray_vec(np.array(s, float), aj)
            I.append(i)
            J.append(j)
            SI.append(out[s])
            VX.append(coef * float(Pv[0]))
            VY.append(coef * float(Pv[1]))
            VZ.append(coef * float(Pv[2]))
    Ia = np.asarray(I, np.int32)
    Ja = np.asarray(J, np.int32)
    SIa = np.asarray(SI, np.int32)
    VXa, VYa, VZa = map(np.asarray, (VX, VY, VZ))
    by_s: dict[int, list[tuple[int, np.ndarray, float]]] = defaultdict(list)
    for a, (k, p, lam) in enumerate(basis):
        by_s[out[tuple(int(x) for x in k)]].append((a, p, lam))

    def apply_L_z(z: np.ndarray) -> np.ndarray:
        Z = np.outer(z, z)
        field = np.zeros((n_out, 3))
        np.add.at(field[:, 0], SIa, Z[Ia, Ja] * VXa)
        np.add.at(field[:, 1], SIa, Z[Ia, Ja] * VYa)
        np.add.at(field[:, 2], SIa, Z[Ia, Ja] * VZa)
        outv = np.zeros(D)
        for si, lst in by_s.items():
            Ns = field[si]
            for a, p, lam in lst:
                outv[a] += math.sqrt(lam) * float(np.dot(p, Ns))
        return outv

    rng = np.random.default_rng(0)
    best = 0.0
    for _ in range(trials):
        z = rng.normal(size=D)
        z /= np.linalg.norm(z) + 1e-30
        for _it in range(50):
            g = apply_L_z(z)
            nrm = np.linalg.norm(g)
            if nrm < 1e-30:
                break
            z = g / nrm
        best = max(best, abs(float(np.dot(z, apply_L_z(z)))))
    return 2.0 * math.sqrt(2.0) * best


@dataclass
class GalerkinBoundL0042:
    lemma_id: str = "L-0042"
    route: str = "B"
    status: str = "proved_restricted"
    evidence_level: str = "N7"
    n: int = 24
    c0007_M: float = 0.0
    C_dagger: float = 0.0
    C_sym_12345: float = 0.0
    C_sym_123456: float = 0.0
    C_triangle_123456: float = 0.0
    C_hybrid_AB: float = 0.0
    C_rank1_lo_123456: float = 0.0
    closes_123456: bool = False
    closes_c0007_all_ic: bool = False
    clay_implication: str = (
        "None. Techo on strengthenings of Sym Shor for {1..6}; "
        "not continuum; not Clay."
    )
    notes: str = ""

    def as_dict(self) -> dict:
        return asdict(self)


def lemma_l0042(n: int = 24) -> GalerkinBoundL0042:
    b7 = lemma_l0027(n=n, empirical=False)
    Cd = b7.C_dagger
    c5 = sym_band_C_shor(n, BAND_12345)["C_shor_sym"]
    c6 = sym_band_C_shor(n, BAND_123456)["C_shor_sym"]
    tri = column_triangle_C(n, BAND_123456)
    hy = ab_block_hybrid_C(n)
    r1 = rank1_C_lower(n, BAND_123456, trials=30)
    assert c6 > Cd
    assert tri["C_triangle"] > c6
    assert hy["C_hybrid"] > Cd
    notes = (
        f"{{1..6}} Sym C={c6:.4g}>C_†; triangle C={tri['C_triangle']:.4g}; "
        f"A/B hybrid C={hy['C_hybrid']:.4g}; rank1_lo C≳{r1:.4g}. "
        f"No strengthening closes {{1..6}}. All-IC: False."
    )
    return GalerkinBoundL0042(
        n=n,
        c0007_M=b7.c0007_M,
        C_dagger=Cd,
        C_sym_12345=float(c5),
        C_sym_123456=float(c6),
        C_triangle_123456=float(tri["C_triangle"]),
        C_hybrid_AB=float(hy["C_hybrid"]),
        C_rank1_lo_123456=float(r1),
        closes_123456=False,
        closes_c0007_all_ic=False,
        notes=notes,
    )


def save_lemma_l0042(
    bound: GalerkinBoundL0042,
    path: str | Path = "conjectures/proved_restricted/L-0042.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(bound.as_dict(), indent=2), encoding="utf-8")
