"""
Sprint A (C-0007 cuello de botella): auditoría del tensor cúbico de stretch +
simetrización total + comparación C_fullsym vs C_Sym.

NO modifica lemmas. NO instala Julia. NO afirma nada continuo/Clay.
Objetivo: certificar que el objeto de la cota Sym (l0040) es exactamente la forma
cúbica del vortex stretching, y comprobar si la matricización totalmente simétrica
(flattening del 3-tensor simetrizado) baja la norma de operador por debajo de C_Sym.

Definición (base real transversal, coords z ∈ R^D):
    f(z) = Σ_{a,i,j} G[a,i,j] z_a z_i z_j
    G[a,i,j] = |k_a| · <p_a, P_{k_i+k_j}(p_j)> · <p_i, k_j>   si  k_a = k_i + k_j
donde P_s = proyector de Leray en s.  Esto reproduce z·(Lmat·svec(zz^T)) de l0040.

C = 2√2 · sup_{||z||=1} |f(z)| ,  y se necesita  C ≤ C_† ≈ 9.562.
Cotas superiores (matricización → norma de operador):
    C_Sym     = 2√2 ||M_sym_ij||_op      (l0040; simétrica solo en (i,j))
    C_fullsym = 2√2 ||M_fullsym||_op      (3-tensor totalmente simétrico)
Ambas mayoran el mismo sup rank-1; C_fullsym ≤ C_Sym siempre.

Evidencia: N7 (identidad de tensor exacta) + N2 (chequeo pseudoespectral).
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
from ns_exploration.conjectures.l0040_sym_shor import _sym_col_index
from ns_exploration.conjectures.l0043_onepol_shellblock import (
    _build_twopol_Lmat,
    sym_onepol_C_shor,
)
from ns_exploration.spectral.dealias import dealias_mask
from ns_exploration.spectral.fourier_conventions import wave_number_grids

BAND_123456 = (1, 2, 3, 4, 5, 6)


@dataclass
class StretchBasis:
    """Real transverse basis over the retained modes of a band."""

    ks: np.ndarray  # (D, 3) int wavevectors
    pols: np.ndarray  # (D, 3) real unit transverse vectors
    lams: np.ndarray  # (D,) |k|^2
    shells: np.ndarray  # (D,) shell index
    D: int


def build_basis(n: int, radii: tuple[int, ...]) -> StretchBasis:
    by = shell_modes_by_r(n)
    ks: list[list[int]] = []
    pols: list[np.ndarray] = []
    lams: list[float] = []
    shells: list[int] = []
    for r in radii:
        if r not in by:
            continue
        for k in by[r]:
            lam = float(k[0] * k[0] + k[1] * k[1] + k[2] * k[2])
            for p in _pol_basis(np.array(k, float)):
                ks.append([int(k[0]), int(k[1]), int(k[2])])
                pols.append(np.asarray(p, float))
                lams.append(lam)
                shells.append(int(r))
    return StretchBasis(
        ks=np.array(ks, dtype=int),
        pols=np.array(pols, dtype=float),
        lams=np.array(lams, dtype=float),
        shells=np.array(shells, dtype=int),
        D=len(ks),
    )


def build_onepol_basis(
    n: int,
    radii: tuple[int, ...],
    branch: str = "first",
) -> StretchBasis:
    """One polarization per mode (matches l0043 sym_onepol mode list)."""
    if branch not in ("first", "second"):
        raise ValueError(f"unknown pol branch {branch!r}")
    by = shell_modes_by_r(n)
    ks: list[list[int]] = []
    pols: list[np.ndarray] = []
    lams: list[float] = []
    shells: list[int] = []
    for r in radii:
        if r not in by:
            continue
        for k in by[r]:
            lam = float(k[0] * k[0] + k[1] * k[1] + k[2] * k[2])
            pb = _pol_basis(np.array(k, float))
            if branch == "first":
                p = pb[0]
            else:
                p = pb[1] if len(pb) > 1 else pb[0]
            ks.append([int(k[0]), int(k[1]), int(k[2])])
            pols.append(np.asarray(p, float))
            lams.append(lam)
            shells.append(int(r))
    return StretchBasis(
        ks=np.array(ks, dtype=int),
        pols=np.array(pols, dtype=float),
        lams=np.array(lams, dtype=float),
        shells=np.array(shells, dtype=int),
        D=len(ks),
    )


def _build_G_from_basis(b: StretchBasis, n: int) -> np.ndarray:
    D = b.D
    mask = dealias_mask(n)
    kx, ky, kz = wave_number_grids(n)
    out: dict[tuple[int, int, int], list[int]] = defaultdict(list)
    for a in range(D):
        out[(int(b.ks[a, 0]), int(b.ks[a, 1]), int(b.ks[a, 2]))].append(a)
    G = np.zeros((D, D, D), dtype=np.float64)
    for i in range(D):
        ki = b.ks[i]
        pi = b.pols[i]
        for j in range(D):
            kj = b.ks[j]
            s = (int(ki[0] + kj[0]), int(ki[1] + kj[1]), int(ki[2] + kj[2]))
            rows = out.get(s)
            if not rows:
                continue
            coef_adv = float(np.dot(pi, kj.astype(float))) / math.sqrt(
                b.lams[i] * b.lams[j]
            )
            if coef_adv == 0.0:
                continue
            Pv = _leray_vec(np.array(s, float), b.pols[j])
            for a in rows:
                G[a, i, j] += (
                    math.sqrt(b.lams[a]) * float(np.dot(b.pols[a], Pv)) * coef_adv
                )
    return G


def build_G(n: int, radii: tuple[int, ...]) -> tuple[np.ndarray, StretchBasis]:
    """Dense cubic tensor G[a,i,j] (not symmetric). f(z)=einsum('aij,a,i,j')."""
    b = build_basis(n, radii)
    return _build_G_from_basis(b, n), b


def build_onepol_G(
    n: int,
    radii: tuple[int, ...],
    branch: str = "first",
) -> tuple[np.ndarray, StretchBasis]:
    """Cubic stretch tensor on one-pol Galerkin subclass (l0043 / C-R-0008)."""
    b = build_onepol_basis(n, radii, branch=branch)
    return _build_G_from_basis(b, n), b


def f_from_G(G: np.ndarray, z: np.ndarray) -> float:
    return float(np.einsum("aij,a,i,j->", G, z, z, z, optimize=True))


def f_from_Lmat(Lmat: np.ndarray, z: np.ndarray) -> float:
    """z·(Lmat·svec(zz^T)) with l0040 convention (diag=z_i^2, off=√2 z_i z_j)."""
    D = Lmat.shape[0]
    x = np.empty(Lmat.shape[1], dtype=np.float64)
    x[:D] = z * z
    c = D
    for i in range(D):
        for jj in range(i + 1, D):
            x[c] = math.sqrt(2.0) * z[i] * z[jj]
            c += 1
    return float(z @ (Lmat @ x))


def matricize(H: np.ndarray) -> np.ndarray:
    """Flatten a 3-tensor (symmetric in its last two indices) to M: R^{D(D+1)/2}→R^D."""
    D = H.shape[0]
    n_sym = D * (D + 1) // 2
    M = np.zeros((D, n_sym), dtype=np.float64)
    M[:, :D] = H[:, np.arange(D), np.arange(D)]
    sq2 = math.sqrt(2.0)
    c = D
    for i in range(D):
        for jj in range(i + 1, D):
            M[:, c] = sq2 * H[:, i, jj]
            c += 1
    return M


def sym_ij(G: np.ndarray) -> np.ndarray:
    return 0.5 * (G + np.transpose(G, (0, 2, 1)))


def full_sym(G: np.ndarray) -> np.ndarray:
    perms = [
        (0, 1, 2),
        (0, 2, 1),
        (1, 0, 2),
        (1, 2, 0),
        (2, 0, 1),
        (2, 1, 0),
    ]
    A = np.zeros_like(G)
    for p in perms:
        A += np.transpose(G, p)
    return A / 6.0


def operator_C(H: np.ndarray) -> float:
    M = matricize(H)
    return 2.0 * math.sqrt(2.0) * float(np.linalg.svd(M, compute_uv=False)[0])


def rank1_lower_C(A: np.ndarray, n_restart: int = 120, iters: int = 100) -> float:
    """N2 lower bound: 2√2·max_{‖z‖=1}|f(z)| via tensor power ascent on sym A.

    grad f(z) ∝ A(z,z,·). Vectorized as A2 @ (z⊗z) with A2 = A reshaped (D, D·D).
    Multi-restart on the sphere.
    """
    D = A.shape[0]
    A2 = A.reshape(D, D * D)
    rng = np.random.default_rng(7)
    best = 0.0
    for _ in range(n_restart):
        z = rng.standard_normal(D)
        z /= np.linalg.norm(z)
        for _ in range(iters):
            g = A2 @ np.outer(z, z).ravel()  # A(z,z,·)
            gn = np.linalg.norm(g)
            if gn < 1e-30:
                break
            z_new = g / gn
            if abs(float(z_new @ z)) > 1.0 - 1e-14:
                z = z_new
                break
            z = z_new
        val = abs(float(z @ (A2 @ np.outer(z, z).ravel())))
        if val > best:
            best = val
    return 2.0 * math.sqrt(2.0) * best


@dataclass
class SprintAResult:
    n: int
    radii: list[int]
    D: int
    C_dagger: float
    max_tensor_vs_lmat_err: float
    max_fullsym_vs_G_err: float
    C_sym_reproduced: float
    C_sym_l0040: float
    C_fullsym: float
    C_rank1_lower: float
    relaxation_gap_ratio: float
    fullsym_le_sym: bool
    fullsym_closes_band: bool
    notes: str

    def as_dict(self) -> dict:
        return asdict(self)


def run(n: int = 24, radii: tuple[int, ...] = BAND_123456, n_probe: int = 1500) -> SprintAResult:
    Cd = lemma_l0027(n=n, empirical=False).C_dagger
    G, b = build_G(n, radii)
    D = b.D
    Lmat, _ = _build_twopol_Lmat(n, radii)

    # Test (A): tensor internal consistency vs established l0040 matricization.
    rng = np.random.default_rng(20260718)
    err_lmat = 0.0
    err_fs = 0.0
    Gsym = sym_ij(G)
    A = full_sym(G)
    for _ in range(n_probe):
        z = rng.standard_normal(D)
        fg = f_from_G(Gsym, z)
        fl = f_from_Lmat(Lmat, z)
        scale = max(1.0, abs(fg), abs(fl))
        err_lmat = max(err_lmat, abs(fg - fl) / scale)
        fa = f_from_G(A, z)
        err_fs = max(err_fs, abs(fa - fg) / scale)

    C_sym_rep = operator_C(Gsym)
    C_fullsym = operator_C(A)
    C_sym_l0040 = 2.0 * math.sqrt(2.0) * float(
        np.linalg.svd(Lmat, compute_uv=False)[0]
    )
    C_r1 = rank1_lower_C(A)
    gap = C_fullsym / C_r1 if C_r1 > 0 else float("inf")
    notes = (
        f"tensor==l0040 err={err_lmat:.2e}; fullsym==G err={err_fs:.2e}; "
        f"C_sym_repro={C_sym_rep:.6f} (l0040 {C_sym_l0040:.6f}); "
        f"C_fullsym={C_fullsym:.6f}; C_rank1_lo={C_r1:.6f}; C_†={Cd:.6f}. "
        f"relaxation gap ≈{gap:.1f}×. "
        f"fullsym{'<=' if C_fullsym <= C_sym_rep + 1e-9 else '>'}sym; "
        f"closes_band={C_fullsym <= Cd}."
    )
    return SprintAResult(
        n=n,
        radii=list(radii),
        D=D,
        C_dagger=Cd,
        max_tensor_vs_lmat_err=err_lmat,
        max_fullsym_vs_G_err=err_fs,
        C_sym_reproduced=C_sym_rep,
        C_sym_l0040=C_sym_l0040,
        C_fullsym=C_fullsym,
        C_rank1_lower=C_r1,
        relaxation_gap_ratio=float(gap),
        fullsym_le_sym=bool(C_fullsym <= C_sym_rep + 1e-9),
        fullsym_closes_band=bool(C_fullsym <= Cd),
        notes=notes,
    )


def main() -> dict:
    res = run()
    out = res.as_dict()
    Path("reports").mkdir(exist_ok=True)
    Path("experiments/exploratory/sprintA_stretch_tensor").mkdir(
        parents=True, exist_ok=True
    )
    Path("experiments/exploratory/sprintA_stretch_tensor/summary.json").write_text(
        json.dumps(out, indent=2), encoding="utf-8"
    )
    print(json.dumps(out, indent=2, ensure_ascii=True))
    return out


if __name__ == "__main__":
    main()
