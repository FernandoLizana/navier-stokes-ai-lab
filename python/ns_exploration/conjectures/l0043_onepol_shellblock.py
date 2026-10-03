"""
L-0043: One-pol Sym on consecutive low shells + shell-block/Ky-Fan techos.

Context (L-0040–42): full two-pol Sym on {1..6} gives C_Sym≈11.344>C_†≈9.562;
triangle / A/B hybrid / shell-block majorants are worse; rank-1 probes ≪C_†.

New positive result: restrict each wavevector to a single frozen polarization
(the first vector of `_pol_basis`). Then Sym Shor closes a larger consecutive
band than two-pol Sym:
  one-pol {1,2,3,4,5,6,8} → C≈8.476 ≤ C_†  (shell |k|²=7 empty on Z³)
  one-pol + shell 9        → C≈9.849 > C_†  (techo)
  second-pol branch also closes {1..8} (C≈8.731) but first is the witness.

Techos packaged here (do not beat Sym on two-pol {1..6}):
  1. Shell-pair→target block majorant (sampled lower bound) ≳14.5 > C_†.
  2. Ky-Fan / SVD sum majorant Σ σ_i‖M_i‖_op already exceeds C_† after 4 terms.

C-R-0008: one-pol Fourier support on shells {1,2,3,4,5,6,8}
⇒ C-0007 via L-0043+L-0026+L-0027.

FINITE dealias Galerkin only. Not continuum. Not Clay.
Evidence: N7 (Sym SVD arithmetic) + N2 (shell-block sampling lower bound).
"""

from __future__ import annotations

import itertools
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
from ns_exploration.conjectures.l0040_sym_shor import _sym_col_index
from ns_exploration.spectral.dealias import dealias_mask
from ns_exploration.spectral.fourier_conventions import wave_number_grids

# Frozen C-R-0008 witness: first pol, consecutive existing shells through 8.
SUPPORT_CR0008 = (1, 2, 3, 4, 5, 6, 8)
BAND_123456 = (1, 2, 3, 4, 5, 6)
BAND_ONEPOL_TECHO = (1, 2, 3, 4, 5, 6, 8, 9)
POL_BRANCH = "first"


def _existing_shells(n: int, radii: tuple[int, ...]) -> tuple[int, ...]:
    by = shell_modes_by_r(n)
    return tuple(r for r in radii if r in by and by[r])


def sym_onepol_C_shor(
    n: int,
    radii: tuple[int, ...] | list[int],
    branch: str = POL_BRANCH,
) -> dict:
    """C_Sym on a single polarization branch per wavevector."""
    if branch not in ("first", "second"):
        raise ValueError(f"unknown pol branch {branch!r}")
    by = shell_modes_by_r(n)
    radii_t = _existing_shells(n, tuple(int(r) for r in radii))
    modes: list[tuple[int, int, int]] = []
    for r in radii_t:
        modes.extend(list(by[r]))
    basis: list[tuple[np.ndarray, np.ndarray, float]] = []
    for k in modes:
        lam = float(k[0] * k[0] + k[1] * k[1] + k[2] * k[2])
        pols = _pol_basis(np.array(k, float))
        if branch == "first":
            p = pols[0]
        else:
            p = pols[1] if len(pols) > 1 else pols[0]
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
        "branch": branch,
        "D": D,
        "n_sym": n_sym,
        "L_op_sym": Lop,
        "C_shor_sym": C,
    }


def _build_twopol_Lmat(n: int, radii: tuple[int, ...]) -> tuple[np.ndarray, np.ndarray]:
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


def shellblock_C_lower(n: int = 24, n_dirichlet: int = 800) -> dict:
    """
    Lower bound on the shell-pair→target block majorant for two-pol {1..6}.

    |f| ≤ Σ_{r≤s} Σ_t op[r,s,t] ‖Z_rs‖_F ‖z_t‖ with ω on the simplex.
    A sampled value > C_† proves the majorant cannot close {1..6}.
    """
    radii = BAND_123456
    Lmat, shells = _build_twopol_Lmat(n, radii)
    D = Lmat.shape[0]
    idx = {r: np.array([i for i in range(D) if shells[i] == r], dtype=int) for r in radii}

    def col_for(i: int, j: int) -> int:
        return i if i == j else _sym_col_index(i, j, D)

    def cols_rs(r: int, s: int) -> list[int]:
        I, J = idx[r], idx[s]
        cols: set[int] = set()
        if r == s:
            for a in I:
                for b in I:
                    aa, bb = int(a), int(b)
                    cols.add(col_for(aa, bb) if aa <= bb else col_for(bb, aa))
        else:
            for a in I:
                for b in J:
                    aa, bb = int(a), int(b)
                    cols.add(col_for(aa, bb) if aa <= bb else col_for(bb, aa))
        return sorted(cols)

    ops: dict[tuple[int, int, int], float] = {}
    for r, s in itertools.combinations_with_replacement(radii, 2):
        cols = cols_rs(r, s)
        for t in radii:
            rows = idx[t]
            if len(rows) == 0 or not cols:
                ops[(r, s, t)] = 0.0
                continue
            M = Lmat[np.ix_(rows, np.asarray(cols, dtype=int))]
            ops[(r, s, t)] = float(np.linalg.svd(M, compute_uv=False)[0])

    def eval_w(w: np.ndarray) -> float:
        tot = 0.0
        for r, s in itertools.combinations_with_replacement(radii, 2):
            if r == s:
                zf = w[radii.index(r)]
            else:
                zf = math.sqrt(2.0) * math.sqrt(w[radii.index(r)] * w[radii.index(s)])
            for t in radii:
                tot += ops[(r, s, t)] * zf * math.sqrt(w[radii.index(t)])
        return tot

    rng = np.random.default_rng(0)
    best = 0.0
    best_w = np.ones(len(radii)) / len(radii)
    for _ in range(n_dirichlet):
        w = rng.dirichlet(np.ones(len(radii)))
        val = eval_w(w)
        if val > best:
            best, best_w = val, w.copy()
    for i in range(len(radii)):
        w = np.zeros(len(radii))
        w[i] = 1.0
        val = eval_w(w)
        if val > best:
            best, best_w = val, w.copy()
    for i, j in itertools.combinations(range(len(radii)), 2):
        for a in np.linspace(0.0, 1.0, 11):
            w = np.zeros(len(radii))
            w[i], w[j] = a, 1.0 - a
            val = eval_w(w)
            if val > best:
                best, best_w = val, w.copy()

    C = 2.0 * math.sqrt(2.0) * best
    return {
        "D": D,
        "max_stretch_unit_lo": best,
        "C_shellblock_lo": C,
        "w_star": best_w.tolist(),
        "n_ops": len(ops),
    }


def kyfan_partial_C(n: int = 24, n_terms: int = 4) -> dict:
    """Partial Ky-Fan Σ_{i<n_terms} σ_i‖M_i‖_op; already > C_† for n_terms≥4 on {1..6}."""
    Lmat, _ = _build_twopol_Lmat(n, BAND_123456)
    D = Lmat.shape[0]
    _u, s, vt = np.linalg.svd(Lmat, full_matrices=False)

    def unsym(v: np.ndarray) -> np.ndarray:
        M = np.zeros((D, D), dtype=np.float64)
        for i in range(D):
            M[i, i] = v[i]
        c = D
        sq = math.sqrt(0.5)
        for i in range(D):
            for j in range(i + 1, D):
                M[i, j] = M[j, i] = sq * v[c]
                c += 1
        return M

    acc = 0.0
    for i in range(min(n_terms, len(s))):
        mop = float(np.linalg.norm(unsym(vt[i]), ord=2))
        acc += float(s[i]) * mop
    C = 2.0 * math.sqrt(2.0) * acc
    return {
        "n_terms": n_terms,
        "partial_sum": acc,
        "C_kyfan_partial": C,
        "sigma0": float(s[0]),
    }


@dataclass
class GalerkinBoundL0043:
    lemma_id: str = "L-0043"
    route: str = "B"
    status: str = "proved_restricted"
    evidence_level: str = "N7"
    n: int = 24
    c0007_M: float = 0.0
    C_dagger: float = 0.0
    Omega_star: float = 0.0
    support_cr0008: list[int] | None = None
    support_cr0008_C: float = 0.0
    support_cr0008_D: int = 0
    onepol_123456_C: float = 0.0
    onepol_techo_12345689_C: float = 0.0
    C_shellblock_lo: float = 0.0
    C_kyfan_partial: float = 0.0
    closes_cr0008: bool = False
    closes_123456_twopol: bool = False
    closes_c0007_all_ic: bool = False
    clay_implication: str = (
        "None. One-pol Sym subclass C-R-0008 + shell-block/Ky-Fan techos; "
        "not all-IC; not continuum; not Clay."
    )
    notes: str = ""

    def as_dict(self) -> dict:
        return asdict(self)


def lemma_l0043(n: int = 24) -> GalerkinBoundL0043:
    b7 = lemma_l0027(n=n, empirical=False)
    b6 = lemma_l0026(n=n, n_grid=21)
    Cd = b7.C_dagger
    cr = sym_onepol_C_shor(n, SUPPORT_CR0008, branch=POL_BRANCH)
    c6 = sym_onepol_C_shor(n, BAND_123456, branch=POL_BRANCH)
    c_techo = sym_onepol_C_shor(n, BAND_ONEPOL_TECHO, branch=POL_BRANCH)
    sb = shellblock_C_lower(n=n)
    kf = kyfan_partial_C(n=n, n_terms=4)
    assert cr["C_shor_sym"] <= Cd + 1e-9
    assert c_techo["C_shor_sym"] > Cd
    assert sb["C_shellblock_lo"] > Cd
    assert kf["C_kyfan_partial"] > Cd
    notes = (
        f"One-pol first: support {list(SUPPORT_CR0008)} C={cr['C_shor_sym']:.4g}≤C_† "
        f"(C-R-0008); {{1..6}} one-pol C={c6['C_shor_sym']:.4g}; "
        f"+shell9 techo C={c_techo['C_shor_sym']:.4g}>C_†; "
        f"shellblock_lo C≳{sb['C_shellblock_lo']:.4g}; "
        f"KyFan4 C≳{kf['C_kyfan_partial']:.4g}. All-IC: False."
    )
    return GalerkinBoundL0043(
        n=n,
        c0007_M=b7.c0007_M,
        C_dagger=Cd,
        Omega_star=b6.Omega_star,
        support_cr0008=list(SUPPORT_CR0008),
        support_cr0008_C=float(cr["C_shor_sym"]),
        support_cr0008_D=int(cr["D"]),
        onepol_123456_C=float(c6["C_shor_sym"]),
        onepol_techo_12345689_C=float(c_techo["C_shor_sym"]),
        C_shellblock_lo=float(sb["C_shellblock_lo"]),
        C_kyfan_partial=float(kf["C_kyfan_partial"]),
        closes_cr0008=True,
        closes_123456_twopol=False,
        closes_c0007_all_ic=False,
        notes=notes,
    )


def cr0008_record(bound: GalerkinBoundL0043) -> dict:
    return {
        "id": "C-R-0008",
        "route": "B",
        "status": "proved_restricted",
        "evidence_level": "N7",
        "lemma": "L-0043+L-0026+L-0027",
        "domain": (
            "Pseudospectral T^3, N<=24, 2/3 dealias, Fourier support in shells "
            f"|k|^2 in {list(SUPPORT_CR0008)}, single frozen polarization "
            f"(branch={POL_BRANCH} of _pol_basis) per wavevector. NOT continuum."
        ),
        "nu": 0.1,
        "energy": 0.5,
        "t_end": 0.02,
        "n_max": bound.n,
        "proposed_bound_M": bound.c0007_M,
        "shells": list(SUPPORT_CR0008),
        "pol_branch": POL_BRANCH,
        "C_shor_sym": bound.support_cr0008_C,
        "C_dagger": bound.C_dagger,
        "D": bound.support_cr0008_D,
        "Omega_star": bound.Omega_star,
        "parent_open": "C-0007",
        "statement": (
            f"For divergence-free fields on N≤{bound.n} with Fourier support only on "
            f"shells |k|²∈{list(SUPPORT_CR0008)} and a single frozen polarization "
            f"branch per k, dealiased Galerkin NS (ν=0.1, E≤0.5) satisfies "
            f"Ω(0.02)≤{bound.c0007_M}. Proof: one-pol Sym SVD Shor gives "
            f"C≤{bound.support_cr0008_C}≤C_† ⇒ L-0027; high slab via L-0026. "
            "FINITE only."
        ),
        "clay_implication": (
            "None. One-pol low-shell subclass only; not all-IC C-0007; not Clay."
        ),
        "notes": (
            "Two-pol Sym on {1..6} fails (L-0040/41). One-pol extends consecutive "
            "band through shell 8; adding shell 9 exceeds C_† (L-0043 techo)."
        ),
    }


def save_lemma_l0043(
    bound: GalerkinBoundL0043,
    path: str | Path = "conjectures/proved_restricted/L-0043.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(bound.as_dict(), indent=2), encoding="utf-8")


def save_cr0008(
    bound: GalerkinBoundL0043,
    path: str | Path = "conjectures/proved_restricted/C-R-0008.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(cr0008_record(bound), indent=2), encoding="utf-8")
