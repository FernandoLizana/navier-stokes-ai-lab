"""
L-0049: Rank-1 / cubic power lower bound on full-dealias physical stretch.

L-0048 showed C_fullsym≈25.93≫C_† (Shor method techo). This module asks whether
the *true* cubic constant C_true = max_{‖z‖=1} 2√2 |A(z,z,z)| already exceeds
C_† (would refute low-slab / C-0007) or sits well below (gap for tighter certs).

Method (N2→N7 arithmetic on frozen sparse M):
  1. Assemble sparse CSR M = flatten(full_sym(G)) on all dealias shells.
  2. Tensor power: z ← normalize(M @ svec(z zᵀ)); C = 2√2 |z·(M svec)|.
  3. Multi-start; also FFT cross-check via stretch_inner on z_to_uhat(z).

FINITE Galerkin only. Not continuum. Not Clay.
"""

from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from ns_exploration.conjectures.l0023_cubic_dissipation import stretch_constant_of
from ns_exploration.conjectures.l0027_low_slab_cubic import lemma_l0027
from ns_exploration.conjectures.l0048_fullsym_techo import C_FULLSYM_FULL_DEALIAS
from ns_exploration.conjectures.l0048_matrixfree_fullsym import (
    all_dealias_radii,
    build_sparse_fullsym_csr,
)
from ns_exploration.experiments.sprintB_physical_tensor import z_to_uhat


def svec_zz_fast(z: np.ndarray, iu: tuple[np.ndarray, np.ndarray] | None = None) -> np.ndarray:
    """Vectorized svec(zzᵀ). Pass cached np.triu_indices(D,1) as iu."""
    D = int(z.shape[0])
    n_sym = D * (D + 1) // 2
    out = np.empty(n_sym, dtype=np.float64)
    out[:D] = z * z
    if iu is None:
        iu = np.triu_indices(D, k=1)
    out[D:] = math.sqrt(2.0) * (z[iu[0]] * z[iu[1]])
    return out


def cubic_power_C_lower(
    M,
    n_starts: int = 64,
    n_iters: int = 80,
    seed: int = 0,
) -> dict:
    """C_lo = 2√2 max |z·M svec(zzᵀ)| over power starts (‖z‖=1)."""
    D = M.shape[0]
    iu = np.triu_indices(D, k=1)
    rng = np.random.default_rng(seed)
    best = 0.0
    best_z = None
    history: list[float] = []
    for s in range(n_starts):
        z = rng.standard_normal(D)
        z /= np.linalg.norm(z) + 1e-30
        val = 0.0
        for _ in range(n_iters):
            w = svec_zz_fast(z, iu)
            g = M @ w  # A(z,z,·)
            val = abs(float(z @ g))
            gn = float(np.linalg.norm(g))
            if gn < 1e-30:
                break
            z = g / gn
        history.append(val)
        if val > best:
            best = val
            best_z = z.copy()
        if (s + 1) % 8 == 0:
            print(
                f"  power start {s+1}/{n_starts}: best_raw={best:.6g} "
                f"C>={2*math.sqrt(2)*best:.6g}",
                flush=True,
            )
    C = 2.0 * math.sqrt(2.0) * best
    return {
        "L_raw_max": best,
        "C_rank1_lower": C,
        "n_starts": n_starts,
        "n_iters": n_iters,
        "best_z": best_z,
        "history_raw": history,
    }


def fft_C_of_z(z: np.ndarray, basis) -> float:
    uh = z_to_uhat(z, basis)
    return float(stretch_constant_of(uh)[0])


@dataclass
class GalerkinBoundL0049:
    lemma_id: str = "L-0049"
    route: str = "B"
    status: str = "exploring"
    evidence_level: str = "N2"
    n: int = 24
    C_dagger: float = 0.0
    C_fullsym_full: float = C_FULLSYM_FULL_DEALIAS
    C_rank1_lower: float = 0.0
    C_fft_at_best_z: float = 0.0
    D: int = 0
    nnz: int = 0
    gap_fullsym_over_rank1: float = 0.0
    rank1_exceeds_Cdagger: bool = False
    closes_c0007_all_ic: bool = False
    refutes_c0007: bool = False
    clay_implication: str = (
        "None. Rank-1 lower bound on physical cubic; finite Galerkin only; not Clay."
    )
    notes: str = ""

    def as_dict(self) -> dict:
        return asdict(self)


def run_l0049(
    n: int = 24,
    n_starts: int = 64,
    n_iters: int = 80,
    M_cache: str | Path = "reports/l0049_full_M.npz",
    rebuild_M: bool = False,
) -> GalerkinBoundL0049:
    from scipy import sparse

    Cd = lemma_l0027(n=n, empirical=False).C_dagger
    cache = Path(M_cache)
    basis = None
    if cache.exists() and not rebuild_M:
        print(f"load M from {cache}", flush=True)
        data = np.load(cache, allow_pickle=False)
        M = sparse.csr_matrix(
            (data["data"], data["indices"], data["indptr"]),
            shape=tuple(int(x) for x in data["shape"]),
        )
        D = int(data["D"])
        nnz = int(M.nnz)
    else:
        print("assemble sparse M (full dealias)...", flush=True)
        M, basis = build_sparse_fullsym_csr(
            n, all_dealias_radii(n), progress_every=500_000
        )
        D = M.shape[0]
        nnz = int(M.nnz)
        cache.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(
            cache,
            data=M.data,
            indices=M.indices,
            indptr=M.indptr,
            shape=np.array(M.shape, dtype=np.int64),
            D=np.int64(D),
        )
        print(f"wrote {cache} nnz={nnz}", flush=True)

    print(f"cubic power D={D} nnz={nnz}...", flush=True)
    pw = cubic_power_C_lower(M, n_starts=n_starts, n_iters=n_iters)
    C_r1 = float(pw["C_rank1_lower"])
    C_fft = 0.0
    if pw["best_z"] is not None:
        if basis is None:
            from ns_exploration.experiments.sprintB_physical_tensor import (
                build_hermitian_basis,
            )

            basis = build_hermitian_basis(n, all_dealias_radii(n))
        C_fft = fft_C_of_z(pw["best_z"], basis)
        print(f"FFT C at best z: {C_fft:.6g}", flush=True)

    gap = (C_FULLSYM_FULL_DEALIAS / C_r1) if C_r1 > 0 else float("inf")
    exceeds = C_r1 > Cd + 1e-9
    notes = (
        f"Full-dealias rank1: C_lo≈{C_r1:.4f} (FFT@{C_fft:.4f}); "
        f"C_fullsym≈{C_FULLSYM_FULL_DEALIAS:.4f}; C_†≈{Cd:.4f}; "
        f"gap_fullsym/rank1≈{gap:.2f}; exceeds_C†={exceeds}."
    )
    return GalerkinBoundL0049(
        n=n,
        C_dagger=Cd,
        C_fullsym_full=C_FULLSYM_FULL_DEALIAS,
        C_rank1_lower=C_r1,
        C_fft_at_best_z=C_fft,
        D=D,
        nnz=nnz,
        gap_fullsym_over_rank1=float(gap),
        rank1_exceeds_Cdagger=exceeds,
        closes_c0007_all_ic=False,
        refutes_c0007=exceeds,  # only if lower bound on C_true > C_†
        notes=notes,
        evidence_level="N2",
        status="techo" if exceeds else "exploring",
    )


def save_lemma_l0049(
    bound: GalerkinBoundL0049,
    path: str | Path = "conjectures/active/L-0049.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    d = bound.as_dict()
    path.write_text(json.dumps(d, indent=2), encoding="utf-8")
