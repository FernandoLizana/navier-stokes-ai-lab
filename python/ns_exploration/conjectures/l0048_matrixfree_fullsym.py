"""
L-0048 helpers: sparse / matrix-free physical fullsym op-norm (all-IC attack).

Dense M at full dealias (D≈6748) is ~0.6 TB. Empirically M is highly sparse
(~0.5–1% on small bands). This module:

1. Streams triad contributions into a COO / dict-of-columns sparse M.
2. Computes C = 2√2 ‖M‖_op via Gram eig on the sparse columns (or power).
3. Provides matrix-free matvecs for N2 lower bounds when sparse assembly OOMs.

FINITE Galerkin only. Not continuum. Not Clay. Evidence ≤ N2 until certified.
"""

from __future__ import annotations

import math

import numpy as np

from ns_exploration.conjectures.l0021_quartic_shell import _leray_vec, shell_modes_by_r
from ns_exploration.conjectures.l0040_sym_shor import _sym_col_index
from ns_exploration.conjectures.l0045_streaming_fullsym import (
    _inner_psi_N,
    _mode_to_basis,
    streaming_C_fullsym,
)
from ns_exploration.experiments.sprintB_physical_tensor import (
    HermitianBasis,
    build_hermitian_basis,
)


def all_dealias_radii(n: int = 24) -> tuple[int, ...]:
    return tuple(sorted(shell_modes_by_r(n)))


def _iter_fullsym_updates(
    b: HermitianBasis,
    out_set: set[tuple[int, int, int]],
    hit: dict,
    progress_every: int = 0,
):
    """Yield (p, col, delta) raw updates BEFORE the global factor-1/2."""
    D = b.D
    sq2 = math.sqrt(2.0)
    n_pairs = D * D
    done = 0
    for a in range(D):
        for bb in range(D):
            Nloc: dict[tuple[int, int, int], np.ndarray] = {}
            for kp, va in b.psis[a]:
                for kq, vb in b.psis[bb]:
                    s = (kp[0] + kq[0], kp[1] + kq[1], kp[2] + kq[2])
                    if s not in out_set:
                        continue
                    coef = 1j * np.dot(va, np.array(kq, float))
                    if s not in Nloc:
                        Nloc[s] = np.zeros(3, dtype=np.complex128)
                    Nloc[s] += coef * vb
            if Nloc:
                for s, vec in list(Nloc.items()):
                    Nloc[s] = -_leray_vec(np.array(s, float), vec)
                ms: set[int] = set()
                for s in Nloc:
                    ms.update(hit.get(s, ()))
                for m in ms:
                    g = _inner_psi_N(b, m, Nloc)
                    if g == 0.0:
                        continue
                    contrib = g / 6.0
                    perms = {
                        (m, a, bb),
                        (m, bb, a),
                        (a, m, bb),
                        (a, bb, m),
                        (bb, m, a),
                        (bb, a, m),
                    }
                    for p, q, r in perms:
                        if q == r:
                            yield p, q, contrib
                        else:
                            lo, hi = (q, r) if q < r else (r, q)
                            yield p, _sym_col_index(lo, hi, D), sq2 * contrib
            done += 1
            if progress_every and done % progress_every == 0:
                print(f"  pairs {done}/{n_pairs}", flush=True)


def build_sparse_fullsym_csr(
    n: int,
    radii: tuple[int, ...],
    progress_every: int = 0,
):
    """Assemble M as scipy.sparse.csr_matrix (D × n_sym), factor-1/2 baked in."""
    from scipy.sparse import coo_matrix

    b = build_hermitian_basis(n, radii)
    by = shell_modes_by_r(n)
    out_set: set[tuple[int, int, int]] = set()
    for r in radii:
        if r in by:
            out_set.update(by[r])
    hit = _mode_to_basis(b)
    D = b.D
    n_sym = D * (D + 1) // 2

    # geometric buffers (int32 indices — D,n_sym < 2^31)
    cap = 1 << 20
    rows = np.empty(cap, dtype=np.int32)
    cols_i = np.empty(cap, dtype=np.int32)
    data = np.empty(cap, dtype=np.float64)
    nent = 0

    def _grow(need: int) -> None:
        nonlocal cap, rows, cols_i, data
        new_cap = max(cap * 2, need)
        nr = np.empty(new_cap, dtype=np.int32)
        nc = np.empty(new_cap, dtype=np.int32)
        nd = np.empty(new_cap, dtype=np.float64)
        nr[:nent] = rows[:nent]
        nc[:nent] = cols_i[:nent]
        nd[:nent] = data[:nent]
        rows, cols_i, data, cap = nr, nc, nd, new_cap

    for p, col, delta in _iter_fullsym_updates(b, out_set, hit, progress_every):
        if nent >= cap:
            _grow(cap * 2)
        rows[nent] = p
        cols_i[nent] = col
        data[nent] = 0.5 * delta
        nent += 1

    print(f"  raw updates={nent} cap={cap}", flush=True)
    M = coo_matrix(
        (data[:nent], (rows[:nent], cols_i[:nent])),
        shape=(D, n_sym),
        dtype=np.float64,
    )
    # free buffers before CSR densification of duplicates
    del rows, cols_i, data
    M = M.tocsr()
    M.sum_duplicates()
    M.eliminate_zeros()
    return M, b


def sparse_gram_opnorm_csr(M) -> float:
    """‖M‖_op via dense D×D Gram from CSR M."""
    # M is (D, n_sym) with D ≪ n_sym → G = M @ M.T is D×D
    G = (M @ M.T).toarray()
    lam = float(np.linalg.eigvalsh(G)[-1])
    return math.sqrt(lam) if lam > 0.0 else 0.0


def sparse_C_fullsym(
    n: int,
    radii: tuple[int, ...],
    progress_every: int = 0,
) -> dict:
    M, b = build_sparse_fullsym_csr(n, radii, progress_every=progress_every)
    nnz = int(M.nnz)
    Lop = sparse_gram_opnorm_csr(M)
    C = 2.0 * math.sqrt(2.0) * Lop
    return {
        "n": n,
        "radii": list(radii),
        "D": b.D,
        "nnz": nnz,
        "L_op": Lop,
        "C_fullsym": C,
        "method": "sparse_csr_gram",
    }


def mf_apply_M(
    b: HermitianBasis,
    out_set: set[tuple[int, int, int]],
    hit: dict,
    v: np.ndarray,
    progress_every: int = 0,
) -> np.ndarray:
    """y = M @ v without storing M (one triad pass)."""
    y = np.zeros(b.D, dtype=np.float64)
    for p, col, delta in _iter_fullsym_updates(b, out_set, hit, progress_every):
        y[p] += 0.5 * delta * float(v[col])
    return y


def mf_apply_MT(
    b: HermitianBasis,
    out_set: set[tuple[int, int, int]],
    hit: dict,
    u: np.ndarray,
    progress_every: int = 0,
) -> np.ndarray:
    """w = M.T @ u without storing M (one triad pass)."""
    n_sym = b.D * (b.D + 1) // 2
    w = np.zeros(n_sym, dtype=np.float64)
    for p, col, delta in _iter_fullsym_updates(b, out_set, hit, progress_every):
        w[col] += 0.5 * delta * float(u[p])
    return w


def mf_power_C_lower(
    n: int,
    radii: tuple[int, ...],
    n_iter: int = 8,
    n_starts: int = 2,
    seed: int = 0,
    progress_every: int = 0,
) -> dict:
    """N2 lower bound on C_fullsym via power iteration on MMᵀ (matrix-free)."""
    b = build_hermitian_basis(n, radii)
    by = shell_modes_by_r(n)
    out_set: set[tuple[int, int, int]] = set()
    for r in radii:
        if r in by:
            out_set.update(by[r])
    hit = _mode_to_basis(b)
    rng = np.random.default_rng(seed)
    best = 0.0
    for s in range(n_starts):
        u = rng.standard_normal(b.D)
        u /= np.linalg.norm(u) + 1e-30
        for it in range(n_iter):
            w = mf_apply_MT(b, out_set, hit, u, progress_every=progress_every)
            y = mf_apply_M(b, out_set, hit, w, progress_every=progress_every)
            ny = float(np.linalg.norm(y))
            if ny == 0.0:
                break
            # Rayleigh for σ: ‖M^T u‖ after first half, or √(u·MM^T u)=‖M^T u‖
            sig = float(np.linalg.norm(w))
            best = max(best, sig)
            u = y / ny
            print(f"  start {s} iter {it}: sigma>={sig:.6g} C>={2*math.sqrt(2)*sig:.6g}", flush=True)
        # final sigma from last u
        w = mf_apply_MT(b, out_set, hit, u, progress_every=progress_every)
        best = max(best, float(np.linalg.norm(w)))
    C = 2.0 * math.sqrt(2.0) * best
    return {
        "n": n,
        "radii": list(radii),
        "D": b.D,
        "L_op_lower": best,
        "C_fullsym_lower": C,
        "n_iter": n_iter,
        "n_starts": n_starts,
        "method": "mf_power_MMT",
        "evidence_level": "N2",
    }


def validate_sparse_vs_dense(n: int = 24, radii: tuple[int, ...] = tuple(range(1, 7))) -> dict:
    dense = streaming_C_fullsym(n, radii, dtype=np.float64)
    sparse = sparse_C_fullsym(n, radii)
    rel = abs(sparse["C_fullsym"] - dense["C_fullsym"]) / max(dense["C_fullsym"], 1e-30)
    return {
        "dense_C": dense["C_fullsym"],
        "sparse_C": sparse["C_fullsym"],
        "rel_err": rel,
        "nnz": sparse["nnz"],
        "D": sparse["D"],
        "ok": rel < 1e-9,
    }
