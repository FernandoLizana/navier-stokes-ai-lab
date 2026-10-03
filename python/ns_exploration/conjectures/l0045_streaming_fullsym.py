"""
L-0045 helpers: streaming fullsym matricization for physical Hermitian stretch.

Avoids materializing G (D³). Accumulates M = flatten(full_sym(G)) of shape
(D, n_sym) directly from sparse triad contributions, then C = 2√2 ‖M‖_op.

Also provides matrix-free power estimates for oversized bands (N2 only).
"""

from __future__ import annotations

import math
from collections import defaultdict

import numpy as np

from ns_exploration.conjectures.l0021_quartic_shell import _leray_vec, shell_modes_by_r
from ns_exploration.conjectures.l0040_sym_shor import _sym_col_index
from ns_exploration.experiments.sprintB_physical_tensor import (
    HermitianBasis,
    build_hermitian_basis,
)


def _mode_to_basis(b: HermitianBasis) -> dict[tuple[int, int, int], list[int]]:
    """Map wavevector -> basis indices whose ψ support includes that mode."""
    out: dict[tuple[int, int, int], list[int]] = defaultdict(list)
    for m, psi in enumerate(b.psis):
        for k, _v in psi:
            if m not in out[k]:
                out[k].append(m)
    return out


def _inner_psi_N(
    b: HermitianBasis, m: int, Nloc: dict[tuple[int, int, int], np.ndarray]
) -> float:
    acc = 0.0
    for km, vm in b.psis[m]:
        if km not in Nloc:
            continue
        lam = float(km[0] * km[0] + km[1] * km[1] + km[2] * km[2])
        acc += lam * float(np.real(np.vdot(vm, Nloc[km])))
    return acc


def streaming_fullsym_M(
    n: int,
    radii: tuple[int, ...],
    progress_every: int = 0,
    dtype: np.dtype | type = np.float64,
) -> tuple[np.ndarray, HermitianBasis]:
    """Build M = flatten(full_sym(G)) without storing G."""
    b = build_hermitian_basis(n, radii)
    by = shell_modes_by_r(n)
    out_set: set[tuple[int, int, int]] = set()
    for r in radii:
        if r in by:
            out_set.update(by[r])
    D = b.D
    n_sym = D * (D + 1) // 2
    M = np.zeros((D, n_sym), dtype=dtype)
    hit = _mode_to_basis(b)
    sq2 = dtype(math.sqrt(2.0))

    def add_A(p: int, q: int, r: int, val: float) -> None:
        """A[p,q,r] += val; update flatten columns (sym in last two of A row p)."""
        if q == r:
            M[p, q] += val
        else:
            # svec off-diag: W_qr contributes with √2, and A[p,q,r]=A[p,r,q]
            # φ(W)_p = Σ_{q<=r} ... standard: column for (q,r) gets √2 A[p,q,r]
            # when A is stored once for q<r as the symmetric value.
            lo, hi = (q, r) if q < r else (r, q)
            M[p, _sym_col_index(lo, hi, D)] += sq2 * val

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
            if not Nloc:
                done += 1
                continue
            for s, vec in list(Nloc.items()):
                Nloc[s] = -_leray_vec(np.array(s, float), vec)
            # sparse m touching Nloc
            ms: set[int] = set()
            for s in Nloc:
                ms.update(hit.get(s, ()))
            for m in ms:
                g = _inner_psi_N(b, m, Nloc)
                if g == 0.0:
                    continue
                contrib = g / 6.0
                for p, q, r in (
                    (m, a, bb),
                    (m, bb, a),
                    (a, m, bb),
                    (a, bb, m),
                    (bb, m, a),
                    (bb, a, m),
                ):
                    add_A(p, q, r, contrib)
            done += 1
            if progress_every and done % progress_every == 0:
                print(f"  pairs {done}/{n_pairs}", flush=True)
    # (q,r) and (r,q) each write √2·val into the same off-diag column → exact factor 2 vs matricize(full_sym).
    M *= 0.5
    return M, b


def _gram_largest_eig(M: np.ndarray, chunk_bytes: float = 3.0e8) -> float:
    """σmax(M) via the D×D Gram matrix G = M·Mᵀ (M is (D, n_sym), D≪n_sym).

    Avoids the huge SVD workspace (GESDD on (D, n_sym) OOMs) by forming a tiny
    D×D matrix. G is accumulated in float64 over column chunks so a float32-built
    M keeps full accuracy without ever promoting all of M to float64.
    Returns √λmax(G) = ‖M‖_op.
    """
    D = int(M.shape[0])
    ncol = int(M.shape[1])
    G = np.zeros((D, D), dtype=np.float64)
    # keep each temporary float64 column block ≲ chunk_bytes
    chunk = max(1, int(chunk_bytes // (8.0 * max(1, D))))
    for c0 in range(0, ncol, chunk):
        blk = np.asarray(M[:, c0 : c0 + chunk], dtype=np.float64)
        G += blk @ blk.T
    lam = float(np.linalg.eigvalsh(G)[-1])
    return math.sqrt(lam) if lam > 0.0 else 0.0


def _largest_singular_value(M: np.ndarray) -> float:
    """‖M‖_op. Small M: direct SVD (most accurate). Large M: Gram-matrix eig."""
    # For modest M the direct SVD is fastest and most accurate.
    if M.nbytes <= 1.5e9:
        try:
            return float(np.linalg.svd(M, compute_uv=False)[0])
        except MemoryError:
            pass
    # Large M: never materialize the SVD workspace — use the D×D Gram matrix.
    try:
        return _gram_largest_eig(M)
    except MemoryError:
        pass
    # Last-resort matrix-free power iteration (N2-grade if ever reached).
    rng = np.random.default_rng(0)
    best = 0.0
    for _ in range(4):
        v = rng.standard_normal(M.shape[1]).astype(M.dtype, copy=False)
        nv = float(np.linalg.norm(v))
        if nv == 0.0:
            continue
        v /= nv
        for _ in range(48):
            u = M @ v
            v = M.T @ u
            nv = float(np.linalg.norm(v))
            if nv == 0.0:
                break
            v /= nv
        best = max(best, float(np.linalg.norm(M @ v)))
    return best


def streaming_C_fullsym(
    n: int,
    radii: tuple[int, ...],
    dtype: np.dtype | type = np.float64,
    *,
    svd_dtype: np.dtype | type | None = None,
) -> dict:
    M, b = streaming_fullsym_M(n, radii, dtype=dtype)
    # Default: keep accumulation dtype for SVD (float32 exploration must not
    # promote M→f64 — that doubles RAM and OOMs near the dense-M cap).
    # Certificates should pass dtype=float64 (or svd_dtype=float64 on small M).
    sd = np.dtype(dtype if svd_dtype is None else svd_dtype)
    if M.dtype != sd:
        M = np.asarray(M, dtype=sd)
    Lop = _largest_singular_value(M)
    C = 2.0 * math.sqrt(2.0) * Lop
    return {
        "n": n,
        "radii": list(radii),
        "D": b.D,
        "L_op": Lop,
        "C_fullsym": C,
        "dtype": str(np.dtype(dtype)),
        "svd_dtype": str(sd),
    }
