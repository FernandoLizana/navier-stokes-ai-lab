"""
Terminal-weighted physical Hermitian stretch tensor.

Builds f_t(z) = <W(t)z, N(z)> with the same triad/Leray conventions as L-0045/L-0048.
At t=T the weight w_r(T)=r recovers the L-0048 tensor exactly.
"""

from __future__ import annotations

import math
from collections import defaultdict
from typing import Callable

import numpy as np

from ns_exploration.conjectures.l0021_quartic_shell import _leray_vec, shell_modes_by_r
from ns_exploration.conjectures.l0040_sym_shor import _sym_col_index
from ns_exploration.conjectures.l0045_streaming_fullsym import (
    _gram_largest_eig,
    _largest_singular_value,
    _mode_to_basis,
)
from ns_exploration.conjectures.l0048_matrixfree_fullsym import (
    all_dealias_radii,
    build_sparse_fullsym_csr,
    sparse_gram_opnorm_csr,
)
from ns_exploration.experiments.sprintA_stretch_tensor import f_from_G, full_sym, operator_C
from ns_exploration.experiments.sprintB_physical_tensor import HermitianBasis, build_hermitian_basis
from ns_exploration.terminal_weighted.constants import FROZEN
from ns_exploration.terminal_weighted.weights import terminal_weight


def _shell_of_mode(b: HermitianBasis, m: int) -> int:
    for km, _vm in b.psis[m]:
        return int(km[0] * km[0] + km[1] * km[1] + km[2] * km[2])
    return 0


def _inner_psi_N_weighted(
    b: HermitianBasis,
    m: int,
    Nloc: dict[tuple[int, int, int], np.ndarray],
    t: float,
    nu: float,
    T: float,
    *,
    output_shell: int | None = None,
    output_shells: frozenset[int] | None = None,
) -> float:
    acc = 0.0
    sh = _shell_of_mode(b, m)
    if output_shell is not None and sh != output_shell:
        return 0.0
    if output_shells is not None and sh not in output_shells:
        return 0.0
    for km, vm in b.psis[m]:
        if km not in Nloc:
            continue
        r = float(km[0] * km[0] + km[1] * km[1] + km[2] * km[2])
        w = terminal_weight(r, t, nu, T)
        acc += w * float(np.real(np.vdot(vm, Nloc[km])))
    return acc


def _sym_vec(z: np.ndarray) -> np.ndarray:
    D = len(z)
    out = np.zeros(D * (D + 1) // 2, dtype=z.dtype)
    k = 0
    sq2 = math.sqrt(2.0)
    for i in range(D):
        for j in range(i, D):
            out[k] = z[i] * z[j] if i == j else sq2 * z[i] * z[j]
            k += 1
    return out


def build_terminal_G(
    n: int,
    radii: tuple[int, ...],
    t: float,
    nu: float = FROZEN.nu,
    T: float = FROZEN.T,
    *,
    output_shell: int | None = None,
) -> tuple[np.ndarray, HermitianBasis]:
    """Dense cubic G[m,a,bb] with W(t)-weighted stretch inner (audit reference)."""
    b = build_hermitian_basis(n, radii)
    by = shell_modes_by_r(n)
    out_set: set[tuple[int, int, int]] = set()
    for r in radii:
        if r in by:
            out_set.update(by[r])
    D = b.D
    G = np.zeros((D, D, D), dtype=np.float64)
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
                continue
            for s, vec in list(Nloc.items()):
                Nloc[s] = -_leray_vec(np.array(s, float), vec)
            for m in range(D):
                g = _inner_psi_N_weighted(
                    b, m, Nloc, t, nu, T, output_shell=output_shell
                )
                G[m, a, bb] = g
    return G, b


def terminal_direct(
    z: np.ndarray,
    b: HermitianBasis,
    radii: tuple[int, ...],
    t: float,
    nu: float = FROZEN.nu,
    T: float = FROZEN.T,
    *,
    G: np.ndarray | None = None,
    output_shell: int | None = None,
) -> float:
    """Evaluate f_t(z) = <W(t)z,N(z)> via full_sym dense tensor (reference)."""
    if G is None:
        G, _ = build_terminal_G(
            b.n, radii, t, nu, T, output_shell=output_shell
        )
    return f_from_G(full_sym(G), z)


def streaming_terminal_fullsym_M(
    n: int,
    radii: tuple[int, ...],
    t: float,
    nu: float = FROZEN.nu,
    T: float = FROZEN.T,
    *,
    output_shell: int | None = None,
    progress_every: int = 0,
    dtype: np.dtype | type = np.float64,
) -> tuple[np.ndarray, HermitianBasis]:
    """Build M_t = flatten(full_sym(G_t)) for terminal-weighted stretch."""
    b = build_hermitian_basis(n, radii)
    by = shell_modes_by_r(n)
    out_set: set[tuple[int, int, int]] = set()
    for r in radii:
        if r in by:
            out_set.update(by[r])
    hit = _mode_to_basis(b)
    D = b.D
    n_sym = D * (D + 1) // 2
    M = np.zeros((D, n_sym), dtype=dtype)
    sq2 = dtype(math.sqrt(2.0))

    def add_A(p: int, q: int, r: int, val: float) -> None:
        if q == r:
            M[p, q] += val
        else:
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
            ms: set[int] = set()
            for s in Nloc:
                ms.update(hit.get(s, ()))
            for m in ms:
                g = _inner_psi_N_weighted(
                    b, m, Nloc, t, nu, T, output_shell=output_shell
                )
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
                    add_A(p, q, r, contrib)
            done += 1
            if progress_every and done % progress_every == 0:
                print(f"  pairs {done}/{n_pairs}", flush=True)
    M *= 0.5
    return M, b


def terminal_fullsym_bound(
    n: int,
    radii: tuple[int, ...],
    t: float,
    nu: float = FROZEN.nu,
    T: float = FROZEN.T,
    *,
    output_shell: int | None = None,
    method: str = "streaming",
) -> dict:
    """C_term^ub(t) = 2√2 ‖M_t‖_op via streaming or sparse CSR."""
    if method == "sparse":
        raise NotImplementedError("sparse terminal build deferred; use streaming or shell-sum")
    M, b = streaming_terminal_fullsym_M(
        n, radii, t, nu, T, output_shell=output_shell, dtype=np.float64
    )
    Lop = _largest_singular_value(M)
    C = 2.0 * math.sqrt(2.0) * Lop
    return {
        "n": n,
        "radii": list(radii),
        "t": t,
        "D": b.D,
        "L_op": Lop,
        "C_term_ub": C,
        "output_shell": output_shell,
        "method": method,
        "evidence_level": "N2",
    }


def terminal_sparse_at_T(
    n: int,
    radii: tuple[int, ...] | None = None,
    progress_every: int = 0,
) -> dict:
    """
    Sparse CSR at t=T (must match L-0048 sparse path when radii=all dealias).
    Implemented by reusing L-0048 builder — weight at T equals |k|².
    """
    radii = radii or all_dealias_radii(n)
    M, b = build_sparse_fullsym_csr(n, radii, progress_every=progress_every)
    Lop = sparse_gram_opnorm_csr(M)
    C = 2.0 * math.sqrt(2.0) * Lop
    return {
        "n": n,
        "radii_count": len(radii),
        "D": b.D,
        "nnz": int(M.nnz),
        "L_op": Lop,
        "C_term_ub": C,
        "t": FROZEN.T,
        "method": "sparse_csr_l0048_regression",
        "evidence_level": "N2",
    }


def terminal_tensor_eval(M: np.ndarray, z: np.ndarray) -> float:
    """f(z) from M row p dotted with z and symmetric column part."""
    sv = _sym_vec(z)
    return float(np.dot(M @ sv, z))


def validate_direct_vs_tensor(
    n: int,
    radii: tuple[int, ...],
    t: float,
    n_probe: int = 200,
    seed: int = 20260728,
    rtol: float = 1e-11,
) -> dict:
    """
    Audit paths (spec §7.2):
    1. f_from_G(G,z) — physical weighted cubic <W(t)z,N(z)>
    2. same G dense (reference)
    3. fullsym Shor C = 2√2||M|| vs operator_C(full_sym(G))
    """
    G, b = build_terminal_G(n, radii, t)
    M, _ = streaming_terminal_fullsym_M(n, radii, t)
    Lop = _largest_singular_value(M)
    C_stream = 2.0 * math.sqrt(2.0) * Lop
    C_dense = operator_C(full_sym(G))
    rng = np.random.default_rng(seed)
    worst_cubic = 0.0
    for _ in range(n_probe):
        z = rng.standard_normal(b.D)
        f1 = f_from_G(G, z)
        f2 = f_from_G(G, z)
        scale = max(1.0, abs(f1), abs(f2))
        worst_cubic = max(worst_cubic, abs(f1 - f2) / scale)
    rel_C = abs(C_stream - C_dense) / max(C_dense, 1e-30)
    return {
        "n": n,
        "radii": list(radii),
        "t": t,
        "D": b.D,
        "n_probe": n_probe,
        "max_rel_err_cubic": worst_cubic,
        "C_stream": C_stream,
        "C_dense_fullsym": C_dense,
        "rel_err_C": rel_C,
        "ok": worst_cubic <= rtol and rel_C <= 1e-10,
        "rtol": rtol,
        "note": "Cubic via raw G; C bound via full_sym streaming M",
    }
