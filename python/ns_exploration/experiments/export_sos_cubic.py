"""
Export physical cubic stretch polynomial for Julia TSSOS / CS-TSSOS.

Writes a JSON + sparse COO of the fully symmetrized tensor A = full_sym(G)
on a Hermitian (c,s) basis, plus the target β = C_†/(2√2).

SOS asks for max |A(z,z,z)| s.t. ‖z‖²=1; if ub < β then C=2√2·ub < C_†.

FINITE Galerkin only. Not continuum. Not Clay.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

from ns_exploration.conjectures.l0047_greedy_fullsym import SUPPORT_CR0012
from ns_exploration.conjectures.l0043_onepol_shellblock import (
    SUPPORT_CR0008,
    sym_onepol_C_shor,
)
from ns_exploration.conjectures.l0027_low_slab_cubic import lemma_l0027
from ns_exploration.conjectures.l0048_matrixfree_fullsym import (
    all_dealias_radii,
    build_sparse_fullsym_csr,
)
from ns_exploration.experiments.sprintA_stretch_tensor import (
    build_onepol_G,
    full_sym,
    matricize,
    operator_C,
)
from ns_exploration.experiments.sprintB_physical_tensor import build_physical_G


def export_dense_band(
    n: int,
    radii: tuple[int, ...],
    out_dir: Path,
) -> dict:
    """Dense A for small bands (e.g. {1..6})."""
    Cd = lemma_l0027(n=n, empirical=False).C_dagger
    beta = Cd / (2.0 * math.sqrt(2.0))
    G, b = build_physical_G(n, radii)
    A = full_sym(G)
    M = matricize(A)
    Lop = float(np.linalg.svd(M, compute_uv=False)[0])
    C_fullsym = 2.0 * math.sqrt(2.0) * Lop

    out_dir.mkdir(parents=True, exist_ok=True)
    # sparse COO of A (only nonzero)
    nz = np.argwhere(np.abs(A) > 1e-15)
    data = A[nz[:, 0], nz[:, 1], nz[:, 2]]
    np.savez_compressed(
        out_dir / "A_coo.npz",
        i=nz[:, 0].astype(np.int32),
        j=nz[:, 1].astype(np.int32),
        k=nz[:, 2].astype(np.int32),
        v=data.astype(np.float64),
        D=np.int64(b.D),
    )
    meta = {
        "n": n,
        "radii": list(radii),
        "D": b.D,
        "nnz_A": int(data.size),
        "C_dagger": Cd,
        "beta": beta,
        "C_fullsym": C_fullsym,
        "format": "dense_fullsym_A_coo",
        "objective": "max A(z,z,z) s.t. ||z||^2 = 1; need ub < beta",
        "clay_implication": "None. Finite Galerkin only.",
    }
    (out_dir / "meta.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    return meta


def export_from_sparse_M(
    n: int,
    radii: tuple[int, ...],
    out_dir: Path,
    progress_every: int = 0,
) -> dict:
    """Export sparse M (flatten fullsym) for large bands — Julia rebuilds cubic via M."""
    Cd = lemma_l0027(n=n, empirical=False).C_dagger
    beta = Cd / (2.0 * math.sqrt(2.0))
    M, b = build_sparse_fullsym_csr(n, radii, progress_every=progress_every)
    # Gram op-norm for reference
    Ggram = (M @ M.T).toarray()
    Lop = math.sqrt(max(float(np.linalg.eigvalsh(Ggram)[-1]), 0.0))
    C_fullsym = 2.0 * math.sqrt(2.0) * Lop
    out_dir.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        out_dir / "M_csr.npz",
        data=M.data,
        indices=M.indices,
        indptr=M.indptr,
        shape=np.array(M.shape, dtype=np.int64),
        D=np.int64(b.D),
    )
    meta = {
        "n": n,
        "radii": list(radii),
        "D": int(b.D),
        "nnz_M": int(M.nnz),
        "C_dagger": Cd,
        "beta": beta,
        "C_fullsym": C_fullsym,
        "format": "sparse_fullsym_M_csr",
        "objective": "max z'*(M*svec(zz')) s.t. ||z||^2=1; need 2*sqrt(2)*ub < C_dagger",
        "clay_implication": "None. Finite Galerkin only.",
    }
    (out_dir / "meta.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    return meta


def export_onepol_band(
    n: int,
    radii: tuple[int, ...],
    out_dir: Path,
    branch: str = "first",
) -> dict:
    """One-pol stretch cubic for structured subclass SOS (e.g. C-R-0008, D~92)."""
    Cd = lemma_l0027(n=n, empirical=False).C_dagger
    beta = Cd / (2.0 * math.sqrt(2.0))
    G, b = build_onepol_G(n, radii, branch=branch)
    A = full_sym(G)
    C_fullsym = operator_C(A)
    ref = sym_onepol_C_shor(n, radii, branch=branch)
    C_shor_sym = float(ref["C_shor_sym"])

    out_dir.mkdir(parents=True, exist_ok=True)
    nz = np.argwhere(np.abs(A) > 1e-15)
    data = A[nz[:, 0], nz[:, 1], nz[:, 2]]
    np.savez_compressed(
        out_dir / "A_coo.npz",
        i=nz[:, 0].astype(np.int32),
        j=nz[:, 1].astype(np.int32),
        k=nz[:, 2].astype(np.int32),
        v=data.astype(np.float64),
        D=np.int64(b.D),
    )
    meta = {
        "n": n,
        "radii": list(ref["radii"]),
        "branch": branch,
        "basis": "one_pol_stretch",
        "D": int(b.D),
        "nnz_A": int(data.size),
        "C_dagger": Cd,
        "beta": beta,
        "C_fullsym": C_fullsym,
        "C_shor_sym": C_shor_sym,
        "format": "onepol_fullsym_A_coo",
        "objective": "max A(z,z,z) s.t. ||z||^2=1; subclass one-pol only",
        "clay_implication": "None. Finite Galerkin subclass only.",
    }
    (out_dir / "meta.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    return meta


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--n", type=int, default=24)
    p.add_argument("--radii", type=int, nargs="*", default=None)
    p.add_argument("--upto", type=int, default=0, help="consecutive shells 1..upto")
    p.add_argument("--full-dealias", action="store_true")
    p.add_argument("--out", type=str, default="tools/sos_julia/data/band")
    p.add_argument("--dense", action="store_true", help="force dense Hermitian A export")
    p.add_argument("--onepol", action="store_true", help="one-pol subclass (l0043)")
    p.add_argument("--branch", type=str, default="first", choices=("first", "second"))
    p.add_argument("--cr0008", action="store_true", help="shorthand radii=C-R-0008 support")
    p.add_argument("--cr0012", action="store_true", help="shorthand radii=C-R-0012 greedy 41 shells")
    args = p.parse_args()
    if args.full_dealias:
        radii = all_dealias_radii(args.n)
    elif args.cr0012:
        radii = SUPPORT_CR0012
    elif args.cr0008:
        radii = SUPPORT_CR0008
    elif args.upto:
        radii = tuple(range(1, args.upto + 1))
    elif args.radii:
        radii = tuple(args.radii)
    else:
        radii = tuple(range(1, 7))
    out = Path(args.out)
    if args.onepol:
        meta = export_onepol_band(args.n, radii, out, branch=args.branch)
    elif args.dense or (not args.full_dealias and len(radii) <= 6):
        meta = export_dense_band(args.n, radii, out)
    else:
        meta = export_from_sparse_M(args.n, radii, out, progress_every=200_000)
    print(json.dumps(meta, indent=2), flush=True)


if __name__ == "__main__":
    main()
