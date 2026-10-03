"""
L-0037: Shell-difference stretch + bi-shell Shor techo (C_† / C-0007).

Energy conservation of the nonlinear term gives Σ_k T_k = 0 with
T_k = Re(û̄_k · N̂_k). The stretch is
  stretch = Σ_k |k|² T_k.
Hence for any reference λ,
  stretch = Σ_k (|k|² - λ) T_k.
In particular:
  (mono-radial) support on a single shell r ⇒ stretch = r Σ T_k = 0.
This is the exact shell-difference cancellation (equal-radius triads drop out).

On a bi-shell support {r1,r2}, stretch = (r1-r2) T_{r1}. Restricting the
Ω-weighted polarization Shor/matricization of L-0036 to that subspace yields
a majorant C_Shor(r1,r2). On N=24 dealias:
  - near pairs can meet C_† (e.g. (1,2): C_Shor ≲ 6.7 ≤ C_†);
  - wide pairs do not (sample max C_Shor ≳ 32 > C_†, e.g. (2,134)).
So shell-diff helps structurally (mono = 0; near bi-shell OK) but bi-shell
Shor still fails uniformly — not enough to close all-IC C-0007 / C_†.

Still missing: SOS / rank-constrained control across many shells.

FINITE dealias Galerkin only. Not continuum. Not Clay.
Evidence: N7 (mono identity) + N2 (bi-shell Shor power iteration).
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
from ns_exploration.spectral.dealias import dealias_mask
from ns_exploration.spectral.fourier_conventions import wave_number_grids


def mono_radial_stretch_identity() -> str:
    """Document the N7 identity; stretch vanishes on mono-radial fields."""
    return (
        "If supp(û) ⊂ {|k|² = r}, then stretch = Σ |k|² T_k = r Σ T_k = 0 "
        "because Σ T_k = ⟨u, N⟩ = 0 for the Leray–Galerkin nonlinearity."
    )


def bishell_shor_C_lower(
    n: int,
    r1: int,
    r2: int,
    iters: int = 20,
    seed: int = 0,
) -> dict:
    """Lower-bound C_Shor on the Ω-weighted bi-shell polarization subspace."""
    by = shell_modes_by_r(n)
    if r1 not in by or r2 not in by:
        raise KeyError(f"shells {r1},{r2} not in dealias mask n={n}")
    modes = list(by[r1]) if r1 == r2 else list(by[r1]) + list(by[r2])
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
    I_a = np.asarray(I, np.int32)
    J_a = np.asarray(J, np.int32)
    SI_a = np.asarray(SI, np.int32)
    VX_a = np.asarray(VX, float)
    VY_a = np.asarray(VY, float)
    VZ_a = np.asarray(VZ, float)

    by_s: dict[int, list[tuple[int, np.ndarray, float]]] = defaultdict(list)
    for a, (k, p, lam) in enumerate(basis):
        by_s[out[tuple(int(x) for x in k)]].append((a, p, lam))

    def apply_L(Z: np.ndarray) -> np.ndarray:
        field = np.zeros((n_out, 3), dtype=np.float64)
        np.add.at(field[:, 0], SI_a, Z[I_a, J_a] * VX_a)
        np.add.at(field[:, 1], SI_a, Z[I_a, J_a] * VY_a)
        np.add.at(field[:, 2], SI_a, Z[I_a, J_a] * VZ_a)
        outv = np.zeros(D, dtype=np.float64)
        for si, lst in by_s.items():
            Ns = field[si]
            for a, p, lam in lst:
                outv[a] += math.sqrt(lam) * float(np.dot(p, Ns))
        return outv

    rng = np.random.default_rng(seed)
    Z = rng.normal(size=(D, D))
    Z /= np.linalg.norm(Z) + 1e-30
    est = 0.0
    for _ in range(int(iters)):
        v = apply_L(Z)
        U = np.zeros((n_out, 3), dtype=np.float64)
        for si, lst in by_s.items():
            for a, p, lam in lst:
                U[si] += v[a] * math.sqrt(lam) * p
        Z2 = np.zeros((D, D), dtype=np.float64)
        np.add.at(Z2, (I_a, J_a), U[SI_a, 0] * VX_a + U[SI_a, 1] * VY_a + U[SI_a, 2] * VZ_a)
        Z = Z2 / (float(np.linalg.norm(Z2)) + 1e-30)
        est = float(np.linalg.norm(apply_L(Z)))
    C_lo = 2.0 * math.sqrt(2.0) * est
    return {
        "n": n,
        "r1": int(r1),
        "r2": int(r2),
        "D": D,
        "n_pairs": int(len(I_a)),
        "L_op_lower": est,
        "C_shor_lower": C_lo,
        "iters": int(iters),
    }


def scan_bishell_shor(
    n: int = 24,
    iters: int = 16,
    max_pairs: int | None = None,
) -> list[dict]:
    """Scan a structured sample of shell pairs; return rows sorted by C desc."""
    by = shell_modes_by_r(n)
    radii = sorted(by.keys())
    # Extreme gaps first so truncation still sees wide pairs.
    pairs: list[tuple[int, int]] = []
    for r1 in radii[:5]:
        for r2 in radii[-5:]:
            pairs.append((r1, r2))
    for i, r1 in enumerate(radii[:20]):
        for r2 in radii[i + 1 : min(i + 15, len(radii))]:
            pairs.append((r1, r2))
    pairs = list(dict.fromkeys(pairs))
    if max_pairs is not None:
        pairs = pairs[: int(max_pairs)]
    rows = [bishell_shor_C_lower(n, r1, r2, iters=iters) for r1, r2 in pairs]
    rows.sort(key=lambda d: float(d["C_shor_lower"]), reverse=True)
    return rows


@dataclass
class GalerkinBoundL0037:
    lemma_id: str = "L-0037"
    route: str = "B"
    status: str = "proved_restricted"
    evidence_level: str = "N7"
    n: int = 24
    c0007_M: float = 0.0
    C_dagger: float = 0.0
    mono_stretch_zero: bool = True
    C_near_pair: float = 0.0
    near_pair: tuple[int, int] | list[int] | None = None
    C_far_pair: float = 0.0
    far_pair: tuple[int, int] | list[int] | None = None
    C_bishell_shor_max_sample: float = 0.0
    n_pairs_scanned: int = 0
    closes_c0007: bool = False
    clay_implication: str = (
        "None. Shell-diff mono stretch=0 and bi-shell Shor sample; "
        "uniform C_† not reached; not continuum; not Clay."
    )
    notes: str = ""

    def as_dict(self) -> dict:
        return asdict(self)


def lemma_l0037(
    n: int = 24,
    iters: int = 16,
    max_pairs: int | None = 80,
) -> GalerkinBoundL0037:
    b7 = lemma_l0027(n=n, empirical=False)
    # Witness pairs: near (1,2) and far (2,134) if present, else from scan
    by = shell_modes_by_r(n)
    near = (1, 2) if 1 in by and 2 in by else None
    far = (2, 134) if 2 in by and 134 in by else None
    c_near = 0.0
    c_far = 0.0
    if near:
        c_near = float(bishell_shor_C_lower(n, near[0], near[1], iters=iters)["C_shor_lower"])
    if far:
        c_far = float(bishell_shor_C_lower(n, far[0], far[1], iters=iters)["C_shor_lower"])
    rows = scan_bishell_shor(n=n, iters=iters, max_pairs=max_pairs)
    c_max = max(
        [c_near, c_far] + [float(r["C_shor_lower"]) for r in rows]
    )
    if rows and float(rows[0]["C_shor_lower"]) >= c_far - 1e-9:
        top = rows[0]
        far = (int(top["r1"]), int(top["r2"]))
        c_far = float(top["C_shor_lower"])
    notes = (
        f"Mono-radial stretch=0 (shell-diff). Bi-shell Shor: near{near} C≳{c_near:.4g}; "
        f"far{far} C≳{c_far:.4g}; sample max C≳{c_max:.4g} over {len(rows)} pairs "
        f"(C_†={b7.C_dagger:.4g}). Closes C-0007: False."
    )
    return GalerkinBoundL0037(
        n=n,
        c0007_M=b7.c0007_M,
        C_dagger=b7.C_dagger,
        mono_stretch_zero=True,
        C_near_pair=c_near,
        near_pair=list(near) if near else None,
        C_far_pair=c_far,
        far_pair=list(far) if far else None,
        C_bishell_shor_max_sample=c_max,
        n_pairs_scanned=len(rows),
        closes_c0007=False,
        notes=notes,
    )


def save_lemma_l0037(
    bound: GalerkinBoundL0037,
    path: str | Path = "conjectures/proved_restricted/L-0037.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(bound.as_dict(), indent=2), encoding="utf-8")
