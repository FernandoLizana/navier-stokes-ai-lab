"""
L-0036: Signed Ω-weighted stretch matricization (Shor) techo for C_† / C-0007.

Write u in an Ω-weighted polarization basis z_α (L-0021-style dealias modes),
‖z‖₂² = 2Ω. The enstrophy stretch is a homogeneous cubic T(z,z,z). Matricize
  T(z,z,z) = z · L(zzᵀ)
with L : Sym → ℝ^D. Then |stretch| ≤ ‖L‖_{F→2} ‖z‖³, hence
  C ≤ C_Shor := 2√2 ‖L‖_{F→2}.
Power iteration on L*L yields a rigorous lower bound on ‖L‖ (N2 float; the
true ‖L‖ is at least this Rayleigh quotient).

On dealias grids:
  N=12:  C_Shor ≳ 42.8  > C_† ≈ 9.56
  N=24:  C_Shor ≳ 141.8 ≫ C_†
(already after a few power iterations). Therefore any C-0007 proof that
routes through this Shor / matricization majorant cannot meet C_†.

The ambient polarization space contains the physical real fields as a
subspace, so C ≤ C_Shor remains a valid (but too loose) majorant on real ICs.
Rank-1 / empirical stretch remains ≪ C_†; the gap is the Shor relaxation,
not the absence of sign in the cubic.

Still missing: tighter than Shor (SOS / chordal / rank-constrained) or a
shell-difference block lemma that keeps cancellation.

FINITE dealias Galerkin only. Not continuum. Not Clay.
Evidence: N7 (construction) + N2 (power iteration lower bound on ‖L‖).
"""

from __future__ import annotations

import json
import math
from collections import defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from ns_exploration.conjectures.l0021_quartic_shell import _leray_vec, _pol_basis
from ns_exploration.conjectures.l0027_low_slab_cubic import lemma_l0027
from ns_exploration.spectral.dealias import dealias_mask
from ns_exploration.spectral.fourier_conventions import wave_number_grids


def build_omega_weighted_basis(n: int) -> list[tuple[np.ndarray, np.ndarray, float]]:
    mask = dealias_mask(n)
    kx, ky, kz = wave_number_grids(n)
    basis: list[tuple[np.ndarray, np.ndarray, float]] = []
    for ix, iy, iz in np.argwhere(mask & ~((kx == 0) & (ky == 0) & (kz == 0))):
        k = (int(kx[ix, iy, iz]), int(ky[ix, iy, iz]), int(kz[ix, iy, iz]))
        lam = float(k[0] * k[0] + k[1] * k[1] + k[2] * k[2])
        for p in _pol_basis(np.array(k, float)):
            basis.append((np.array(k, int), np.asarray(p, float), lam))
    return basis


def signed_flattening_C_lower(
    n: int = 12,
    iters: int = 12,
    seed: int = 0,
) -> dict:
    """
    Lower-bound C_Shor = 2√2 ‖L‖_{F→2} by power iteration.
    Returns dict with D, n_pairs, L_op_lower, C_shor_lower.
    """
    basis = build_omega_weighted_basis(n)
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
        nrm = float(np.linalg.norm(Z2))
        Z = Z2 / (nrm + 1e-30)
        est = float(np.linalg.norm(apply_L(Z)))
    C_lo = 2.0 * math.sqrt(2.0) * est
    return {
        "n": n,
        "D": D,
        "n_pairs": int(len(I_a)),
        "L_op_lower": est,
        "C_shor_lower": C_lo,
        "iters": int(iters),
    }


@dataclass
class GalerkinBoundL0036:
    lemma_id: str = "L-0036"
    route: str = "B"
    status: str = "proved_restricted"
    evidence_level: str = "N7"
    n: int = 12
    c0007_M: float = 0.0
    C_dagger: float = 0.0
    D: int = 0
    n_pairs: int = 0
    L_op_lower: float = 0.0
    C_shor_lower: float = 0.0
    C_shor_lower_N24: float = 0.0
    closes_c0007: bool = False
    clay_implication: str = (
        "None. Signed Shor/matricization stretch majorant already exceeds C_†; "
        "not continuum; not Clay."
    )
    notes: str = ""

    def as_dict(self) -> dict:
        return asdict(self)


def lemma_l0036(
    n: int = 12,
    iters: int = 12,
    include_n24_probe: bool = False,
    n24_iters: int = 3,
) -> GalerkinBoundL0036:
    b7 = lemma_l0027(n=24, empirical=False)
    d = signed_flattening_C_lower(n=n, iters=iters)
    c24 = 0.0
    if include_n24_probe:
        d24 = signed_flattening_C_lower(n=24, iters=n24_iters)
        c24 = float(d24["C_shor_lower"])
    notes = (
        f"Signed Shor stretch: N={n} D={d['D']}, ‖L‖≳{d['L_op_lower']:.4g}, "
        f"C_Shor≳{d['C_shor_lower']:.4g} > C_†={b7.C_dagger:.4g}. "
        + (f"N=24 probe C_Shor≳{c24:.4g}. " if include_n24_probe else "")
        + "Closes C-0007: False."
    )
    return GalerkinBoundL0036(
        n=int(d["n"]),
        c0007_M=b7.c0007_M,
        C_dagger=b7.C_dagger,
        D=int(d["D"]),
        n_pairs=int(d["n_pairs"]),
        L_op_lower=float(d["L_op_lower"]),
        C_shor_lower=float(d["C_shor_lower"]),
        C_shor_lower_N24=c24,
        closes_c0007=False,
        notes=notes,
    )


def save_lemma_l0036(
    bound: GalerkinBoundL0036,
    path: str | Path = "conjectures/proved_restricted/L-0036.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(bound.as_dict(), indent=2), encoding="utf-8")
