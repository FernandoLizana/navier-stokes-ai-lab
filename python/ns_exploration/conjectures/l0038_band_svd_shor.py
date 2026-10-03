"""
L-0038: Exact SVD Shor on low-shell bands + C-R-0004.

On a fixed set of radial shells R, build the Ω-weighted polarization basis
(L-0021 / L-0036 style) and the stretch matricization L : Sym → ℝ^D with
  stretch = z · L(zzᵀ),  ‖z‖₂² = 2Ω.
Then C ≤ C_Shor(R) := 2√2 ‖L‖_{op}. For small |R| we form the dense
matrix of L (size D × D²) and compute ‖L‖_op by SVD (exact in float64).

On N=24 dealias:
  R={1,2}     → C_Shor ≈ 6.633 ≤ C_†
  R={1,2,3}   → C_Shor ≈ 9.165 ≤ C_†
  R={1,2,5}   → C_Shor ≈ 9.550 ≤ C_†
  R={1,2,3,4} → C_Shor ≈ 10.066 > C_†  (techo for consecutive K≥4)

Together with L-0026 (high Ω0) and L-0027 (low slab needs C≤C_†), this proves
C-0007 on the IC subclass of fields supported on shells {1,2,3}
(C-R-0004). All-IC C-0007 remains open: larger bands fail the Shor test.

FINITE dealias Galerkin only. Not continuum. Not Clay.
Evidence: N7 (SVD arithmetic on finite matrix) + L-0026/L-0027.
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
from ns_exploration.conjectures.l0026_c0007_techo import lemma_l0026
from ns_exploration.conjectures.l0027_low_slab_cubic import lemma_l0027
from ns_exploration.spectral.dealias import dealias_mask
from ns_exploration.spectral.fourier_conventions import wave_number_grids


def exact_band_C_shor(n: int, radii: tuple[int, ...] | list[int]) -> dict:
    """Exact (float64 SVD) C_Shor on the Ω-weighted span of the given shells."""
    by = shell_modes_by_r(n)
    radii_t = tuple(int(r) for r in radii)
    for r in radii_t:
        if r not in by:
            raise KeyError(f"shell {r} missing on n={n}")
    modes: list[tuple[int, int, int]] = []
    for r in radii_t:
        modes.extend(list(by[r]))
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
    by_s: dict[int, list[tuple[int, np.ndarray, float]]] = defaultdict(list)
    for a, (k, p, lam) in enumerate(basis):
        by_s[out[tuple(int(x) for x in k)]].append((a, p, lam))

    Lmat = np.zeros((D, D * D), dtype=np.float64)
    n_pairs = 0
    for i, (ki, ai, li) in enumerate(basis):
        for j, (kj, aj, lj) in enumerate(basis):
            s = (int(ki[0] + kj[0]), int(ki[1] + kj[1]), int(ki[2] + kj[2]))
            if s not in out:
                continue
            n_pairs += 1
            coef = float(np.dot(ai, kj.astype(float))) / math.sqrt(li * lj)
            Pv = _leray_vec(np.array(s, float), aj)
            vec = coef * Pv
            col = i * D + j
            for a, p, lam in by_s[out[s]]:
                Lmat[a, col] += math.sqrt(lam) * float(np.dot(p, vec))
    Lop = float(np.linalg.svd(Lmat, compute_uv=False)[0]) if D > 0 else 0.0
    C = 2.0 * math.sqrt(2.0) * Lop
    return {
        "n": n,
        "radii": list(radii_t),
        "D": D,
        "n_pairs": n_pairs,
        "L_op": Lop,
        "C_shor": C,
    }


@dataclass
class GalerkinBoundL0038:
    lemma_id: str = "L-0038"
    route: str = "B"
    status: str = "proved_restricted"
    evidence_level: str = "N7"
    n: int = 24
    c0007_M: float = 0.0
    C_dagger: float = 0.0
    band_123_C: float = 0.0
    band_1234_C: float = 0.0
    band_12_C: float = 0.0
    band_125_C: float = 0.0
    meets_Cdagger_bands: list | None = None
    fails_Cdagger_bands: list | None = None
    cr0004_shells: list | None = None
    closes_c0007_all_ic: bool = False
    closes_cr0004: bool = True
    clay_implication: str = (
        "None. Exact SVD Shor on low-shell bands; C-R-0004 subclass only; "
        "not all-IC; not continuum; not Clay."
    )
    notes: str = ""

    def as_dict(self) -> dict:
        return asdict(self)


def lemma_l0038(n: int = 24) -> GalerkinBoundL0038:
    b7 = lemma_l0027(n=n, empirical=False)
    bands = {
        "12": exact_band_C_shor(n, (1, 2)),
        "123": exact_band_C_shor(n, (1, 2, 3)),
        "1234": exact_band_C_shor(n, (1, 2, 3, 4)),
        "125": exact_band_C_shor(n, (1, 2, 5)),
    }
    Cd = b7.C_dagger
    meets = []
    fails = []
    for key, d in bands.items():
        entry = {"radii": d["radii"], "C_shor": d["C_shor"], "D": d["D"]}
        if d["C_shor"] <= Cd + 1e-9:
            meets.append(entry)
        else:
            fails.append(entry)
    notes = (
        f"Exact SVD Shor: {{1,2}} C={bands['12']['C_shor']:.4g}; "
        f"{{1,2,3}} C={bands['123']['C_shor']:.4g}≤C_†; "
        f"{{1,2,5}} C={bands['125']['C_shor']:.4g}≤C_†; "
        f"{{1,2,3,4}} C={bands['1234']['C_shor']:.4g}>C_†={Cd:.4g}. "
        f"C-R-0004 on {{1,2,3}} via L-0026+L-0027. All-IC: False."
    )
    return GalerkinBoundL0038(
        n=n,
        c0007_M=b7.c0007_M,
        C_dagger=Cd,
        band_123_C=float(bands["123"]["C_shor"]),
        band_1234_C=float(bands["1234"]["C_shor"]),
        band_12_C=float(bands["12"]["C_shor"]),
        band_125_C=float(bands["125"]["C_shor"]),
        meets_Cdagger_bands=meets,
        fails_Cdagger_bands=fails,
        cr0004_shells=[1, 2, 3],
        closes_c0007_all_ic=False,
        closes_cr0004=True,
        notes=notes,
    )


def build_cr0004(n: int = 24) -> dict:
    b7 = lemma_l0027(n=n, empirical=False)
    b26 = lemma_l0026(n=n, n_grid=21)
    d = exact_band_C_shor(n, (1, 2, 3))
    assert d["C_shor"] <= b7.C_dagger + 1e-9
    return {
        "id": "C-R-0004",
        "route": "B",
        "status": "proved_restricted",
        "evidence_level": "N7",
        "lemma": "L-0038+L-0026+L-0027",
        "domain": (
            "Pseudospectral T^3, N<=24, 2/3 dealias, Fourier support in shells "
            "|k|^2 in {1,2,3}. NOT continuum."
        ),
        "nu": 0.1,
        "energy": 0.5,
        "t_end": 0.02,
        "n_max": n,
        "proposed_bound_M": b7.c0007_M,
        "shells": [1, 2, 3],
        "C_shor": d["C_shor"],
        "C_dagger": b7.C_dagger,
        "Omega_star": b26.Omega_star,
        "parent_open": "C-0007",
        "statement": (
            f"For divergence-free fields on N≤{n} with Fourier support only on "
            f"shells |k|²∈{{1,2,3}}, dealiased Galerkin NS (ν=0.1, E≤0.5) satisfies "
            f"Ω(0.02)≤{b7.c0007_M:.8f}. Proof: exact SVD Shor gives C≤{d['C_shor']:.8f}≤C_† "
            f"⇒ low slab via L-0027; high slab Ω0≥Ω★ via L-0026/L-0024. FINITE only."
        ),
        "clay_implication": (
            "None. Low-shell-band subclass only; not all-IC C-0007; not Clay."
        ),
        "notes": (
            "L-0038 SVD on D="
            f"{d['D']} polarization modes. Larger consecutive band {{1,2,3,4}} fails C_†."
        ),
    }


def save_lemma_l0038(
    bound: GalerkinBoundL0038,
    path: str | Path = "conjectures/proved_restricted/L-0038.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(bound.as_dict(), indent=2), encoding="utf-8")


def save_cr0004(
    d: dict,
    path: str | Path = "conjectures/proved_restricted/C-R-0004.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(d, indent=2), encoding="utf-8")
