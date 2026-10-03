"""
L-0044: Physical Hermitian stretch tensor + full symmetrization (+ C-R-0009).

Sprint A/B audit finding:
  The l0040 Galerkin tensor uses real amplitudes on each complex mode and omits
  the Fourier derivative factor i. The physical stretch ⟨ω, curl N⟩ for real
  (Hermitian) fields uses N_s = -P_s Σ_{p+q=s} i (û_p·q) û_q.

On the real Hermitian basis (c,s) per half-space mode×pol, with ‖z‖²=2Ω,
  stretch_inner(uh(z)) = f(z) = Σ G[m,a,b] z_m z_a z_b
matches the FFT path to ~1e-14.

Full symmetrization A = avg_6π G preserves f and yields the Shor majorant
  C_fullsym = 2√2 ‖flatten(A)‖_op.
On consecutive shells {1,2,3,4,5,6} (N=24):
  C_sym (i,j only) ≈ 11.344 > C_†
  C_fullsym        ≈ 2.993  ≤ C_† ≈ 9.562
hence the low-slab cubic of L-0027 closes on that Fourier support.

C-R-0009: Hermitian fields with Fourier support in |k|²∈{1,2,3,4,5,6}
⇒ Ω(0.02)≤M via L-0044+L-0026+L-0027.

FINITE dealias Galerkin only. Not continuum. Not Clay.
Evidence: N7 (SVD of exact triad tensor) + N2 (FFT cross-check).
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path

from ns_exploration.conjectures.l0026_c0007_techo import lemma_l0026
from ns_exploration.conjectures.l0027_low_slab_cubic import lemma_l0027
from ns_exploration.experiments.sprintA_stretch_tensor import (
    f_from_G,
    full_sym,
    operator_C,
    sym_ij,
)
from ns_exploration.experiments.sprintB_physical_tensor import (
    build_physical_G,
    z_to_uhat,
)
from ns_exploration.conjectures.l0023_cubic_dissipation import stretch_inner
import numpy as np

SUPPORT_CR0009 = (1, 2, 3, 4, 5, 6)


def physical_band_C_fullsym(n: int, radii: tuple[int, ...]) -> dict:
    G, basis = build_physical_G(n, radii)
    A = full_sym(G)
    Gs = sym_ij(G)
    return {
        "n": n,
        "radii": list(radii),
        "D": basis.D,
        "C_sym": operator_C(Gs),
        "C_fullsym": operator_C(A),
        "G": G,
        "A": A,
        "basis": basis,
    }


def fft_tensor_max_err(n: int, radii: tuple[int, ...], n_probe: int = 100) -> float:
    d = physical_band_C_fullsym(n, radii)
    G, basis = d["G"], d["basis"]
    rng = np.random.default_rng(44)
    err = 0.0
    for _ in range(n_probe):
        z = rng.standard_normal(basis.D)
        st = stretch_inner(z_to_uhat(z, basis))
        ft = f_from_G(G, z)
        scale = max(1.0, abs(st), abs(ft))
        err = max(err, abs(st - ft) / scale)
    return err


@dataclass
class GalerkinBoundL0044:
    lemma_id: str = "L-0044"
    route: str = "B"
    status: str = "proved_restricted"
    evidence_level: str = "N7"
    n: int = 24
    c0007_M: float = 0.0
    C_dagger: float = 0.0
    Omega_star: float = 0.0
    support_cr0009: list[int] | None = None
    support_cr0009_D: int = 0
    C_sym_123456: float = 0.0
    C_fullsym_123456: float = 0.0
    fft_tensor_err: float = 0.0
    closes_cr0009: bool = False
    closes_c0007_all_ic: bool = False
    clay_implication: str = (
        "None. Physical Hermitian fullsym Shor on shells {1..6}; C-R-0009 only; "
        "not all-IC; not continuum; not Clay."
    )
    notes: str = ""

    def as_dict(self) -> dict:
        return asdict(self)


def lemma_l0044(n: int = 24) -> GalerkinBoundL0044:
    b7 = lemma_l0027(n=n, empirical=False)
    b6 = lemma_l0026(n=n, n_grid=21)
    Cd = b7.C_dagger
    d = physical_band_C_fullsym(n, SUPPORT_CR0009)
    err = fft_tensor_max_err(n, SUPPORT_CR0009, n_probe=80)
    assert d["C_fullsym"] <= Cd + 1e-9
    assert d["C_sym"] > Cd
    assert err < 1e-11
    notes = (
        f"Physical Hermitian stretch on {list(SUPPORT_CR0009)}: "
        f"C_fullsym={d['C_fullsym']:.4g}≤C_†; C_sym={d['C_sym']:.4g}>C_†; "
        f"FFT↔tensor err={err:.2e}. C-R-0009. All-IC: False."
    )
    return GalerkinBoundL0044(
        n=n,
        c0007_M=b7.c0007_M,
        C_dagger=Cd,
        Omega_star=b6.Omega_star,
        support_cr0009=list(SUPPORT_CR0009),
        support_cr0009_D=int(d["D"]),
        C_sym_123456=float(d["C_sym"]),
        C_fullsym_123456=float(d["C_fullsym"]),
        fft_tensor_err=float(err),
        closes_cr0009=True,
        closes_c0007_all_ic=False,
        notes=notes,
    )


def cr0009_record(bound: GalerkinBoundL0044) -> dict:
    return {
        "id": "C-R-0009",
        "route": "B",
        "status": "proved_restricted",
        "evidence_level": "N7",
        "lemma": "L-0044+L-0026+L-0027",
        "domain": (
            "Pseudospectral T^3, N<=24, 2/3 dealias, real (Hermitian) Fourier "
            f"fields with support in shells |k|^2 in {list(SUPPORT_CR0009)}. "
            "NOT continuum."
        ),
        "nu": 0.1,
        "energy": 0.5,
        "t_end": 0.02,
        "n_max": bound.n,
        "proposed_bound_M": bound.c0007_M,
        "shells": list(SUPPORT_CR0009),
        "C_fullsym": bound.C_fullsym_123456,
        "C_sym": bound.C_sym_123456,
        "C_dagger": bound.C_dagger,
        "D": bound.support_cr0009_D,
        "Omega_star": bound.Omega_star,
        "parent_open": "C-0007",
        "statement": (
            f"For real divergence-free fields on N≤{bound.n} with Fourier support "
            f"only on shells |k|²∈{list(SUPPORT_CR0009)}, dealiased Galerkin NS "
            f"(ν=0.1, E≤0.5) satisfies Ω(0.02)≤{bound.c0007_M}. Proof: physical "
            f"Hermitian stretch tensor + full symmetrization gives "
            f"C≤{bound.C_fullsym_123456}≤C_† ⇒ L-0027; high slab via L-0026. "
            "FINITE only."
        ),
        "clay_implication": (
            "None. Low-shell Hermitian subclass only; not all-IC C-0007; not Clay."
        ),
        "notes": (
            "Uses the i-factor physical triad tensor (matches stretch_inner). "
            "Exploratory: C_fullsym remains ≤C_† at least through shells r≤20 "
            "(D=776, C≈7.70); dense G OOMs beyond. All-IC still open."
        ),
    }


def save_lemma_l0044(
    bound: GalerkinBoundL0044,
    path: str | Path = "conjectures/proved_restricted/L-0044.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(bound.as_dict(), indent=2), encoding="utf-8")


def save_cr0009(
    bound: GalerkinBoundL0044,
    path: str | Path = "conjectures/proved_restricted/C-R-0009.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(cr0009_record(bound), indent=2), encoding="utf-8")
