"""
L-0029: Frobenius Jacobian stretch embedding + techo vs C_†.

The enstrophy production satisfies
  |⟨ω, curl N⟩| = |∫ ωᵀ (∇u) ω| ≤ ‖∇u‖_{∞,F} ‖ω‖₂² = 2 Ω ‖∇u‖_{∞,F},
where ‖∇u‖_{∞,F} := ess sup_x ‖Jacobian(u)(x)‖_Frobenius.

Fourier triangle inequality modewise:
  ‖J(x)‖_F ≤ Σ_k |k| ‖û_k‖₂,
and Cauchy–Schwarz over the M nonzero dealias modes,
  Σ_k |k| ‖û_k‖ ≤ √M √(Σ |k|² ‖û_k‖²) = √M √(2Ω).
Hence
  |stretch| ≤ 2 √(2M) Ω^{3/2},
  C_F := 2 √(2M).

This removes the √3 from L-0002's componentwise ℓ¹ bound
  C_{L0002} = 2√2 √(3M) = 2√(6M) = √3 · C_F.

On N=24 dealias, M=3374 ⇒ C_F ≈ 164.29 ≪ C_{L0002}≈284.56, but still
C_F ≫ C_†≈9.562 (L-0027). Meeting C_† via this embedding requires
  M ≤ M_† := C_†² / 8 ≈ 11.43.
The only radial hard truncation with M≤11 is the |k|²=1 shell (M=6),
whose maximal Ω/E equals 1 ≪ Ω★/E0≈37.75, so it cannot probe the
L-0027 low slab edge. Any Galerkin subspace that can realize Ω/E≥Ω★/E0
must include a mode with |k|²≥38 and — on the full dealias mask —
still uses M=3374 for this embedding.

Conclusion: the Frobenius improvement is real but structurally insufficient
for L-0027 / C-0007; a triad-aware cubic bound remains necessary.

FINITE dealias Galerkin only. Not continuum. Not Clay.
Evidence: N7 (embedding arithmetic).
"""

from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

from ns_exploration.conjectures.l0002_viscous_galerkin import stretch_constant
from ns_exploration.conjectures.l0027_low_slab_cubic import lemma_l0027
from ns_exploration.conjectures.l0028_embedding_techo import modes_upto_radius


def stretch_constant_frobenius(M: int) -> float:
    """C_F = 2 √(2M) in |stretch| ≤ C_F Ω^{3/2}."""
    if M <= 0:
        return 0.0
    return 2.0 * math.sqrt(2.0 * float(M))


def M_star_for_Cdagger_frobenius(C_dagger: float) -> float:
    """M with 2√(2M) = C_†  ⇒  M = C_†² / 8."""
    return float(C_dagger) ** 2 / 8.0


def max_lambda_upto(n: int, R: int) -> int:
    """Max |k|² among nonzero dealias modes with |k|² ≤ R."""
    return int(R)


@dataclass
class GalerkinBoundL0029:
    lemma_id: str = "L-0029"
    route: str = "B"
    status: str = "proved_restricted"
    evidence_level: str = "N7"
    n: int = 24
    C_dagger: float = 0.0
    Omega_star: float = 0.0
    lambda_star: float = 0.0
    M_full: int = 0
    C_F_full: float = 0.0
    C_L0002_full: float = 0.0
    improvement_factor: float = 0.0
    M_star_F: float = 0.0
    M_shell_r1: int = 0
    C_F_shell_r1: float = 0.0
    shell_r1_meets_Cdagger: bool = False
    shell_r1_reaches_low_slab_edge: bool = False
    frobenius_closes_c0007: bool = False
    clay_implication: str = (
        "None. Frobenius stretch embedding techo only; C_† still unproved "
        "on the full dealias class; not continuum; not Clay."
    )
    notes: str = ""

    def as_dict(self) -> dict:
        return asdict(self)


def lemma_l0029(n: int = 24) -> GalerkinBoundL0029:
    b7 = lemma_l0027(n=n, empirical=False)
    Cd = b7.C_dagger
    om_star = b7.Omega_star
    E0 = b7.E0
    lam_star = om_star / E0
    M_full = modes_upto_radius(n, 10**9)
    C_F = stretch_constant_frobenius(M_full)
    C_old = stretch_constant(M_full)
    Mstar = M_star_for_Cdagger_frobenius(Cd)
    M_r1 = modes_upto_radius(n, 1)
    C_r1 = stretch_constant_frobenius(M_r1)
    meets = C_r1 <= Cd + 1e-12
    reaches = 1.0 >= lam_star - 1e-12  # False
    notes = (
        f"C_F=2√(2M): full M={M_full} ⇒ C_F={C_F:.6g} "
        f"(L-0002 C={C_old:.6g}, factor √3≈{math.sqrt(3):.4f}). "
        f"C_†={Cd:.6g} needs M≤M_†^F={Mstar:.4f}. "
        f"r=1: M={M_r1}, C_F={C_r1:.6g} meets C_†={meets}, "
        f"λ_max=1 vs λ★={lam_star:.4f} reaches edge={reaches}. "
        f"Full-mask Frobenius closes C-0007: False."
    )
    return GalerkinBoundL0029(
        n=n,
        C_dagger=Cd,
        Omega_star=om_star,
        lambda_star=lam_star,
        M_full=M_full,
        C_F_full=C_F,
        C_L0002_full=C_old,
        improvement_factor=math.sqrt(3.0),
        M_star_F=Mstar,
        M_shell_r1=M_r1,
        C_F_shell_r1=C_r1,
        shell_r1_meets_Cdagger=meets,
        shell_r1_reaches_low_slab_edge=reaches,
        frobenius_closes_c0007=False,
        notes=notes,
    )


def save_lemma_l0029(
    bound: GalerkinBoundL0029,
    path: str | Path = "conjectures/proved_restricted/L-0029.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(bound.as_dict(), indent=2), encoding="utf-8")
