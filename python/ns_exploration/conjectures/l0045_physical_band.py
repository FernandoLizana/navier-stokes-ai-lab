"""
L-0045: Streaming physical fullsym — maximal consecutive band (+ C-R-0010).

L-0044 closed Hermitian shells {1..6} with dense G. Here we accumulate
M = flatten(full_sym(G)) without storing G (D³), enabling larger bands.

On N=24 dealias, physical Hermitian stretch:
  r_max=25 → C_fullsym ≈ 9.119 ≤ C_† ≈ 9.562  (closes)
  r_max=26 → C_fullsym ≈ 9.779 > C_†           (techo)

C-R-0010: real fields with Fourier support in all shells |k|² ≤ 25
(existing Z³ shells: 22 radii, D=1028).

All-IC (full dealias, D≈6748) still open — M would be ~1 TB dense.

FINITE dealias Galerkin only. Not continuum. Not Clay.
Evidence: N7 (SVD of streamed fullsym matricization) + L-0044 FFT audit.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path

from ns_exploration.conjectures.l0021_quartic_shell import shell_modes_by_r
from ns_exploration.conjectures.l0026_c0007_techo import lemma_l0026
from ns_exploration.conjectures.l0027_low_slab_cubic import lemma_l0027
from ns_exploration.conjectures.l0045_streaming_fullsym import streaming_C_fullsym

# Frozen witness: all existing shells with |k|² ≤ 25 on N=24 dealias.
SUPPORT_CR0010 = (
    1,
    2,
    3,
    4,
    5,
    6,
    8,
    9,
    10,
    11,
    12,
    13,
    14,
    16,
    17,
    18,
    19,
    20,
    21,
    22,
    24,
    25,
)
TECHO_RADII = SUPPORT_CR0010 + (26,)  # first consecutive failure


def radii_upto(n: int, rmax: int) -> tuple[int, ...]:
    by = shell_modes_by_r(n)
    return tuple(r for r in sorted(by.keys()) if r <= rmax)


@dataclass
class GalerkinBoundL0045:
    lemma_id: str = "L-0045"
    route: str = "B"
    status: str = "proved_restricted"
    evidence_level: str = "N7"
    n: int = 24
    c0007_M: float = 0.0
    C_dagger: float = 0.0
    Omega_star: float = 0.0
    support_cr0010: list[int] | None = None
    support_cr0010_D: int = 0
    C_fullsym_rmax25: float = 0.0
    C_fullsym_rmax26: float = 0.0
    closes_cr0010: bool = False
    closes_c0007_all_ic: bool = False
    clay_implication: str = (
        "None. Physical fullsym on shells |k|²≤25; C-R-0010 only; "
        "not all-IC; not continuum; not Clay."
    )
    notes: str = ""

    def as_dict(self) -> dict:
        return asdict(self)


def lemma_l0045(n: int = 24) -> GalerkinBoundL0045:
    b7 = lemma_l0027(n=n, empirical=False)
    b6 = lemma_l0026(n=n, n_grid=21)
    Cd = b7.C_dagger
    assert tuple(SUPPORT_CR0010) == radii_upto(n, 25)
    d25 = streaming_C_fullsym(n, SUPPORT_CR0010)
    d26 = streaming_C_fullsym(n, radii_upto(n, 26))
    assert d25["C_fullsym"] <= Cd + 1e-9
    assert d26["C_fullsym"] > Cd
    notes = (
        f"Streaming physical fullsym: r≤25 C={d25['C_fullsym']:.4g}≤C_† "
        f"(D={d25['D']}); r≤26 C={d26['C_fullsym']:.4g}>C_†. "
        f"C-R-0010 on {len(SUPPORT_CR0010)} shells. All-IC: False."
    )
    return GalerkinBoundL0045(
        n=n,
        c0007_M=b7.c0007_M,
        C_dagger=Cd,
        Omega_star=b6.Omega_star,
        support_cr0010=list(SUPPORT_CR0010),
        support_cr0010_D=int(d25["D"]),
        C_fullsym_rmax25=float(d25["C_fullsym"]),
        C_fullsym_rmax26=float(d26["C_fullsym"]),
        closes_cr0010=True,
        closes_c0007_all_ic=False,
        notes=notes,
    )


def cr0010_record(bound: GalerkinBoundL0045) -> dict:
    return {
        "id": "C-R-0010",
        "route": "B",
        "status": "proved_restricted",
        "evidence_level": "N7",
        "lemma": "L-0045+L-0026+L-0027",
        "domain": (
            "Pseudospectral T^3, N<=24, 2/3 dealias, real (Hermitian) Fourier "
            f"fields with support in shells |k|^2 in {list(SUPPORT_CR0010)} "
            "(|k|²≤25). NOT continuum."
        ),
        "nu": 0.1,
        "energy": 0.5,
        "t_end": 0.02,
        "n_max": bound.n,
        "proposed_bound_M": bound.c0007_M,
        "shells": list(SUPPORT_CR0010),
        "r_max": 25,
        "C_fullsym": bound.C_fullsym_rmax25,
        "C_dagger": bound.C_dagger,
        "D": bound.support_cr0010_D,
        "Omega_star": bound.Omega_star,
        "parent_open": "C-0007",
        "statement": (
            f"For real divergence-free fields on N≤{bound.n} with Fourier support "
            f"only on shells |k|²≤25 ({list(SUPPORT_CR0010)}), dealiased Galerkin NS "
            f"(ν=0.1, E≤0.5) satisfies Ω(0.02)≤{bound.c0007_M}. Proof: streaming "
            f"physical fullsym Shor gives C≤{bound.C_fullsym_rmax25}≤C_† ⇒ L-0027; "
            "high slab via L-0026. FINITE only."
        ),
        "clay_implication": (
            "None. Shells |k|²≤25 subclass only; not all-IC C-0007; not Clay."
        ),
        "notes": (
            f"Techo: r≤26 gives C_fullsym={bound.C_fullsym_rmax26:.4g}>C_†. "
            "Full dealias (D≈6748) still open."
        ),
    }


def save_lemma_l0045(
    bound: GalerkinBoundL0045,
    path: str | Path = "conjectures/proved_restricted/L-0045.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(bound.as_dict(), indent=2), encoding="utf-8")


def save_cr0010(
    bound: GalerkinBoundL0045,
    path: str | Path = "conjectures/proved_restricted/C-R-0010.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(cr0010_record(bound), indent=2), encoding="utf-8")
