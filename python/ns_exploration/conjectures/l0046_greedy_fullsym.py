"""
L-0046: Greedy high-shell enlargement of physical fullsym (+ C-R-0011).

L-0045 closed consecutive |k|²≤25 (C≈9.119). Streaming fullsym + greedy
addition of sparse high shells (smallest D first, M≲8GB) yields a 30-shell
support including the spectral top shell 147:

  C_fullsym ≈ 9.317 ≤ C_† ≈ 9.562.

Frozen witness SUPPORT_CR0011. Further shells skipped by memory cap (M>8GB).

C-R-0011: Hermitian fields with Fourier support in SUPPORT_CR0011
⇒ Ω(0.02)≤M via L-0046+L-0026+L-0027.

All-IC still open (full dealias D≈6748).

FINITE dealias Galerkin only. Not continuum. Not Clay.
Evidence: N7 (streamed fullsym SVD) + L-0044/45.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path

from ns_exploration.conjectures.l0026_c0007_techo import lemma_l0026
from ns_exploration.conjectures.l0027_low_slab_cubic import lemma_l0027
from ns_exploration.conjectures.l0045_physical_band import SUPPORT_CR0010
from ns_exploration.conjectures.l0045_streaming_fullsym import streaming_C_fullsym

# Greedy enlargement of SUPPORT_CR0010 (frozen witness).
SUPPORT_CR0011 = (
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
    32,
    37,
    40,
    48,
    72,
    98,
    108,
    147,
)


@dataclass
class GalerkinBoundL0046:
    lemma_id: str = "L-0046"
    route: str = "B"
    status: str = "proved_restricted"
    evidence_level: str = "N7"
    n: int = 24
    c0007_M: float = 0.0
    C_dagger: float = 0.0
    Omega_star: float = 0.0
    support_cr0011: list[int] | None = None
    support_cr0011_D: int = 0
    C_fullsym_cr0011: float = 0.0
    C_fullsym_cr0010: float = 0.0
    closes_cr0011: bool = False
    closes_c0007_all_ic: bool = False
    clay_implication: str = (
        "None. Physical fullsym on 30-shell support incl. 147; C-R-0011 only; "
        "not all-IC; not continuum; not Clay."
    )
    notes: str = ""

    def as_dict(self) -> dict:
        return asdict(self)


def lemma_l0046(n: int = 24) -> GalerkinBoundL0046:
    b7 = lemma_l0027(n=n, empirical=False)
    b6 = lemma_l0026(n=n, n_grid=21)
    Cd = b7.C_dagger
    d10 = streaming_C_fullsym(n, SUPPORT_CR0010)
    d11 = streaming_C_fullsym(n, SUPPORT_CR0011)
    assert set(SUPPORT_CR0010).issubset(set(SUPPORT_CR0011))
    assert 147 in SUPPORT_CR0011
    assert d11["C_fullsym"] <= Cd + 1e-9
    assert d11["C_fullsym"] + 1e-9 >= d10["C_fullsym"]
    notes = (
        f"Greedy fullsym: C-R-0010 C={d10['C_fullsym']:.4g}; "
        f"C-R-0011 ({len(SUPPORT_CR0011)} shells, incl. 147) "
        f"C={d11['C_fullsym']:.4g}≤C_† (D={d11['D']}). All-IC: False."
    )
    return GalerkinBoundL0046(
        n=n,
        c0007_M=b7.c0007_M,
        C_dagger=Cd,
        Omega_star=b6.Omega_star,
        support_cr0011=list(SUPPORT_CR0011),
        support_cr0011_D=int(d11["D"]),
        C_fullsym_cr0011=float(d11["C_fullsym"]),
        C_fullsym_cr0010=float(d10["C_fullsym"]),
        closes_cr0011=True,
        closes_c0007_all_ic=False,
        notes=notes,
    )


def cr0011_record(bound: GalerkinBoundL0046) -> dict:
    return {
        "id": "C-R-0011",
        "route": "B",
        "status": "proved_restricted",
        "evidence_level": "N7",
        "lemma": "L-0046+L-0026+L-0027",
        "domain": (
            "Pseudospectral T^3, N<=24, 2/3 dealias, real (Hermitian) Fourier "
            f"fields with support in shells |k|^2 in {list(SUPPORT_CR0011)}. "
            "NOT continuum."
        ),
        "nu": 0.1,
        "energy": 0.5,
        "t_end": 0.02,
        "n_max": bound.n,
        "proposed_bound_M": bound.c0007_M,
        "shells": list(SUPPORT_CR0011),
        "C_fullsym": bound.C_fullsym_cr0011,
        "C_dagger": bound.C_dagger,
        "D": bound.support_cr0011_D,
        "Omega_star": bound.Omega_star,
        "parent_open": "C-0007",
        "statement": (
            f"For real divergence-free fields on N≤{bound.n} with Fourier support "
            f"only on shells |k|²∈{list(SUPPORT_CR0011)}, dealiased Galerkin NS "
            f"(ν=0.1, E≤0.5) satisfies Ω(0.02)≤{bound.c0007_M}. Proof: streaming "
            f"physical fullsym Shor gives C≤{bound.C_fullsym_cr0011}≤C_† ⇒ L-0027; "
            "high slab via L-0026. FINITE only."
        ),
        "clay_implication": (
            "None. 30-shell Hermitian subclass incl. top shell 147; "
            "not all-IC C-0007; not Clay."
        ),
        "notes": (
            "Grown greedily from C-R-0010 by adding smallest high shells under "
            "an 8GB dense-M memory cap. All-IC still open."
        ),
    }


def save_lemma_l0046(
    bound: GalerkinBoundL0046,
    path: str | Path = "conjectures/proved_restricted/L-0046.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(bound.as_dict(), indent=2), encoding="utf-8")


def save_cr0011(
    bound: GalerkinBoundL0046,
    path: str | Path = "conjectures/proved_restricted/C-R-0011.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(cr0011_record(bound), indent=2), encoding="utf-8")
