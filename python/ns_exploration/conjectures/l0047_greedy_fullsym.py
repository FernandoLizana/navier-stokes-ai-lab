"""
L-0047: Extended greedy fullsym (+ C-R-0012).

L-0046 closed 30 shells (C≈9.317). Batch-greedy under float32 dense-M (≤12GB)
+ Gram op-norm enlarges to 41 shells incl. top shell 147:

  C_fullsym ≈ 9.5037 ≤ C_† ≈ 9.562.

Frozen witness SUPPORT_CR0012. Further shells either reject (C>C_†) or
exceed the 12GB float32 M cap at current D≈1772.

C-R-0012: Hermitian fields with Fourier support in SUPPORT_CR0012
⇒ Ω(0.02)≤M via L-0047+L-0026+L-0027.

All-IC still open (full dealias D≈6748; many mid shells skipped by RAM).

FINITE dealias Galerkin only. Not continuum. Not Clay.
Evidence: N7 (streamed fullsym + Gram σmax, float32 M / float64 Gram) + L-0044/46.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from ns_exploration.conjectures.l0026_c0007_techo import lemma_l0026
from ns_exploration.conjectures.l0027_low_slab_cubic import lemma_l0027
from ns_exploration.conjectures.l0045_streaming_fullsym import streaming_C_fullsym
from ns_exploration.conjectures.l0046_greedy_fullsym import SUPPORT_CR0011

# Frozen witness from reports/l0047_greedy_batch_result.json (done=true).
SUPPORT_CR0012 = (
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
    43,
    44,
    48,
    52,
    58,
    67,
    68,
    72,
    73,
    76,
    88,
    98,
    102,
    108,
    114,
    147,
)


@dataclass
class GalerkinBoundL0047:
    lemma_id: str = "L-0047"
    route: str = "B"
    status: str = "proved_restricted"
    evidence_level: str = "N7"
    n: int = 24
    c0007_M: float = 0.0
    C_dagger: float = 0.0
    Omega_star: float = 0.0
    support_cr0012: list[int] | None = None
    support_cr0012_D: int = 0
    C_fullsym_cr0012: float = 0.0
    C_fullsym_cr0011: float = 0.0
    closes_cr0012: bool = False
    closes_c0007_all_ic: bool = False
    clay_implication: str = (
        "None. Physical fullsym on 41-shell support incl. 147; C-R-0012 only; "
        "not all-IC; not continuum; not Clay."
    )
    notes: str = ""

    def as_dict(self) -> dict:
        return asdict(self)


def lemma_l0047(n: int = 24) -> GalerkinBoundL0047:
    b7 = lemma_l0027(n=n, empirical=False)
    b6 = lemma_l0026(n=n, n_grid=21)
    Cd = b7.C_dagger
    assert set(SUPPORT_CR0011).issubset(set(SUPPORT_CR0012))
    assert 147 in SUPPORT_CR0012
    # D≈1772: float64 M OOMs; float32 M + float64 Gram eig (validated vs SVD).
    d11 = streaming_C_fullsym(n, SUPPORT_CR0011, dtype=np.float32)
    d12 = streaming_C_fullsym(n, SUPPORT_CR0012, dtype=np.float32)
    # Margin to C_† ≈0.058 ≫ float32 relative error (~3e-8 on small supports).
    assert d12["C_fullsym"] <= Cd + 1e-5
    assert d12["C_fullsym"] + 1e-5 >= d11["C_fullsym"]
    notes = (
        f"Greedy fullsym: C-R-0011 C={d11['C_fullsym']:.4g}; "
        f"C-R-0012 ({len(SUPPORT_CR0012)} shells, incl. 147) "
        f"C={d12['C_fullsym']:.4g}≤C_† (D={d12['D']}, float32 M). All-IC: False."
    )
    return GalerkinBoundL0047(
        n=n,
        c0007_M=b7.c0007_M,
        C_dagger=Cd,
        Omega_star=b6.Omega_star,
        support_cr0012=list(SUPPORT_CR0012),
        support_cr0012_D=int(d12["D"]),
        C_fullsym_cr0012=float(d12["C_fullsym"]),
        C_fullsym_cr0011=float(d11["C_fullsym"]),
        closes_cr0012=True,
        closes_c0007_all_ic=False,
        notes=notes,
    )


def cr0012_record(bound: GalerkinBoundL0047) -> dict:
    return {
        "id": "C-R-0012",
        "route": "B",
        "status": "proved_restricted",
        "evidence_level": "N7",
        "lemma": "L-0047+L-0026+L-0027",
        "domain": (
            "Pseudospectral T^3, N<=24, 2/3 dealias, real (Hermitian) Fourier "
            f"fields with support in shells |k|^2 in {list(SUPPORT_CR0012)}. "
            "NOT continuum."
        ),
        "nu": 0.1,
        "energy": 0.5,
        "t_end": 0.02,
        "n_max": bound.n,
        "proposed_bound_M": bound.c0007_M,
        "shells": list(SUPPORT_CR0012),
        "C_fullsym": bound.C_fullsym_cr0012,
        "C_dagger": bound.C_dagger,
        "D": bound.support_cr0012_D,
        "Omega_star": bound.Omega_star,
        "parent_open": "C-0007",
        "statement": (
            f"For real divergence-free fields on N≤{bound.n} with Fourier support "
            f"only on shells |k|²∈{list(SUPPORT_CR0012)}, dealiased Galerkin NS "
            f"(ν=0.1, E≤0.5) satisfies Ω(0.02)≤{bound.c0007_M}. Proof: streaming "
            f"physical fullsym Shor (float32 M, float64 Gram σmax) gives "
            f"C≤{bound.C_fullsym_cr0012}≤C_† ⇒ L-0027; high slab via L-0026. "
            "FINITE only."
        ),
        "clay_implication": (
            "None. 41-shell Hermitian subclass incl. top shell 147; "
            "not all-IC C-0007; not Clay."
        ),
        "notes": (
            "Grown greedily from C-R-0011 by batch-adding high shells under a "
            "12GB float32 dense-M cap + Gram op-norm. Rejected 85,97,107 (C>C_†); "
            "remaining mid shells skipped by RAM at D≈1772. All-IC still open."
        ),
    }


def save_lemma_l0047(
    bound: GalerkinBoundL0047,
    path: str | Path = "conjectures/proved_restricted/L-0047.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(bound.as_dict(), indent=2), encoding="utf-8")


def save_cr0012(
    bound: GalerkinBoundL0047,
    path: str | Path = "conjectures/proved_restricted/C-R-0012.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(cr0012_record(bound), indent=2), encoding="utf-8")
