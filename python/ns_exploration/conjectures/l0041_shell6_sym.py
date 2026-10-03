"""
L-0041: Shell-6 Sym catalog + consecutive-{1..6} techo + C-R-0007.

L-0040 closed consecutive shells {1..5} by Sym-restricted Shor (C_Sym≤C_†).
The next consecutive band {1..6} has
  C_Sym({1..6}) ≈ 11.344 > C_† ≈ 9.562
so Sym Shor cannot prove C-0007 on that full band (techo). Empirically
rank-1 stretch on {1..6} remains ≪ C_† (~0.26), so the gap is still the
Sym relaxation.

Leave-one-out on {1..6} shows several 5-shell subsets that *include* shell 6
and still meet C_†, e.g.
  {1,2,3,4,6} → C_Sym ≈ 8.406
  {2,3,4,5,6} → C_Sym ≈ 8.639
  {1,3,4,5,6} → C_Sym ≈ 9.374
Greedy Sym enlargement of {1,2,3,4,6} yields the 8-shell witness
  A6 = {1,2,3,4,6,11,12,25} with C_Sym ≈ 9.414 ≤ C_†,
proving C-0007 on that IC subclass (C-R-0007).

FINITE dealias Galerkin only. Not continuum. Not Clay.
Evidence: N7 (Sym SVD) + L-0026/L-0027.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path

from ns_exploration.conjectures.l0026_c0007_techo import lemma_l0026
from ns_exploration.conjectures.l0027_low_slab_cubic import lemma_l0027
from ns_exploration.conjectures.l0040_sym_shor import sym_band_C_shor

# Consecutive techo.
BAND_123456 = (1, 2, 3, 4, 5, 6)
# 5-shell seeds that include shell 6 and meet C_†.
SEED_12346 = (1, 2, 3, 4, 6)
SEED_23456 = (2, 3, 4, 5, 6)
SEED_13456 = (1, 3, 4, 5, 6)
# Greedy enlargement of SEED_12346 (frozen witness for C-R-0007).
SUPPORT_CR0007 = (1, 2, 3, 4, 6, 11, 12, 25)


@dataclass
class GalerkinBoundL0041:
    lemma_id: str = "L-0041"
    route: str = "B"
    status: str = "proved_restricted"
    evidence_level: str = "N7"
    n: int = 24
    c0007_M: float = 0.0
    C_dagger: float = 0.0
    band_123456_C: float = 0.0
    seed_12346_C: float = 0.0
    seed_23456_C: float = 0.0
    seed_13456_C: float = 0.0
    support_cr0007: list | None = None
    support_cr0007_C: float = 0.0
    support_cr0007_D: int = 0
    closes_cr0007: bool = True
    closes_c0007_all_ic: bool = False
    consecutive_6_techo: bool = True
    clay_implication: str = (
        "None. Sym Shor shell-6 catalog / consecutive-6 techo; C-R-0007 only; "
        "not all-IC; not continuum; not Clay."
    )
    notes: str = ""

    def as_dict(self) -> dict:
        return asdict(self)


def lemma_l0041(n: int = 24) -> GalerkinBoundL0041:
    b7 = lemma_l0027(n=n, empirical=False)
    Cd = b7.C_dagger
    d6 = sym_band_C_shor(n, BAND_123456)
    s1 = sym_band_C_shor(n, SEED_12346)
    s2 = sym_band_C_shor(n, SEED_23456)
    s3 = sym_band_C_shor(n, SEED_13456)
    cr = sym_band_C_shor(n, SUPPORT_CR0007)
    assert d6["C_shor_sym"] > Cd
    assert s1["C_shor_sym"] <= Cd + 1e-9
    assert s2["C_shor_sym"] <= Cd + 1e-9
    assert s3["C_shor_sym"] <= Cd + 1e-9
    assert cr["C_shor_sym"] <= Cd + 1e-8
    notes = (
        f"Sym: {{1..6}} C={d6['C_shor_sym']:.4g}>C_† (techo); "
        f"seeds with 6: {{1,2,3,4,6}}={s1['C_shor_sym']:.4g}, "
        f"{{2,3,4,5,6}}={s2['C_shor_sym']:.4g}, {{1,3,4,5,6}}={s3['C_shor_sym']:.4g}; "
        f"C-R-0007 support {list(SUPPORT_CR0007)} C={cr['C_shor_sym']:.4g}. All-IC: False."
    )
    return GalerkinBoundL0041(
        n=n,
        c0007_M=b7.c0007_M,
        C_dagger=Cd,
        band_123456_C=float(d6["C_shor_sym"]),
        seed_12346_C=float(s1["C_shor_sym"]),
        seed_23456_C=float(s2["C_shor_sym"]),
        seed_13456_C=float(s3["C_shor_sym"]),
        support_cr0007=list(SUPPORT_CR0007),
        support_cr0007_C=float(cr["C_shor_sym"]),
        support_cr0007_D=int(cr["D"]),
        closes_cr0007=True,
        closes_c0007_all_ic=False,
        consecutive_6_techo=True,
        notes=notes,
    )


def build_cr0007(n: int = 24) -> dict:
    b7 = lemma_l0027(n=n, empirical=False)
    b26 = lemma_l0026(n=n, n_grid=21)
    d = sym_band_C_shor(n, SUPPORT_CR0007)
    assert d["C_shor_sym"] <= b7.C_dagger + 1e-8
    return {
        "id": "C-R-0007",
        "route": "B",
        "status": "proved_restricted",
        "evidence_level": "N7",
        "lemma": "L-0041+L-0026+L-0027",
        "domain": (
            "Pseudospectral T^3, N<=24, 2/3 dealias, Fourier support in shells "
            f"{list(SUPPORT_CR0007)} (includes |k|^2=6). NOT continuum."
        ),
        "nu": 0.1,
        "energy": 0.5,
        "t_end": 0.02,
        "n_max": n,
        "proposed_bound_M": b7.c0007_M,
        "shells": list(SUPPORT_CR0007),
        "C_shor_sym": d["C_shor_sym"],
        "C_dagger": b7.C_dagger,
        "D": d["D"],
        "Omega_star": b26.Omega_star,
        "parent_open": "C-0007",
        "statement": (
            f"For divergence-free fields on N≤{n} with Fourier support in shells "
            f"{list(SUPPORT_CR0007)}, dealiased Galerkin NS (ν=0.1, E≤0.5) satisfies "
            f"Ω(0.02)≤{b7.c0007_M:.8f}. Proof: Sym SVD Shor gives "
            f"C≤{d['C_shor_sym']:.8f}≤C_† ⇒ L-0027; high slab via L-0026. FINITE only."
        ),
        "clay_implication": (
            "None. Shell-support subclass including r=6 only; not all-IC; not Clay."
        ),
        "notes": (
            "Consecutive {1..6} fails Sym Shor (L-0041 techo). "
            "Witness grown from seed {1,2,3,4,6}."
        ),
    }


def save_lemma_l0041(
    bound: GalerkinBoundL0041,
    path: str | Path = "conjectures/proved_restricted/L-0041.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(bound.as_dict(), indent=2), encoding="utf-8")


def save_cr0007(
    d: dict,
    path: str | Path = "conjectures/proved_restricted/C-R-0007.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(d, indent=2), encoding="utf-8")
