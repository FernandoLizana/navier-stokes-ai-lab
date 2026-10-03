"""
L-0018: Spectral enstrophy envelope closes C-S-0002 on N≤16.

On the 2/3-dealiased Galerkin space of size N, every retained mode satisfies
|k|² ≤ K²(N) := max{|k|² : k in dealias mask}. With the project conventions
  E = (1/2) Σ_k ‖û_k‖²,   Ω = (1/2) Σ_k |k|² ‖û_k‖²
(div-free: Ω = (1/2)‖ω‖₂²), one has the Parseval envelope
  Ω ≤ K²(N) E.
For unforced Navier–Stokes Galerkin, E(t) ≤ E(0) = E0, hence
  Ω(t) ≤ K²(N) E0   for all t ≥ 0.

For N≤16: K²(16)=75 is maximal among N∈{8,…,16} on the dealias mask, so
  Ω(t) ≤ 75 · E0 = 37.5 < M_CS0002 ≈ 48.159.
This proves C-S-0002 (shell IC |k|_∞≤4 is a subclass) at evidence N7,
with K² certified by exact integer max over the mask (N5).

Remark. The L-0012–L-0017 comparison bounds used Ω ≤ K_IC² E0 + K_full² E_H,
which allows E_L=E0 and E_H>0 at once (valid but loose). The envelope avoids
that double-budget and is strong enough to close the parent at T=0.02 for all t.

FINITE Galerkin / dealias only. Not continuum. Not Clay. Evidence: N7 (+ N5).
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path

from ns_exploration.validation.l0003_certificate import full_dealias_exact_stats


def dealias_K2(n: int) -> int:
    """Exact max |k|² on the 2/3 dealias mask."""
    K2, _ = full_dealias_exact_stats(n)
    return int(K2)


def envelope_Omega(n: int, E0: float) -> float:
    """Uniform-in-time enstrophy cap Ω ≤ K²(N) E0."""
    return float(dealias_K2(n)) * float(E0)


def max_K2_upto(n_max: int) -> tuple[int, int]:
    """Return (worst K², argmax N) over even N in [8, n_max] (project grid sizes)."""
    best_K = -1
    best_n = n_max
    for n in range(8, n_max + 1, 2):
        K2 = dealias_K2(n)
        if K2 >= best_K:
            best_K = K2
            best_n = n
    return int(best_K), int(best_n)


@dataclass
class GalerkinBoundL0018:
    lemma_id: str = "L-0018"
    route: str = "B"
    status: str = "proved_restricted"
    evidence_level: str = "N7"
    n_max: int = 16
    K2_worst: int = 0
    n_worst: int = 0
    E0: float = 0.5
    nu: float = 0.1
    Omega_cap: float = 0.0
    cs0002_M: float = 48.15928260103266
    closes_cs0002: bool = False
    margin: float = 0.0
    clay_implication: str = (
        "None. Finite dealiased Galerkin spectral envelope only; "
        "does not address continuum Clay routes."
    )
    notes: str = ""

    def as_dict(self) -> dict:
        return asdict(self)


def lemma_l0018(
    n_max: int = 16,
    E0: float = 0.5,
    nu: float = 0.1,
    cs0002_M: float = 48.15928260103266,
) -> GalerkinBoundL0018:
    K2, n_w = max_K2_upto(n_max)
    cap = float(K2) * E0
    closes = cap <= cs0002_M + 1e-15
    notes = (
        f"Ω ≤ K² E0 with K²({n_w})={K2} ⇒ Ω_cap={cap:.6g} "
        f"{'≤' if closes else '>'} M={cs0002_M:.6g} (margin={cs0002_M - cap:.6g}). "
        f"Closes C-S-0002 for all t≥0 on dealiased Galerkin N≤{n_max}."
    )
    return GalerkinBoundL0018(
        n_max=n_max,
        K2_worst=K2,
        n_worst=n_w,
        E0=E0,
        nu=nu,
        Omega_cap=cap,
        cs0002_M=cs0002_M,
        closes_cs0002=closes,
        margin=cs0002_M - cap,
        notes=notes,
    )


def save_lemma_l0018(
    bound: GalerkinBoundL0018,
    path: str | Path = "conjectures/proved_restricted/L-0018.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(bound.as_dict(), indent=2), encoding="utf-8")
