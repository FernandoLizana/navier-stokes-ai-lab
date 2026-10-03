"""
L-0055: Greedy C-R-0012 SOS feasibility scan (D≈1772).

Estimates SDP scale / RAM for order-2 TSSOS on the 41-shell greedy witness
without launching the solve. Uses calibrated laptop runs (D≤356) + basis counts.

FINITE Galerkin only. Not Clay.
"""

from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from ns_exploration.conjectures.l0027_low_slab_cubic import lemma_l0027
from ns_exploration.conjectures.l0047_greedy_fullsym import SUPPORT_CR0012
from ns_exploration.conjectures.l0051_structured_ladder import C_GREEDY_41_SHELLS
from ns_exploration.experiments.sprintB_physical_tensor import build_hermitian_basis

# Calibrated TSSOS+Clarabel/COSMO (order=2, TS=MD): (D, sdp_variables, solve_seconds)
SDP_CALIB = (
    (52, 5803, 20.0),
    (64, 8203, 68.0),
    (92, 18972, 552.0),
    (356, 305_000, 4521.0),  # vars: power extrap from D=52; time: measured COSMO
)
LAPTOP_D_CLARABEL = 160
LAPTOP_D_COSMO = 356
RAM_DENSE_M_CAP_GB = 12.0  # L-0047 float32 M cap


@dataclass
class PrefixRow:
    n_shells: int
    radii: list[int]
    D: int
    est_sdp_vars: float
    est_solve_s: float
    laptop_clarabel_ok: bool
    laptop_cosmo_ok: bool


@dataclass
class GalerkinBoundL0055:
    lemma_id: str = "L-0055"
    route: str = "B"
    status: str = "techo_laptop"
    evidence_level: str = "N7"
    n: int = 24
    C_dagger: float = 0.0
    C_greedy_shor: float = C_GREEDY_41_SHELLS
    greedy_D: int = 0
    n_greedy_shells: int = 41
    ratio_median_sos: float = 0.567
    C_ub_sos_extrap: float = 0.0
    extrap_closes_cr0012: bool = False
    est_sdp_variables: float = 0.0
    est_solve_hours: float = 0.0
    est_dense_A_gb: float = 0.0
    est_gram_export_gb: float = 0.0
    laptop_feasible: bool = False
    verdict: str = ""
    prefix_rows: list[PrefixRow] | None = None
    clay_implication: str = (
        "None. Feasibility scan only; not a SOS certificate; not Clay."
    )
    notes: str = ""

    def as_dict(self) -> dict:
        d = asdict(self)
        if self.prefix_rows:
            d["prefix_rows"] = [asdict(r) for r in self.prefix_rows]
        return d


def _fit_power_law(
    points: tuple[tuple[float, float], ...] | list[tuple[float, float]],
) -> tuple[float, float]:
    """Return (a, b) with y ≈ exp(a + b*log(x))."""
    pts = list(points)
    xs = np.array([math.log(p[0]) for p in pts], dtype=float)
    ys = np.array([math.log(p[1]) for p in pts], dtype=float)
    b, a = np.polyfit(xs, ys, 1)
    return float(a), float(b)


def _est_at_D(D: int, a_vars: float, b_vars: float, a_time: float, b_time: float) -> tuple[float, float]:
    ld = math.log(max(D, 1))
    vars_est = math.exp(a_vars + b_vars * ld)
    time_est = math.exp(a_time + b_time * ld)
    return vars_est, time_est


def _est_nnz_A(D: int) -> float:
    """Heuristic nnz(fullsym A) from band calibrations."""
    # band_123: 2208/52, band_1234: 3072/64 → ~42*D
    return 42.0 * float(D)


def _prefix_scan(n: int = 24) -> list[PrefixRow]:
    a_v, b_v = _fit_power_law(((p[0], p[1]) for p in SDP_CALIB))
    a_t, b_t = _fit_power_law(((p[0], p[2]) for p in SDP_CALIB))
    rows: list[PrefixRow] = []
    shells = list(SUPPORT_CR0012)
    for k in range(1, len(shells) + 1):
        radii = tuple(shells[:k])
        D = build_hermitian_basis(n, radii).D
        ev, et = _est_at_D(D, a_v, b_v, a_t, b_t)
        rows.append(
            PrefixRow(
                n_shells=k,
                radii=list(radii),
                D=D,
                est_sdp_vars=ev,
                est_solve_s=et,
                laptop_clarabel_ok=D <= LAPTOP_D_CLARABEL,
                laptop_cosmo_ok=D <= LAPTOP_D_COSMO,
            )
        )
    return rows


def lemma_l0055(n: int = 24, ratio_median: float = 0.567) -> GalerkinBoundL0055:
    Cd = lemma_l0027(n=n, empirical=False).C_dagger
    rows = _prefix_scan(n=n)
    D_g = build_hermitian_basis(n, SUPPORT_CR0012).D
    a_v, b_v = _fit_power_law(((p[0], p[1]) for p in SDP_CALIB))
    a_t, b_t = _fit_power_law(((p[0], p[2]) for p in SDP_CALIB))
    est_vars, est_time = _est_at_D(D_g, a_v, b_v, a_t, b_t)
    nnz = _est_nnz_A(D_g)
    dense_A_gb = nnz * 24 / (1024**3)  # i,j,k int32 + float64
    # band_123 gram json ~110k lines for D=52; scale ~ (D/52)^3 rough
    gram_gb = 0.11 * (D_g / 52.0) ** 3
    C_ub_extrap = ratio_median * C_GREEDY_41_SHELLS
    closes = C_ub_extrap < Cd - 1e-6
    laptop_ok = D_g <= LAPTOP_D_COSMO and est_vars < 5e6 and gram_gb < 2.0
    if laptop_ok:
        verdict = "marginal_laptop"
    elif closes and est_time < 86400 * 7:
        verdict = "cluster_candidate"
    else:
        verdict = "laptop_infeasible"
    max_cosmo_prefix = max((r for r in rows if r.laptop_cosmo_ok), key=lambda r: r.D, default=None)
    notes = (
        f"Greedy C-R-0012: D={D_g}, C_Shor≈{C_GREEDY_41_SHELLS:.4f}. "
        f"Extrap SOS C_ub≈{C_ub_extrap:.3f} "
        f"({'<' if closes else '>'} C_†≈{Cd:.3f}). "
        f"Est SDP vars≈{est_vars/1e6:.1f}M, solve≈{est_time/3600:.1f}h, "
        f"Gram export≈{gram_gb:.1f}GB. Verdict: {verdict}. "
        f"Largest COSMO-safe prefix: {max_cosmo_prefix.n_shells if max_cosmo_prefix else 0} shells "
        f"(D={max_cosmo_prefix.D if max_cosmo_prefix else 0})."
    )
    return GalerkinBoundL0055(
        n=n,
        C_dagger=Cd,
        greedy_D=D_g,
        ratio_median_sos=ratio_median,
        C_ub_sos_extrap=C_ub_extrap,
        extrap_closes_cr0012=closes,
        est_sdp_variables=est_vars,
        est_solve_hours=est_time / 3600.0,
        est_dense_A_gb=dense_A_gb,
        est_gram_export_gb=gram_gb,
        laptop_feasible=laptop_ok,
        verdict=verdict,
        prefix_rows=rows,
        notes=notes,
    )


def save_lemma_l0055(
    bound: GalerkinBoundL0055,
    path: str | Path = "conjectures/active/L-0055.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(bound.as_dict(), indent=2), encoding="utf-8")
