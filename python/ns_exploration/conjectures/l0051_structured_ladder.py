"""
L-0051: Structured ladder for C-0007 stretch constant (post SOS scaling techo).

After L-0050 scaled TSSOS through D=356 (band {1..12}), the SOS/Shor ratio
stabilizes at ≈0.56–0.57. Extrapolation to full-dealias Shor C≈25.93 gives
C_ub,SOS ≈14.5 ≫ C_† ≈9.56 — brute-force band SOS cannot close all-IC.

This module packages the three-tier evidence ladder and ranks lightweight
structured routes (no D≳500 SDP on a laptop).

FINITE Galerkin only. Not continuum. Not Clay.
"""

from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

from ns_exploration.conjectures.l0027_low_slab_cubic import lemma_l0027
from ns_exploration.conjectures.l0043_onepol_shellblock import (
    SUPPORT_CR0008,
    lemma_l0043,
    sym_onepol_C_shor,
)
from ns_exploration.conjectures.l0048_fullsym_techo import C_FULLSYM_FULL_DEALIAS

# Frozen C-R-0012 greedy witness (L-0047)
C_GREEDY_41_SHELLS = 9.5037

# Frozen L-0049 rank-1 lower (full dealias)
C_RANK1_FULL = 3.4463930313736206

SCALE_JSON = Path("tools/sos_julia/data/scale_results.json")
ONEPOL_SOS_JSON = Path("tools/sos_julia/data/onepol_cr0008/tssos_result.json")


@dataclass
class BandSosRow:
    radii: list[int]
    D: int
    C_fullsym: float
    C_ub: float
    ratio: float


@dataclass
class StructuredRoute:
    id: str
    target: str
    D_est: int | None
    C_bound_est: float | None
    resource: str
    closes_all_ic: bool
    status: str
    note: str


@dataclass
class GalerkinBoundL0051:
    lemma_id: str = "L-0051"
    route: str = "B"
    status: str = "exploring"
    evidence_level: str = "N2"
    n: int = 24
    C_dagger: float = 0.0
    C_rank1_full: float = C_RANK1_FULL
    C_fullsym_full: float = C_FULLSYM_FULL_DEALIAS
    C_greedy_41_shells: float = C_GREEDY_41_SHELLS
    sos_scale_rows: list[BandSosRow] | None = None
    ratio_median_large: float = 0.0
    C_ub_sos_extrapolated_full: float = 0.0
    brute_sos_closes_all_ic: bool = False
    structured_routes: list[StructuredRoute] | None = None
    clay_implication: str = (
        "None. Synthesis of finite Galerkin bounds; not continuum; not Clay."
    )
    notes: str = ""

    def as_dict(self) -> dict:
        d = asdict(self)
        if self.sos_scale_rows:
            d["sos_scale_rows"] = [asdict(r) for r in self.sos_scale_rows]
        if self.structured_routes:
            d["structured_routes"] = [asdict(r) for r in self.structured_routes]
        return d


def _load_scale_rows(path: Path = SCALE_JSON) -> list[BandSosRow]:
    if not path.is_file():
        return []
    data = json.loads(path.read_text(encoding="utf-8"))
    rows: list[BandSosRow] = []
    for b in data.get("bands", []):
        cu = b.get("C_ub")
        cf = b.get("C_fullsym")
        if cu is None or cf is None:
            continue
        r = b.get("ratio_ub_over_fullsym")
        if r is None and cf > 0:
            r = float(cu) / float(cf)
        rows.append(
            BandSosRow(
                radii=list(b["radii"]),
                D=int(b["D"]),
                C_fullsym=float(cf),
                C_ub=float(cu),
                ratio=float(r),
            )
        )
    return rows


def _median(xs: list[float]) -> float:
    ys = sorted(xs)
    m = len(ys) // 2
    return ys[m] if len(ys) % 2 else 0.5 * (ys[m - 1] + ys[m])


def _load_onepol_sos(path: Path = ONEPOL_SOS_JSON) -> dict | None:
    if not path.is_file():
        return None
    return json.loads(path.read_text(encoding="utf-8"))


def lemma_l0051(n: int = 24, scale_path: Path = SCALE_JSON) -> GalerkinBoundL0051:
    Cd = lemma_l0027(n=n, empirical=False).C_dagger
    rows = _load_scale_rows(scale_path)
    large = [r for r in rows if r.D >= 112]
    ratio_med = _median([r.ratio for r in large]) if large else 0.57
    C_ub_extrap = ratio_med * C_FULLSYM_FULL_DEALIAS
    brute_closes = C_ub_extrap < Cd - 1e-6

    l43 = lemma_l0043(n=n)
    op = sym_onepol_C_shor(n, SUPPORT_CR0008)
    onepol_sos = _load_onepol_sos()
    C_onepol_sos = float(onepol_sos["C_ub"]) if onepol_sos else None
    C_onepol_fullsym = (
        float(onepol_sos["C_fullsym"]) if onepol_sos else float(op["C_shor_sym"])
    )
    onepol_ratio = (
        C_onepol_sos / C_onepol_fullsym
        if C_onepol_sos is not None and C_onepol_fullsym > 0
        else None
    )

    routes = [
        StructuredRoute(
            id="greedy-sos-scan",
            target="greedy C-R-0012 SOS feasibility (D~1772)",
            D_est=1772,
            C_bound_est=C_GREEDY_41_SHELLS * ratio_med,
            resource="cluster est. ~M vars, ~100h+; laptop infeasible",
            closes_all_ic=False,
            status="techo_laptop",
            note=f"L-0055: extrap SOS~{C_GREEDY_41_SHELLS * ratio_med:.2f}<C_dagger but SDP too large for laptop.",
        ),
        StructuredRoute(
            id="greedy-41",
            target="C-R-0012 greedy 41 shells (physical fullsym Shor)",
            D_est=None,
            C_bound_est=C_GREEDY_41_SHELLS,
            resource="Python streaming, ~GB RAM",
            closes_all_ic=False,
            status="proved_restricted",
            note="C≈9.50≤C_†; largest closed Hermitian subclass without SOS.",
        ),
        StructuredRoute(
            id="onepol-cr0008",
            target=f"one-pol Sym shells {list(SUPPORT_CR0008)} (C-R-0008)",
            D_est=int(op["D"]),
            C_bound_est=float(op["C_shor_sym"]),
            resource="Python SVD, seconds",
            closes_all_ic=False,
            status="proved_restricted",
            note="C≈8.48≤C_†; fixed polarization — structured subclass.",
        ),
        StructuredRoute(
            id="rank1-full",
            target="full-dealias rank-1 power (L-0049)",
            D_est=6748,
            C_bound_est=C_RANK1_FULL,
            resource="Python sparse matvec, minutes",
            closes_all_ic=False,
            status="exploring",
            note="Lower bound only; does not refute C-0007.",
        ),
        StructuredRoute(
            id="sos-tiny",
            target="TSSOS order-2 on {1,2,3} Hermitian (L-0050)",
            D_est=52,
            C_bound_est=0.843822489307523,
            resource="Julia Clarabel, ~1 min",
            closes_all_ic=False,
            status="validated",
            note="Beats Shor; N5 Gram cert L-0052 (band {1,2,3}).",
        ),
        StructuredRoute(
            id="sos-onepol-d92",
            target=f"one-pol SOS on C-R-0008 support (D={op['D']})",
            D_est=int(op["D"]),
            C_bound_est=C_onepol_sos,
            resource="Julia Clarabel, ~5 min (302s measured)",
            closes_all_ic=False,
            status="validated" if C_onepol_sos is not None else "planned",
            note=(
                f"SOS C≈{C_onepol_sos:.3f} ≪ Shor≈{C_onepol_fullsym:.2f} "
                f"(ratio≈{onepol_ratio:.3f}); N5 Gram L-0053."
                if C_onepol_sos is not None and onepol_ratio is not None
                else "Structured SOS: half the Hermitian D; may tighten one-pol Shor≈8.48."
            ),
        ),
        StructuredRoute(
            id="sos-brute-scale",
            target="consecutive Hermitian bands to full dealias",
            D_est=6748,
            C_bound_est=C_ub_extrap,
            resource="Julia COSMO, 10+ GB RAM, hours",
            closes_all_ic=False,
            status="techo",
            note=(
                f"Ratio≈{ratio_med:.3f} on D≥112 ⇒ extrapolated C_ub≈{C_ub_extrap:.2f}"
                f" > C_†≈{Cd:.2f}. Abandon laptop brute scaling."
            ),
        ),
    ]

    notes = (
        f"Ladder: C_rank1≈{C_RANK1_FULL:.2f} ≪ C_ub,SOS,band ≪ C_fullsym,full≈"
        f"{C_FULLSYM_FULL_DEALIAS:.2f}. Median SOS/Shor ratio≈{ratio_med:.3f} on "
        f"D≥112 ⇒ extrapolated all-IC SOS≈{C_ub_extrap:.2f} "
        f"{'<' if brute_closes else '>'} C_†≈{Cd:.2f}. "
        f"Pivot: proved subclasses (C-R-0008 one-pol Shor≈{l43.support_cr0008_C:.2f}"
        + (
            f", SOS≈{C_onepol_sos:.2f}" if C_onepol_sos is not None else ""
        )
        + "), targeted low-D SOS, N5 rationalization — not D=500+ brute."
    )
    return GalerkinBoundL0051(
        n=n,
        C_dagger=Cd,
        sos_scale_rows=rows,
        ratio_median_large=ratio_med,
        C_ub_sos_extrapolated_full=C_ub_extrap,
        brute_sos_closes_all_ic=brute_closes,
        structured_routes=routes,
        notes=notes,
    )


def save_lemma_l0051(
    bound: GalerkinBoundL0051,
    path: str | Path = "conjectures/active/L-0051.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(bound.as_dict(), indent=2), encoding="utf-8")
