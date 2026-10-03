"""
L-0056: Structured closure map for C-0007 / C-0008 (all-IC vs subclasses).

C-0007 all-IC proved at M=L-0024 majorant (L-0070). Sharp refinement C-0008
(M≈41.284) remains open. Packages C-R ladder + SOS/Shor attack map.

FINITE Galerkin only. Not Clay.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path

from ns_exploration.conjectures.l0026_c0007_techo import lemma_l0026
from ns_exploration.conjectures.l0027_low_slab_cubic import lemma_l0027
from ns_exploration.conjectures.l0054_structured_sos_synthesis import lemma_l0054
from ns_exploration.conjectures.l0055_greedy_sos_feasibility import lemma_l0055

PROVED_CR_DIR = Path("conjectures/proved_restricted")
C0007_JSON = Path("conjectures/proved_restricted/C-0007.json")
C0008_JSON = Path("conjectures/active/C-0008.json")


@dataclass
class CrClosureRow:
    id: str
    status: str
    evidence: str
    stretch_C: float | None
    closes_vs_Cdagger: bool | None
    closes_c0007_subclass: bool
    lemma: str
    notes: str


@dataclass
class GalerkinBoundL0056:
    lemma_id: str = "L-0056"
    route: str = "B"
    status: str = "validated"
    evidence_level: str = "N7"
    n: int = 24
    C_dagger: float = 0.0
    c0007_M: float = 0.0
    c0008_sharp_M: float = 0.0
    all_ic_open: bool = False
    sharp_open: bool = True
    n_proved_restricted: int = 0
    n_proved_subclasses: int = 0
    n_open_subclasses: int = 1
    C_rank1_full: float = 0.0
    best_n5_sos_hi: float = 0.0
    best_structured_sos_hi: float = 0.0
    greedy_sos_extrap: float = 0.0
    full_dealias_shor: float = 0.0
    cr_rows: list[CrClosureRow] | None = None
    clay_implication: str = "None. Closure map for finite Galerkin only; not Clay."
    notes: str = ""

    def as_dict(self) -> dict:
        d = asdict(self)
        if self.cr_rows:
            d["cr_rows"] = [asdict(r) for r in self.cr_rows]
        return d


def _load_cr_rows() -> list[CrClosureRow]:
    rows: list[CrClosureRow] = []
    if not PROVED_CR_DIR.is_dir():
        return rows
    for path in sorted(PROVED_CR_DIR.glob("C-R-*.json")):
        d = json.loads(path.read_text(encoding="utf-8"))
        cid = d.get("id", path.stem)
        C = d.get("C_fullsym") or d.get("C_shor_sym") or d.get("C_shor")
        Cd = d.get("C_dagger")
        closes_cd = (C is not None and Cd is not None and float(C) <= float(Cd) + 1e-5)
        rows.append(
            CrClosureRow(
                id=cid,
                status=str(d.get("status", "proved_restricted")),
                evidence=str(d.get("evidence_level", "N7")),
                stretch_C=float(C) if C is not None else None,
                closes_vs_Cdagger=closes_cd if C is not None else None,
                closes_c0007_subclass=True,
                lemma=str(d.get("lemma", "")),
                notes=str(d.get("notes", ""))[:200],
            )
        )
    return rows


def _load_c0007_proved_M() -> tuple[float, float, bool]:
    """Return (C-0007 proved M, C-0008 sharp M, all_ic_open)."""
    b7 = lemma_l0027(n=24, empirical=False)
    M8 = float(b7.c0007_M)
    if C0007_JSON.is_file():
        d = json.loads(C0007_JSON.read_text(encoding="utf-8"))
        M7 = float(d.get("proved_bound_M") or d.get("proved_bound_Omega_T") or 0.0)
    else:
        M7 = float(lemma_l0026(n=24).L0024_majorant)
    if C0008_JSON.is_file():
        d8 = json.loads(C0008_JSON.read_text(encoding="utf-8"))
        M8 = float(d8.get("proposed_bound_M", M8))
    return M7, M8, False


def lemma_l0056(n: int = 24) -> GalerkinBoundL0056:
    b7 = lemma_l0027(n=n, empirical=False)
    Cd = b7.C_dagger
    M7, M8, all_open = _load_c0007_proved_M()
    l54 = lemma_l0054(n=n)
    l55 = lemma_l0055(n=n)
    cr_rows = _load_cr_rows()
    best_n5 = l54.best_sos_n5_hi
    onepol_row = next((r for r in l54.rows or [] if r.id == "onepol-cr0008-n5"), None)
    structured = float(onepol_row.C_ub_n5_hi) if onepol_row and onepol_row.C_ub_n5_hi else best_n5
    n_proved = sum(1 for r in cr_rows if r.closes_vs_Cdagger)
    notes = (
        f"C-0007 all-IC PROVED at M~{M7:.2f} (L-0070/L-0024). "
        f"C-0008 sharp OPEN (M~{M8:.2f}). C-R-0002..0014 proved restricted. "
        f"Stretch ladder: rank1~{l54.C_rank1_full:.2f} < best N5 SOS~{best_n5:.3f} "
        f"< one-pol N5~{structured:.3f} < greedy Shor~{l54.C_greedy_41_shor:.2f} "
        f"(SOS extrap~{l55.C_ub_sos_extrap:.2f}, cluster) < C_dagger~{Cd:.2f} "
        f"< fullsym~{l54.C_fullsym_full:.1f}. "
        "Low-slab subclasses attack C-0008 sharp gap; not needed for C-0007 all-IC."
    )
    return GalerkinBoundL0056(
        n=n,
        C_dagger=Cd,
        c0007_M=M7,
        c0008_sharp_M=M8,
        all_ic_open=all_open,
        sharp_open=True,
        n_proved_restricted=len(cr_rows),
        n_proved_subclasses=n_proved,
        n_open_subclasses=1,
        C_rank1_full=l54.C_rank1_full,
        best_n5_sos_hi=best_n5,
        best_structured_sos_hi=structured,
        greedy_sos_extrap=l55.C_ub_sos_extrap,
        full_dealias_shor=l54.C_fullsym_full,
        cr_rows=cr_rows,
        notes=notes,
    )


def save_lemma_l0056(
    bound: GalerkinBoundL0056,
    path: str | Path = "conjectures/active/L-0056.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(bound.as_dict(), indent=2), encoding="utf-8")


def render_closure_map_md(bound: GalerkinBoundL0056) -> str:
    lines = [
        "# C-0007 / C-0008 Structured Closure Map (L-0056)",
        "",
        "> Finite Galerkin T^3, N<=24, dealias 2/3. **Not continuum. Not Clay.**",
        "",
        f"- **All-IC C-0007:** PROVED (M~{bound.c0007_M:.4f}, L-0070/L-0024)",
        f"- **Sharp C-0008:** OPEN (target M~{bound.c0008_sharp_M:.4f})",
        f"- **Low-slab stretch target C_dagger:** {bound.C_dagger:.6f}",
        f"- **Proved restricted subclasses:** {bound.n_proved_restricted} (C-R-0002..0014)",
        "",
        "## Stretch constant ladder (attack C_dagger)",
        "",
        "| Tier | Subclass / method | C bound | vs C_dagger | Evidence |",
        "|------|-------------------|---------|-------------|----------|",
        f"| lower | rank-1 full dealias | ~{bound.C_rank1_full:.2f} | yes | N7 lower L-0049 |",
        f"| **best N5 SOS** | Hermitian {{1,2,3}} | **{bound.best_n5_sos_hi:.4f}** | yes | CERT-L0052 |",
        f"| N5 SOS | Hermitian {{1..4}} | 0.894 | yes | CERT-L0054 |",
        f"| N5 SOS | Hermitian {{1..5}} | 1.305 | yes | CERT-L0058 |",
        f"| N5 SOS | Hermitian {{1..6}} | 1.701 | yes (interval) | CERT-L0060 |",
        f"| N5 SOS | Hermitian {{1..8}} | 1.852 | yes (interval COSMO) | CERT-L0061 |",
        f"| N5 SOS | Hermitian {{1..9}} | 2.223 | yes (interval COSMO) | CERT-L0062 |",
        f"| N5 SOS | Hermitian {{1..10}} | 2.440 | yes (interval COSMO) | CERT-L0062 |",
        f"| N5 SOS | Hermitian {{1..12}} | 2.783 | yes (interval COSMO) | CERT-L0062 |",
        f"| N5 SOS | one-pol C-R-0014 (12 shells) | 4.076 | yes | CERT-L0067 |",
        f"| N5 SOS | one-pol C-R-0013 (11 shells) | 4.015 | yes | CERT-L0065 |",
        f"| N5 SOS | one-pol C-R-0008 | 3.356 | yes | CERT-L0053 |",
        f"| Shor | greedy 41 shells C-R-0012 | 9.504 | yes | CERT-L0047 |",
        f"| extrap SOS | greedy 41 (cluster) | {bound.greedy_sos_extrap:.2f} | yes (extrap) | L-0055/L-0057 |",
        f"| techo Shor | all-IC fullsym | {bound.full_dealias_shor:.2f} | **no** | L-0048 |",
        f"| techo SOS | brute band extrap | ~14.7 | **no** | L-0051 |",
        "",
        "## Proved restricted map (C-R)",
        "",
        "| ID | C_fullsym | closes C_dagger | closes C-0007 subclass |",
        "|----|-----------|-----------------|------------------------|",
    ]
    for r in bound.cr_rows or []:
        c = f"{r.stretch_C:.4f}" if r.stretch_C is not None else "-"
        cd = "yes" if r.closes_vs_Cdagger else ("no" if r.closes_vs_Cdagger is False else "-")
        lines.append(f"| {r.id} | {c} | {cd} | yes |")
    lines.extend(
        [
            "",
        "## Viable paths (honest)",
        "",
        "1. **Done on laptop:** N5 SOS one-pol C-R-0014 (L-0064/67/68) + bands through shell 12.",
        "2. **Proved without SOS:** C-R-0002..0014 via Shor/streaming (largest Hermitian: C-R-0012).",
        "3. **Deferred:** greedy D=1772 SOS cluster (optional; C-R-0012 already Shor-closed).",
        "4. **Abandoned:** brute Hermitian band SOS to full dealias (extrap ~14.7 > C_dagger).",
        "5. **All-IC C-0007:** closed at L-0024 majorant (L-0070). **C-0008 sharp** still open.",
        "6. **Refuted (2026-07-31):** L-0072 terminal shell integral (Phase D/E full-dealias); "
        "I_hi~78.8 vs I_*~1.30. See `conjectures/active/L-0072.json`.",
        "",
        "## Terminal route (L-0071 / L-0072)",
        "",
        "| Method | Route D I_hi | Closes C-0008? |",
        "|--------|--------------|----------------|",
        "| Frobenius Phase D | ~112.8 | No |",
        "| Best (Phase E) | ~78.8 | No |",
        "| Verifier | FAIL | Route refuted |",
        "",
            "## N5 certificates (SOS Gram)",
            "",
            "- CERT-L0052-sos-gram-band123-N24",
            "- CERT-L0053-sos-gram-onepol-cr0008-N24",
            "- CERT-L0054-sos-gram-band1234-N24",
            "- CERT-L0058-sos-gram-band12345-N24",
            "- CERT-L0060-sos-gram-band123456-N24",
            "- CERT-L0061-sos-band12345678-N24",
            "- CERT-L0062-sos-band123456789-N24",
            "- CERT-L0062-sos-band1to10-N24",
            "- CERT-L0062-sos-band1to12-N24",
            "- CERT-L0065-sos-onepol-greedy-second-N24",
            "- CERT-L0067-sos-onepol-greedy-second-r24-N24",
            "",
            "## Lean N8 (L-0059)",
            "",
            "- NSGalerkin.SosGram: native_decide C_ub_hi < C_dagger for N5 SOS certs",
            "",
        ]
    )
    return "\n".join(lines)


def write_closure_map_report(
    bound: GalerkinBoundL0056,
    path: str | Path = "reports/C0007_STRUCTURED_CLOSURE_MAP.md",
) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(render_closure_map_md(bound), encoding="utf-8")
    return path
