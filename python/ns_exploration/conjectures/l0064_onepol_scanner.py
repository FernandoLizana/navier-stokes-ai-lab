"""
L-0064: Fast one-pol Shor scanner for structured C-0007 subclasses.

Scans consecutive bands, C-R-0008 extensions, and greedy one-pol growth
using sym_onepol_C_shor (seconds per candidate on laptop). Flags SOS-feasible
candidates (D≤160 Clarabel, D≤356 COSMO) without launching SDP.

FINITE Galerkin only. Not all-IC. Not Clay.
"""

from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

from ns_exploration.conjectures.l0027_low_slab_cubic import lemma_l0027
from ns_exploration.conjectures.l0043_onepol_shellblock import (
    SUPPORT_CR0008,
    sym_onepol_C_shor,
)
from ns_exploration.conjectures.l0048_matrixfree_fullsym import all_dealias_radii

# Measured one-pol SOS / Shor on C-R-0008 (L-0053).
ONEPOL_SOS_RATIO = 3.356 / 8.4764
LAPTOP_D_CLARABEL = 160
LAPTOP_D_COSMO = 356
DEFAULT_MAX_CONSECUTIVE = 15
DEFAULT_MAX_EXTEND_SHELL = 20
DEFAULT_MAX_GREEDY_SHELL = 20


@dataclass
class OnepolScanRow:
    scan_id: str
    pattern: str
    radii: list[int]
    branch: str
    D: int
    C_shor_sym: float
    closes_C_dagger: bool
    margin: float
    est_C_ub_sos: float
    est_closes_sos: bool
    sos_tier: str
    notes: str = ""


@dataclass
class GalerkinBoundL0064:
    lemma_id: str = "L-0064"
    route: str = "B"
    status: str = "exploring"
    evidence_level: str = "N7"
    n: int = 24
    C_dagger: float = 0.0
    onepol_sos_ratio: float = ONEPOL_SOS_RATIO
    n_candidates: int = 0
    n_closing: int = 0
    best_closing: OnepolScanRow | None = None
    greedy_growth_first: list[int] | None = None
    greedy_growth_second: list[int] | None = None
    rows: list[OnepolScanRow] | None = None
    sos_queue: list[OnepolScanRow] | None = None
    clay_implication: str = (
        "None. One-pol Shor scan for structured subclasses; not all-IC; not Clay."
    )
    notes: str = ""

    def as_dict(self) -> dict:
        d = asdict(self)
        if self.rows:
            d["rows"] = [asdict(r) for r in self.rows]
        if self.best_closing:
            d["best_closing"] = asdict(self.best_closing)
        if self.sos_queue:
            d["sos_queue"] = [asdict(r) for r in self.sos_queue]
        return d


def _sos_tier(D: int) -> str:
    if D <= LAPTOP_D_CLARABEL:
        return "clarabel_laptop"
    if D <= LAPTOP_D_COSMO:
        return "cosmo_laptop"
    return "cluster_or_skip"


def _eval_row(
    n: int,
    radii: tuple[int, ...] | list[int],
    branch: str,
    scan_id: str,
    pattern: str,
    Cd: float,
    ratio: float,
    cache: dict[tuple[tuple[int, ...], str], dict] | None = None,
    notes: str = "",
) -> OnepolScanRow:
    key = (tuple(int(r) for r in radii), branch)
    if cache is not None and key in cache:
        ref = cache[key]
    else:
        ref = sym_onepol_C_shor(n, radii, branch=branch)
        if cache is not None:
            cache[key] = ref
    C = float(ref["C_shor_sym"])
    rad = list(ref["radii"])
    D = int(ref["D"])
    closes = C <= Cd + 1e-9
    est_ub = ratio * C
    return OnepolScanRow(
        scan_id=scan_id,
        pattern=pattern,
        radii=rad,
        branch=branch,
        D=D,
        C_shor_sym=C,
        closes_C_dagger=closes,
        margin=Cd - C,
        est_C_ub_sos=est_ub,
        est_closes_sos=est_ub < Cd - 1e-6,
        sos_tier=_sos_tier(D),
        notes=notes,
    )


def _scan_consecutive(
    n: int,
    Cd: float,
    ratio: float,
    max_k: int,
    cache: dict[tuple[tuple[int, ...], str], dict],
) -> list[OnepolScanRow]:
    rows: list[OnepolScanRow] = []
    for k in range(1, max_k + 1):
        radii = tuple(range(1, k + 1))
        for branch in ("first", "second"):
            rows.append(
                _eval_row(
                    n,
                    radii,
                    branch,
                    scan_id=f"consec-1to{k}-{branch}",
                    pattern=f"consecutive {{1..{k}}}",
                    Cd=Cd,
                    ratio=ratio,
                    cache=cache,
                )
            )
    return rows


def _scan_extend_base(
    n: int,
    base: tuple[int, ...],
    Cd: float,
    ratio: float,
    base_id: str,
    max_shell: int,
    cache: dict[tuple[tuple[int, ...], str], dict],
) -> list[OnepolScanRow]:
    rows: list[OnepolScanRow] = []
    base_set = set(base)
    for shell in all_dealias_radii(n):
        if shell > max_shell or shell in base_set:
            continue
        radii = tuple(sorted(base_set | {shell}))
        for branch in ("first", "second"):
            rows.append(
                _eval_row(
                    n,
                    radii,
                    branch,
                    scan_id=f"{base_id}+s{shell}-{branch}",
                    pattern=f"{base_id} + shell {shell}",
                    Cd=Cd,
                    ratio=ratio,
                    cache=cache,
                )
            )
    return rows


def _greedy_growth(
    n: int,
    branch: str,
    Cd: float,
    max_shell: int,
    cache: dict[tuple[tuple[int, ...], str], dict],
) -> list[int]:
    support: set[int] = set()
    for shell in all_dealias_radii(n):
        if shell > max_shell:
            continue
        trial = tuple(sorted(support | {shell}))
        key = (trial, branch)
        if key in cache:
            ref = cache[key]
        else:
            ref = sym_onepol_C_shor(n, trial, branch=branch)
            cache[key] = ref
        if float(ref["C_shor_sym"]) <= Cd + 1e-9:
            support.add(shell)
    return sorted(support)


def _pick_best_closing(rows: list[OnepolScanRow]) -> OnepolScanRow | None:
    closing = [r for r in rows if r.closes_C_dagger]
    if not closing:
        return None
    return max(closing, key=lambda r: (r.D, r.margin, -len(r.radii)))


def _sos_queue(rows: list[OnepolScanRow], known: set[tuple[int, ...]]) -> list[OnepolScanRow]:
    """New closing candidates worth SOS on laptop, largest D first."""
    out: list[OnepolScanRow] = []
    seen: set[tuple[int, ...]] = set()
    for r in sorted(
        [x for x in rows if x.closes_C_dagger and x.est_closes_sos],
        key=lambda x: (-x.D, x.branch),
    ):
        key = (tuple(r.radii), r.branch)
        if key in known or key in seen:
            continue
        if r.sos_tier in ("clarabel_laptop", "cosmo_laptop"):
            out.append(r)
            seen.add(key)
    return out


def lemma_l0064(
    n: int = 24,
    max_consecutive: int = DEFAULT_MAX_CONSECUTIVE,
    max_extend_shell: int = DEFAULT_MAX_EXTEND_SHELL,
    max_greedy_shell: int = DEFAULT_MAX_GREEDY_SHELL,
    onepol_sos_ratio: float = ONEPOL_SOS_RATIO,
) -> GalerkinBoundL0064:
    Cd = lemma_l0027(n=n, empirical=False).C_dagger
    cache: dict[tuple[tuple[int, ...], str], dict] = {}
    rows: list[OnepolScanRow] = []
    rows.extend(_scan_consecutive(n, Cd, onepol_sos_ratio, max_consecutive, cache))
    rows.extend(
        _scan_extend_base(
            n, SUPPORT_CR0008, Cd, onepol_sos_ratio, "cr0008", max_extend_shell, cache
        )
    )
    g_first = _greedy_growth(n, "first", Cd, max_greedy_shell, cache)
    g_second = _greedy_growth(n, "second", Cd, max_greedy_shell, cache)
    for branch, support in (("first", g_first), ("second", g_second)):
        rows.append(
            _eval_row(
                n,
                support,
                branch,
                scan_id=f"greedy-{branch}",
                pattern="greedy one-pol growth",
                Cd=Cd,
                ratio=onepol_sos_ratio,
                cache=cache,
                notes=f"{len(support)} shells (r≤{max_greedy_shell})",
            )
        )
    known = {tuple(SUPPORT_CR0008)}
    best = _pick_best_closing(rows)
    queue = _sos_queue(rows, known)
    closing = [r for r in rows if r.closes_C_dagger]
    notes = (
        f"Scanned {len(rows)} one-pol candidates (consecutive 1..{max_consecutive}, "
        f"cr0008 extensions, greedy growth). {len(closing)} close C_†≈{Cd:.3f}. "
        f"Best D={best.D if best else 0} C≈{best.C_shor_sym:.4f} ({best.pattern if best else 'n/a'}). "
        f"Greedy first: {len(g_first)} shells; second: {len(g_second)}. "
        f"SOS laptop queue: {len(queue)} new candidates."
    )
    return GalerkinBoundL0064(
        n=n,
        C_dagger=Cd,
        onepol_sos_ratio=onepol_sos_ratio,
        n_candidates=len(rows),
        n_closing=len(closing),
        best_closing=best,
        greedy_growth_first=g_first,
        greedy_growth_second=g_second,
        rows=rows,
        sos_queue=queue,
        notes=notes,
    )


def markdown_report(bound: GalerkinBoundL0064) -> str:
    lines = [
        "# One-pol structured scan (L-0064)",
        "",
        "> Finite Galerkin T³, N≤24. **Not all-IC. Not Clay.**",
        "",
        f"- **C_dagger:** {bound.C_dagger:.6f}",
        f"- **One-pol SOS ratio (L-0053):** {bound.onepol_sos_ratio:.4f}",
        f"- **Candidates:** {bound.n_candidates} | **Closing:** {bound.n_closing}",
        "",
    ]
    if bound.best_closing:
        b = bound.best_closing
        lines.extend(
            [
                "## Best closing witness",
                "",
                f"| Field | Value |",
                f"|-------|-------|",
                f"| Pattern | {b.pattern} |",
                f"| Branch | {b.branch} |",
                f"| Shells | `{b.radii}` |",
                f"| D | {b.D} |",
                f"| C_Shor | {b.C_shor_sym:.6f} |",
                f"| Est SOS C_ub | {b.est_C_ub_sos:.4f} |",
                f"| SOS tier | {b.sos_tier} |",
                "",
            ]
        )
    lines.extend(
        [
            "## Greedy one-pol growth",
            "",
            f"- **first:** `{bound.greedy_growth_first}` ({len(bound.greedy_growth_first or [])} shells)",
            f"- **second:** `{bound.greedy_growth_second}` ({len(bound.greedy_growth_second or [])} shells)",
            "",
        ]
    )
    if bound.sos_queue:
        lines.extend(["## SOS laptop queue (new, D≤356)", ""])
        lines.append("| scan_id | D | C_Shor | est SOS | tier | radii |")
        lines.append("|---------|---|--------|---------|------|-------|")
        for r in bound.sos_queue[:15]:
            rad = ",".join(str(x) for x in r.radii)
            lines.append(
                f"| {r.scan_id} | {r.D} | {r.C_shor_sym:.4f} | {r.est_C_ub_sos:.3f} | "
                f"{r.sos_tier} | `{rad}` |"
            )
        lines.append("")
    lines.extend(["## Consecutive ladder (first branch)", ""])
    lines.append("| k | D | C_Shor | closes | est SOS |")
    lines.append("|---|-----|--------|--------|---------|")
    for r in bound.rows or []:
        if r.pattern.startswith("consecutive") and r.branch == "first":
            k = r.radii[-1] if r.radii else 0
            lines.append(
                f"| {k} | {r.D} | {r.C_shor_sym:.4f} | {'yes' if r.closes_C_dagger else 'no'} | "
                f"{r.est_C_ub_sos:.3f} |"
            )
    lines.append("")
    lines.extend(["## Next SOS commands (laptop)", ""])
    lines.append("```powershell")
    lines.append("$env:PYTHONPATH='python'")
    for r in (bound.sos_queue or [])[:3]:
        rad = " ".join(str(x) for x in r.radii)
        slug = r.scan_id.replace("+", "_").replace("-", "_")
        lines.append(
            f"python -m ns_exploration.experiments.export_sos_cubic --n 24 --onepol "
            f"--branch {r.branch} --radii {rad} --out tools/sos_julia/data/onepol_{slug}"
        )
    lines.append("```")
    lines.append("")
    return "\n".join(lines)


def save_lemma_l0064(
    bound: GalerkinBoundL0064,
    path: str | Path = "conjectures/active/L-0064.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(bound.as_dict(), indent=2), encoding="utf-8")


def save_report(
    bound: GalerkinBoundL0064,
    path: str | Path = "reports/ONEPOL_SCAN_L0064.md",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(markdown_report(bound), encoding="utf-8")
