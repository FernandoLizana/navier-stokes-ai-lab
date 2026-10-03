"""
L-0054: Structured SOS evidence synthesis (post L-0052/L-0053 N5 Gram certs).

Packages validated structured routes: Shor subclasses, N2 TSSOS bands,
N5 rational Gram certificates. All-IC remains open; brute band SOS is techo.

FINITE Galerkin only. Not continuum. Not Clay.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path

from ns_exploration.conjectures.l0027_low_slab_cubic import lemma_l0027
from ns_exploration.conjectures.l0043_onepol_shellblock import SUPPORT_CR0008
from ns_exploration.conjectures.l0051_structured_ladder import C_GREEDY_41_SHELLS, C_RANK1_FULL
from ns_exploration.conjectures.l0052_sos_gram_n5 import (
    C_FULLSYM_BAND_123,
    C_UB_BAND_123,
    certified_C_ub_hi,
    load_gram_export,
    verify_gram_blocks,
)
from ns_exploration.conjectures.l0048_fullsym_techo import C_FULLSYM_FULL_DEALIAS
from ns_exploration.conjectures.l0053_onepol_sos_gram_n5 import (
    C_FULLSYM_ONEPOL,
    C_UB_ONEPOL,
)

SCALE_JSON = Path("tools/sos_julia/data/scale_results.json")
CERT_L0052 = Path("certificates/CERT-L0052-sos-gram-band123-N24.json")
CERT_L0053 = Path("certificates/CERT-L0053-sos-gram-onepol-cr0008-N24.json")
GRAM_1234 = Path("tools/sos_julia/data/band_1234/gram_export.json")


@dataclass
class StructuredSosRow:
    id: str
    subclass: str
    D: int | None
    C_shor: float | None
    C_ub_sos: float | None
    C_ub_n5_hi: float | None
    evidence: str
    cert_id: str | None
    closes_vs_Cdagger: bool


@dataclass
class GalerkinBoundL0054:
    lemma_id: str = "L-0054"
    route: str = "B"
    status: str = "validated"
    evidence_level: str = "N5"
    n: int = 24
    C_dagger: float = 0.0
    C_rank1_full: float = C_RANK1_FULL
    C_fullsym_full: float = C_FULLSYM_FULL_DEALIAS
    C_greedy_41_shor: float = C_GREEDY_41_SHELLS
    best_sos_n5_hi: float = 0.0
    best_sos_subclass: str = ""
    rows: list[StructuredSosRow] | None = None
    all_ic_open: bool = True
    clay_implication: str = (
        "None. Synthesis of finite structured SOS/Shor bounds; not all-IC; not Clay."
    )
    notes: str = ""

    def as_dict(self) -> dict:
        d = asdict(self)
        if self.rows:
            d["rows"] = [asdict(r) for r in self.rows]
        return d


def _load_scale() -> dict:
    if not SCALE_JSON.is_file():
        return {}
    return json.loads(SCALE_JSON.read_text(encoding="utf-8"))


def _n5_hi_from_gram(path: Path) -> float | None:
    if not path.is_file():
        return None
    gram = load_gram_export(path)
    if not all(verify_gram_blocks(gram).values()):
        return None
    _, hi = certified_C_ub_hi(float(gram["ub_raw_abs"]))
    return hi


def _interval_hi_from_tssos(path: Path) -> float | None:
    if not path.is_file():
        return None
    d = json.loads(path.read_text(encoding="utf-8"))
    ub_raw = d.get("ub_raw_abs")
    if ub_raw is None:
        return None
    _, hi = certified_C_ub_hi(float(ub_raw))
    return hi


def lemma_l0054(n: int = 24) -> GalerkinBoundL0054:
    Cd = lemma_l0027(n=n, empirical=False).C_dagger
    scale = _load_scale()

    rows: list[StructuredSosRow] = [
        StructuredSosRow(
            id="rank1-full",
            subclass="full-dealias rank-1",
            D=6748,
            C_shor=None,
            C_ub_sos=C_RANK1_FULL,
            C_ub_n5_hi=None,
            evidence="N7 lower",
            cert_id=None,
            closes_vs_Cdagger=True,
        ),
        StructuredSosRow(
            id="band-123-n5",
            subclass="Hermitian {1,2,3}",
            D=52,
            C_shor=C_FULLSYM_BAND_123,
            C_ub_sos=C_UB_BAND_123,
            C_ub_n5_hi=_n5_hi_from_gram(Path("tools/sos_julia/data/band_123/gram_export.json")),
            evidence="N5 Gram",
            cert_id="CERT-L0052-sos-gram-band123-N24" if CERT_L0052.is_file() else None,
            closes_vs_Cdagger=C_UB_BAND_123 < Cd,
        ),
        StructuredSosRow(
            id="band-1234-n5",
            subclass="Hermitian {1,2,3,4}",
            D=64,
            C_shor=1.36986977843755,
            C_ub_sos=0.8940401324113963,
            C_ub_n5_hi=_n5_hi_from_gram(GRAM_1234),
            evidence="N5 Gram" if GRAM_1234.is_file() else "N2 TSSOS",
            cert_id="CERT-L0054-sos-gram-band1234-N24" if GRAM_1234.is_file() else None,
            closes_vs_Cdagger=True,
        ),
        StructuredSosRow(
            id="band-12345-n5",
            subclass="Hermitian {1,2,3,4,5}",
            D=112,
            C_shor=2.148097061972368,
            C_ub_sos=1.3052615184180565,
            C_ub_n5_hi=_n5_hi_from_gram(Path("tools/sos_julia/data/band_12345/gram_export.json")),
            evidence="N5 Gram",
            cert_id="CERT-L0058-sos-gram-band12345-N24"
            if Path("certificates/CERT-L0058-sos-gram-band12345-N24.json").is_file()
            else None,
            closes_vs_Cdagger=True,
        ),
        StructuredSosRow(
            id="band-123456-n5",
            subclass="Hermitian {1,2,3,4,5,6}",
            D=160,
            C_shor=2.9933259094191538,
            C_ub_sos=1.7009322214824112,
            C_ub_n5_hi=_n5_hi_from_gram(Path("tools/sos_julia/data/band_123456/gram_export.json"))
            or _interval_hi_from_tssos(Path("tools/sos_julia/data/band_123456/tssos_result.json")),
            evidence="N5 Gram"
            if Path("tools/sos_julia/data/band_123456/gram_export.json").is_file()
            else "N5 interval",
            cert_id="CERT-L0060-sos-gram-band123456-N24"
            if Path("certificates/CERT-L0060-sos-gram-band123456-N24.json").is_file()
            else None,
            closes_vs_Cdagger=True,
        ),
        StructuredSosRow(
            id="band-12345678-n5",
            subclass="Hermitian {1,2,3,4,5,6,7,8}",
            D=184,
            C_shor=3.2876874682501214,
            C_ub_sos=1.8518343630712857,
            C_ub_n5_hi=_n5_hi_from_gram(Path("tools/sos_julia/data/band_12345678/gram_export.json"))
            or _interval_hi_from_tssos(
                Path("tools/sos_julia/data/band_12345678/cosmo_tssos_result.json")
            ),
            evidence="N5 Gram"
            if Path("tools/sos_julia/data/band_12345678/gram_export.json").is_file()
            else "N5 interval COSMO",
            cert_id="CERT-L0061-sos-band12345678-N24"
            if Path("certificates/CERT-L0061-sos-band12345678-N24.json").is_file()
            else None,
            closes_vs_Cdagger=True,
        ),
        StructuredSosRow(
            id="band-123456789-n5",
            subclass="Hermitian {1..9}",
            D=244,
            C_shor=3.991247855265094,
            C_ub_sos=2.223426008475949,
            C_ub_n5_hi=_interval_hi_from_tssos(
                Path("tools/sos_julia/data/band_123456789/cosmo_tssos_result.json")
            ),
            evidence="N5 interval COSMO",
            cert_id="CERT-L0062-sos-band123456789-N24"
            if Path("certificates/CERT-L0062-sos-band123456789-N24.json").is_file()
            else None,
            closes_vs_Cdagger=True,
        ),
        StructuredSosRow(
            id="band-1to10-n5",
            subclass="Hermitian {1..10}",
            D=292,
            C_shor=4.260132170346936,
            C_ub_sos=2.4401954034218707,
            C_ub_n5_hi=_interval_hi_from_tssos(
                Path("tools/sos_julia/data/band_1to10/cosmo_tssos_result.json")
            ),
            evidence="N5 interval COSMO",
            cert_id="CERT-L0062-sos-band1to10-N24"
            if Path("certificates/CERT-L0062-sos-band1to10-N24.json").is_file()
            else None,
            closes_vs_Cdagger=True,
        ),
        StructuredSosRow(
            id="band-1to12-n5",
            subclass="Hermitian {1..12}",
            D=356,
            C_shor=4.912140950058888,
            C_ub_sos=2.782805118332653,
            C_ub_n5_hi=_interval_hi_from_tssos(
                Path("tools/sos_julia/data/band_1to12/cosmo_tssos_result.json")
            ),
            evidence="N5 interval COSMO",
            cert_id="CERT-L0062-sos-band1to12-N24"
            if Path("certificates/CERT-L0062-sos-band1to12-N24.json").is_file()
            else None,
            closes_vs_Cdagger=True,
        ),
        StructuredSosRow(
            id="onepol-cr0008-n5",
            subclass=f"one-pol C-R-0008 {list(SUPPORT_CR0008)}",
            D=92,
            C_shor=C_FULLSYM_ONEPOL,
            C_ub_sos=C_UB_ONEPOL,
            C_ub_n5_hi=_n5_hi_from_gram(Path("tools/sos_julia/data/onepol_cr0008/gram_export.json")),
            evidence="N5 Gram",
            cert_id="CERT-L0053-sos-gram-onepol-cr0008-N24" if CERT_L0053.is_file() else None,
            closes_vs_Cdagger=C_UB_ONEPOL < Cd,
        ),
        StructuredSosRow(
            id="greedy-sos-extrap",
            subclass="C-R-0012 greedy 41 shells (SOS extrap)",
            D=1772,
            C_shor=C_GREEDY_41_SHELLS,
            C_ub_sos=5.3885979,
            C_ub_n5_hi=None,
            evidence="N7 feasibility L-0055",
            cert_id=None,
            closes_vs_Cdagger=True,
        ),
        StructuredSosRow(
            id="greedy-41-shor",
            subclass="C-R-0012 greedy 41 shells",
            D=1772,
            C_shor=C_GREEDY_41_SHELLS,
            C_ub_sos=None,
            C_ub_n5_hi=None,
            evidence="N7 Shor",
            cert_id="CERT-L0047-greedy-fullsym-C0007-N24",
            closes_vs_Cdagger=C_GREEDY_41_SHELLS < Cd,
        ),
        StructuredSosRow(
            id="full-dealias-shor",
            subclass="all-IC physical fullsym",
            D=6748,
            C_shor=C_FULLSYM_FULL_DEALIAS,
            C_ub_sos=None,
            C_ub_n5_hi=None,
            evidence="N7 Shor techo",
            cert_id=None,
            closes_vs_Cdagger=False,
        ),
    ]

    extrap = scale.get("extrapolation", {})
    if extrap.get("C_ub_sos_extrapolated_full"):
        rows.append(
            StructuredSosRow(
                id="brute-sos-extrap",
                subclass="Hermitian bands → full dealias (extrapolated)",
                D=6748,
                C_shor=C_FULLSYM_FULL_DEALIAS,
                C_ub_sos=float(extrap["C_ub_sos_extrapolated_full"]),
                C_ub_n5_hi=None,
                evidence="N2 techo",
                cert_id=None,
                closes_vs_Cdagger=False,
            )
        )

    n5_rows = [r for r in rows if r.C_ub_n5_hi is not None]
    best = min(n5_rows, key=lambda r: r.C_ub_n5_hi or float("inf")) if n5_rows else rows[0]
    best_hi = best.C_ub_n5_hi if best.C_ub_n5_hi is not None else float("nan")

    notes = (
        f"Structured SOS synthesis: best N5 bound C_ub_hi≈{best_hi:.4f} "
        f"on {best.subclass}; greedy Shor≈{C_GREEDY_41_SHELLS:.2f}≤C_†; "
        f"fullsym all-IC≈{C_FULLSYM_FULL_DEALIAS:.2f}≫C_†; "
        f"brute SOS extrap≈{extrap.get('C_ub_sos_extrapolated_full', '?')} > C_†."
    )
    return GalerkinBoundL0054(
        n=n,
        C_dagger=Cd,
        best_sos_n5_hi=float(best.C_ub_n5_hi or 0.0),
        best_sos_subclass=best.subclass,
        rows=rows,
        notes=notes,
    )


def save_lemma_l0054(
    bound: GalerkinBoundL0054,
    path: str | Path = "conjectures/active/L-0054.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(bound.as_dict(), indent=2), encoding="utf-8")
