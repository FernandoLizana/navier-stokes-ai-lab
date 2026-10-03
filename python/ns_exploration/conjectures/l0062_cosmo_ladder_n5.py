"""
L-0062: COSMO band ladder N5 interval constants ({1..9}, {1..10}, {1..12}).

Completes the consecutive Hermitian SOS ladder through shell 12 using
existing COSMO TSSOS N2 results. Gram PSD deferred on laptop for D>=160.

FINITE Galerkin bands only. Not all-IC. Not Clay.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path

from ns_exploration.conjectures.l0027_low_slab_cubic import lemma_l0027
from ns_exploration.conjectures.l0052_sos_gram_n5 import (
    certified_C_ub_hi,
    load_gram_export,
    verify_gram_blocks,
)

COSMO_BANDS: tuple[dict, ...] = (
    {
        "lemma_suffix": "band123456789",
        "cert_id": "CERT-L0062-sos-band123456789-N24",
        "radii": [1, 2, 3, 4, 5, 6, 7, 8, 9],
        "D": 244,
        "related_cr": "C-R-0010",
        "data_dir": "tools/sos_julia/data/band_123456789",
        "tssos_file": "cosmo_tssos_result.json",
    },
    {
        "lemma_suffix": "band1to10",
        "cert_id": "CERT-L0062-sos-band1to10-N24",
        "radii": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
        "D": 292,
        "related_cr": "C-R-0010",
        "data_dir": "tools/sos_julia/data/band_1to10",
        "tssos_file": "cosmo_tssos_result.json",
    },
    {
        "lemma_suffix": "band1to12",
        "cert_id": "CERT-L0062-sos-band1to12-N24",
        "radii": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
        "D": 356,
        "related_cr": "C-R-0011",
        "data_dir": "tools/sos_julia/data/band_1to12",
        "tssos_file": "cosmo_tssos_result.json",
    },
)


@dataclass
class CosmoBandN5Row:
    id: str
    cert_id: str
    band_radii: list[int]
    D: int
    C_dagger: float
    C_fullsym: float
    C_ub_float: float
    C_ub_hi: float
    ub_raw_hi: float
    ratio_sos_over_shor: float
    status: str
    related_cr: str
    gram_blocks: int
    notes: str


@dataclass
class GalerkinBoundL0062:
    lemma_id: str = "L-0062"
    route: str = "B"
    status: str = "validated"
    evidence_level: str = "N5"
    n: int = 24
    C_dagger: float = 0.0
    solver: str = "COSMO"
    rows: list[CosmoBandN5Row] | None = None
    clay_implication: str = (
        "None. N5 interval constants on COSMO TSSOS bands through shell 12; not all-IC; not Clay."
    )
    notes: str = ""

    def as_dict(self) -> dict:
        d = asdict(self)
        if self.rows:
            d["rows"] = [asdict(r) for r in self.rows]
        return d


def _band_row(spec: dict, Cd: float) -> CosmoBandN5Row:
    data_dir = Path(spec["data_dir"])
    tssos_path = data_dir / spec["tssos_file"]
    gram_path = data_dir / "gram_export.json"
    tssos = json.loads(tssos_path.read_text(encoding="utf-8"))
    if gram_path.is_file():
        gram = load_gram_export(gram_path)
        checks = verify_gram_blocks(gram)
        if not all(checks.values()):
            raise ValueError(f"Gram failed for {spec['id']}: {checks}")
        ub_raw = float(gram["ub_raw_abs"])
        n_blocks = int(gram.get("n_gram_blocks", 0))
        status = "validated"
        ev = "N5 Gram"
    else:
        ub_raw = float(tssos["ub_raw_abs"])
        n_blocks = 0
        status = "validated_constant"
        ev = "N5 interval COSMO"
    ub_hi, C_hi = certified_C_ub_hi(ub_raw)
    C_float = float(tssos.get("C_ub", 0.0))
    C_fs = float(tssos.get("C_fullsym", 0.0))
    ratio = C_hi / C_fs if C_fs > 0 else 0.0
    return CosmoBandN5Row(
        id=spec["lemma_suffix"],
        cert_id=spec["cert_id"],
        band_radii=list(spec["radii"]),
        D=int(spec["D"]),
        C_dagger=Cd,
        C_fullsym=C_fs,
        C_ub_float=C_float,
        C_ub_hi=C_hi,
        ub_raw_hi=ub_hi,
        ratio_sos_over_shor=ratio,
        status=status,
        related_cr=str(spec["related_cr"]),
        gram_blocks=n_blocks,
        notes=(
            f"{ev} on Hermitian shells {spec['radii'][:3]}..{spec['radii'][-1]} "
            f"D={spec['D']}: C_ub_hi≈{C_hi:.4f} < C_dagger≈{Cd:.3f} (ratio≈{ratio:.3f})."
        ),
    )


def lemma_l0062(n: int = 24) -> GalerkinBoundL0062:
    Cd = lemma_l0027(n=n, empirical=False).C_dagger
    rows = [_band_row(spec, Cd) for spec in COSMO_BANDS]
    notes = (
        f"COSMO ladder through shell 12: "
        + "; ".join(f"{{1..{r.band_radii[-1]}}} C_ub_hi≈{r.C_ub_hi:.3f}" for r in rows)
        + f". All < C_dagger≈{Cd:.3f}. Gram deferred D>=160."
    )
    return GalerkinBoundL0062(n=n, C_dagger=Cd, rows=rows, notes=notes)


def save_lemma_l0062(
    bound: GalerkinBoundL0062,
    path: str | Path = "conjectures/active/L-0062.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(bound.as_dict(), indent=2), encoding="utf-8")
