"""N5 certificate for band {1,2,3,4} Gram (L-0054 extension)."""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass
from pathlib import Path

from ns_exploration.conjectures.l0052_sos_gram_n5 import (
    certified_C_ub_hi,
    load_gram_export,
    verify_gram_blocks,
)

GRAM_JSON = Path("tools/sos_julia/data/band_1234/gram_export.json")
C_FULLSYM = 1.36986977843755
C_UB = 0.8940401324113963
D_BAND = 64


@dataclass
class L0054Band1234Certificate:
    certificate_id: str
    lemma_ref: str
    route: str
    evidence_level: str
    n: int
    band_radii: list[int]
    D: int
    C_dagger: float
    C_fullsym: float
    C_ub_float: float
    C_ub_hi: float
    gram_denom: int
    n_gram_blocks: int
    min_psd_eig: float
    clay_implication: str
    payload_sha256: str

    def as_dict(self) -> dict:
        return asdict(self)


def build_l0054_band1234_certificate(n: int = 24) -> L0054Band1234Certificate:
    gram = load_gram_export(GRAM_JSON)
    if not all(verify_gram_blocks(gram).values()):
        raise ValueError("Gram verification failed for band_1234")
    _, C_hi = certified_C_ub_hi(float(gram["ub_raw_abs"]))
    payload = {"D": D_BAND, "C_ub_hi": C_hi, "radii": [1, 2, 3, 4]}
    digest = hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    return L0054Band1234Certificate(
        certificate_id="CERT-L0054-sos-gram-band1234-N24",
        lemma_ref="L-0054 extension: TSSOS + rational Gram on Hermitian {1,2,3,4}",
        route="B",
        evidence_level="N5",
        n=n,
        band_radii=[1, 2, 3, 4],
        D=D_BAND,
        C_dagger=float(gram["C_dagger"]),
        C_fullsym=C_FULLSYM,
        C_ub_float=float(gram.get("C_ub", C_UB)),
        C_ub_hi=C_hi,
        gram_denom=int(gram["denom"]),
        n_gram_blocks=int(gram["n_gram_blocks"]),
        min_psd_eig=float(gram["min_psd_eig_all_blocks"]),
        clay_implication="None. Band {1..4} only; not all-IC; not Clay.",
        payload_sha256=digest,
    )


def save_l0054_band1234_certificate(
    cert: L0054Band1234Certificate,
    path: str | Path = "certificates/CERT-L0054-sos-gram-band1234-N24.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(cert.as_dict(), indent=2), encoding="utf-8")
