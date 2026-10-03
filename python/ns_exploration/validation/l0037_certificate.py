"""Certificate for L-0037 shell-diff / bi-shell Shor techo."""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path

from ns_exploration.conjectures.l0037_shell_diff_stretch import (
    bishell_shor_C_lower,
    lemma_l0037,
    mono_radial_stretch_identity,
)


@dataclass
class CertificateL0037:
    cert_id: str
    lemma_id: str
    n: int
    c0007_M: float
    C_dagger: float
    mono_stretch_zero: bool
    mono_identity: str
    C_near_pair: float
    near_pair: list[int] | None
    C_far_pair: float
    far_pair: list[int] | None
    C_bishell_shor_max_sample: float
    closes_c0007: bool
    clay_implication: str

    def as_dict(self) -> dict:
        return asdict(self)


def build_l0037_certificate(n: int = 24) -> CertificateL0037:
    # Compact scan for cert build; witnesses (1,2) and (2,134) explicit
    b = lemma_l0037(n=n, iters=14, max_pairs=40)
    return CertificateL0037(
        cert_id="CERT-L0037-shell-diff-C0007-N24",
        lemma_id="L-0037",
        n=b.n,
        c0007_M=b.c0007_M,
        C_dagger=b.C_dagger,
        mono_stretch_zero=b.mono_stretch_zero,
        mono_identity=mono_radial_stretch_identity(),
        C_near_pair=b.C_near_pair,
        near_pair=list(b.near_pair) if b.near_pair else None,
        C_far_pair=b.C_far_pair,
        far_pair=list(b.far_pair) if b.far_pair else None,
        C_bishell_shor_max_sample=b.C_bishell_shor_max_sample,
        closes_c0007=b.closes_c0007,
        clay_implication=b.clay_implication,
    )


def verify_l0037_certificate(d: dict) -> tuple[bool, dict]:
    near = d.get("near_pair") or [1, 2]
    far = d.get("far_pair") or [2, 134]
    c_near = float(bishell_shor_C_lower(int(d["n"]), int(near[0]), int(near[1]), iters=12)["C_shor_lower"])
    c_far = float(bishell_shor_C_lower(int(d["n"]), int(far[0]), int(far[1]), iters=12)["C_shor_lower"])
    checks = {
        "mono_zero": d["mono_stretch_zero"] is True,
        "near_le_Cdagger": c_near <= float(d["C_dagger"]) * 1.05,  # slack for iters
        "far_above_Cdagger": c_far > float(d["C_dagger"]),
        "sample_max_above_Cdagger": float(d["C_bishell_shor_max_sample"]) > float(d["C_dagger"]),
        "closes_false": d["closes_c0007"] is False,
    }
    return all(checks.values()), checks


def save_l0037_certificate(
    cert: CertificateL0037,
    path: str | Path = "certificates/CERT-L0037-shell-diff-C0007-N24.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(cert.as_dict(), indent=2), encoding="utf-8")
