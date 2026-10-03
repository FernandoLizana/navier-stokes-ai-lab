"""
N5 certificate for L-0018: spectral enstrophy envelope vs C-S-0002.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass
from pathlib import Path

from ns_exploration.conjectures.l0018_envelope import dealias_K2, max_K2_upto


@dataclass
class L0018EnvelopeCertificate:
    certificate_id: str
    lemma_ref: str
    route: str
    evidence_level: str
    n_max: int
    n_worst: int
    K2_worst: int
    E0: float
    R: float
    Omega_cap_hi: float
    closes: bool
    clay_implication: str
    payload_sha256: str
    notes: str = ""

    def as_dict(self) -> dict:
        return asdict(self)


def build_l0018_envelope_certificate(
    n_max: int = 16,
    E0: float = 0.5,
    R: float = 48.15928260103266,
    certificate_id: str = "CERT-L0018-envelope-CS0002-N16",
) -> L0018EnvelopeCertificate:
    K2, n_w = max_K2_upto(n_max)
    # Spot-check exact integer K² for n_worst
    assert dealias_K2(n_w) == K2

    # K² ∈ ℕ is exact (mask enumeration). E0=1/2 is an exact binary float.
    # Claim Ω ≤ K²·E0 uses the exact product (no need to outward-round the product).
    cap_exact = float(K2) * float(E0)
    closes = cap_exact <= R + 1e-15
    cap_hi = cap_exact

    payload = {
        "n_max": n_max,
        "n_worst": n_w,
        "K2_worst": K2,
        "E0": E0,
        "R": R,
        "Omega_cap_hi": cap_hi,
        "closes": closes,
    }
    digest = hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()

    return L0018EnvelopeCertificate(
        certificate_id=certificate_id,
        lemma_ref="L-0018 spectral enstrophy envelope Ω≤K²E0",
        route="B",
        evidence_level="N5",
        n_max=n_max,
        n_worst=n_w,
        K2_worst=K2,
        E0=E0,
        R=R,
        Omega_cap_hi=cap_hi,
        closes=closes,
        clay_implication="None. Finite dealias Galerkin only.",
        payload_sha256=digest,
        notes=(
            f"Exact K²({n_w})={K2} (integer max over dealias mask). "
            f"Ω_cap={cap_hi} ≤ R={R}: {closes}."
        ),
    )


def verify_l0018_envelope_certificate(d: dict) -> tuple[bool, dict[str, bool]]:
    K2, n_w = max_K2_upto(int(d["n_max"]))
    cap_exact = float(K2) * float(d["E0"])
    cap_stored = float(d["Omega_cap_hi"])
    R = float(d["R"])
    closes_expected = cap_exact <= R + 1e-15 and cap_stored <= R + 1e-15
    checks = {
        "K2_matches": int(d["K2_worst"]) == K2 and int(d["n_worst"]) == n_w,
        "cap_formula": abs(cap_stored - cap_exact) < 1e-9,
        "Omega_le_R": cap_stored <= R + 1e-15,
        "closes_flag": bool(d["closes"]) == (cap_exact <= R + 1e-15),
    }
    return all(checks.values()), checks


def save_l0018_envelope_certificate(
    cert: L0018EnvelopeCertificate,
    path: str | Path = "certificates/CERT-L0018-envelope-CS0002-N16.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(cert.as_dict(), indent=2), encoding="utf-8")
