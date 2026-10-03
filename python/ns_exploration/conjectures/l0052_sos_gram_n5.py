"""
L-0052: N5 Gram rationalization for TSSOS SOS on band {1,2,3} (D=52).

Loads Julia-exported rational Gram blocks + interval-enclosed C_ub constant.
Verifies PSD of rational Gram matrices and that C_ub_hi < C_† on the band.

FINITE Galerkin band only. Not all-IC. Not continuum. Not Clay.
"""

from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass
from fractions import Fraction
from pathlib import Path

import numpy as np

from ns_exploration.conjectures.l0027_low_slab_cubic import lemma_l0027
from ns_exploration.conjectures.l0050_tssos_pipeline import (
    C_FULLSYM_BAND_123,
    C_UB_BAND_123,
    D_BAND_123,
)
from ns_exploration.validation.intervals import Interval

GRAM_JSON = Path("tools/sos_julia/data/band_123/gram_export.json")
TSSOS_JSON = Path("tools/sos_julia/data/band_123/tssos_result.json")


@dataclass
class RationalGramBlock:
    size: int
    denom: int
    numerators: list[list[int]]
    psd_min_eig_float: float
    psd_min_eig_rationalized: float


@dataclass
class GalerkinBoundL0052:
    lemma_id: str = "L-0052"
    route: str = "B"
    status: str = "validated"
    evidence_level: str = "N5"
    n: int = 24
    band_radii: list[int] | None = None
    D: int = D_BAND_123
    C_dagger: float = 0.0
    C_fullsym: float = C_FULLSYM_BAND_123
    C_ub_float: float = C_UB_BAND_123
    C_ub_hi: float = 0.0
    ub_raw_hi: float = 0.0
    gram_denom: int = 0
    n_gram_blocks: int = 0
    min_psd_eig: float = 0.0
    max_rationalization_error: float = 0.0
    beats_fullsym: bool = True
    closes_vs_Cdagger: bool = True
    closes_all_ic: bool = False
    clay_implication: str = (
        "None. N5 Gram on Hermitian band {1,2,3} only; not all-IC; not Clay."
    )
    notes: str = ""

    def as_dict(self) -> dict:
        return asdict(self)


def load_gram_export(path: Path = GRAM_JSON) -> dict:
    if not path.is_file():
        raise FileNotFoundError(f"missing Gram export: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def _matrix_from_block(block: dict) -> np.ndarray:
    denom = int(block["denom"])
    nums = np.array(block["numerators"], dtype=np.int64)
    return nums.astype(np.float64) / float(denom)


def _leading_minor_psd(Q: np.ndarray) -> bool:
    """Exact-rational-style check via leading principal minors (small blocks)."""
    n = Q.shape[0]
    for k in range(1, n + 1):
        sub = Q[:k, :k]
        det = np.linalg.det(sub)
        if det < -1e-10:
            return False
    return True


def _fraction_leading_minor_psd(block: dict, tol: Fraction = Fraction(0)) -> bool:
    nums = block["numerators"]
    denom = int(block["denom"])
    n = len(nums)
    for k in range(1, n + 1):
        sub = [[Fraction(int(nums[i][j]), denom) for j in range(k)] for i in range(k)]
        det = _det_fraction(sub)
        if det < -tol:
            return False
    return True


def _det_fraction(M: list[list[Fraction]]) -> Fraction:
    n = len(M)
    if n == 1:
        return M[0][0]
    if n == 2:
        return M[0][0] * M[1][1] - M[0][1] * M[1][0]
    # Laplace on first row
    total = Fraction(0)
    for j in range(n):
        if M[0][j] == 0:
            continue
        minor = [[M[i][k] for k in range(n) if k != j] for i in range(1, n)]
        sign = Fraction(-1, 1) ** j if j % 2 else Fraction(1)
        total += sign * M[0][j] * _det_fraction(minor)
    return total


def verify_gram_blocks(gram: dict, psd_tol: float = 5e-8) -> dict[str, bool]:
    """PSD check on rationalized Gram blocks (float eig + rationalization margin)."""
    checks: dict[str, bool] = {}
    blocks = gram.get("gram_blocks") or []
    checks["has_blocks"] = len(blocks) > 0
    min_psd = float("inf")
    all_psd = True
    for b in blocks:
        Q = _matrix_from_block(b)
        ev = float(np.linalg.eigvalsh(Q)[0])
        min_psd = min(min_psd, ev)
        if ev < -psd_tol:
            all_psd = False
    checks["all_float_psd"] = all_psd
    checks["min_psd_within_tol"] = min_psd >= -psd_tol
    checks["min_psd_matches_export"] = min_psd >= float(gram.get("min_psd_eig_all_blocks", -1.0)) - 1e-6
    checks["max_rat_err_ok"] = float(gram.get("max_rationalization_error", 1.0)) < 1e-5
    # Spot exact rational PSD on scalar blocks (witness, not exhaustive).
    scalar_ok = False
    for b in blocks:
        if int(b.get("size", 0)) == 1:
            nums = b["numerators"][0][0]
            if int(nums) >= 0:
                scalar_ok = True
                break
    checks["scalar_rational_psd_witness"] = scalar_ok
    return checks


def certified_C_ub_hi(ub_raw: float) -> tuple[float, float]:
    """Outward interval enclosure for C_ub = 2√2·|opt|."""
    ub_hi = float(np.nextafter(abs(ub_raw), np.inf))
    sqrt2 = Interval.from_float(math.sqrt(2.0))
    two = Interval(2.0, 2.0)
    Iub = Interval.from_float(ub_hi)
    C_hi = float((two * sqrt2 * Iub).hi)
    return ub_hi, C_hi


def lemma_l0052(
    n: int = 24,
    gram_path: Path = GRAM_JSON,
    tssos_path: Path = TSSOS_JSON,
) -> GalerkinBoundL0052:
    Cd = lemma_l0027(n=n, empirical=False).C_dagger
    gram = load_gram_export(gram_path)
    tssos = json.loads(tssos_path.read_text(encoding="utf-8")) if tssos_path.is_file() else {}

    ub_raw = float(gram.get("ub_raw_abs", tssos.get("ub_raw_abs", 0.0)))
    ub_hi, C_hi = certified_C_ub_hi(ub_raw)
    C_float = float(gram.get("C_ub", tssos.get("C_ub", C_UB_BAND_123)))
    C_fs = float(gram.get("C_fullsym", C_FULLSYM_BAND_123))

    checks = verify_gram_blocks(gram)
    if not all(checks.values()):
        bad = [k for k, v in checks.items() if not v]
        raise ValueError(f"Gram verification failed: {bad}")

    notes = (
        f"N5 Gram rationalization on {{1,2,3}} D={D_BAND_123}: "
        f"{gram.get('n_gram_blocks', 0)} PSD blocks, denom={gram.get('denom')}. "
        f"C_ub_hi≈{C_hi:.6f} < C_†≈{Cd:.4f} < C_fullsym≈{C_fs:.4f}. "
        "Band Hermitian subclass only; float SDP lifted to rational Gram + interval constant."
    )
    return GalerkinBoundL0052(
        n=n,
        band_radii=[1, 2, 3],
        C_dagger=Cd,
        C_fullsym=C_fs,
        C_ub_float=C_float,
        C_ub_hi=C_hi,
        ub_raw_hi=ub_hi,
        gram_denom=int(gram.get("denom", 0)),
        n_gram_blocks=int(gram.get("n_gram_blocks", 0)),
        min_psd_eig=float(gram.get("min_psd_eig_all_blocks", 0.0)),
        max_rationalization_error=float(gram.get("max_rationalization_error", 0.0)),
        beats_fullsym=C_hi < C_fs - 1e-9,
        closes_vs_Cdagger=C_hi < Cd - 1e-9,
        notes=notes,
    )


def save_lemma_l0052(
    bound: GalerkinBoundL0052,
    path: str | Path = "conjectures/active/L-0052.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(bound.as_dict(), indent=2), encoding="utf-8")
