"""Frozen constants for C-0008 — exact decimal strings for MPFR certificates."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class ExactDecimalConstants:
    """Authoritative string forms for interval/MPFR arithmetic."""

    NU: str = "0.1"
    E0: str = "0.5"
    T: str = "0.02"
    M_TARGET: str = "41.283999509509286"
    M_C0007_PROVED: str = "41.74337593971423"
    STOKES_FLOOR_FROZEN: str = "40.82462307930434"
    DELTA_M: str = "0.4593764302049479"
    C_DAGGER: str = "9.562009938690146"
    L0048_FULLSYM_T: str = "25.925922025377965"
    INTEGRAL_THRESHOLD: str = "1.2993127556607498"


EXACT = ExactDecimalConstants()


@dataclass(frozen=True)
class FrozenConstants:
    """Legacy float accessors — audit code only; certificates use EXACT strings."""

    nu: float = 0.1
    E0: float = 0.5
    T: float = 0.02
    n_max: int = 24
    K2: int = 147
    r_next: int = 134
    gap: int = 13
    stokes_floor: float = 40.82462307930434
    c0007_proved_M: float = 41.74337593971423
    c0008_target_M: float = 41.283999509509286
    delta_M: float = 0.4593764302049479
    C_dagger: float = 9.562009938690146
    l0048_fullsym_T: float = 25.925922025377965
    integral_threshold: float = 1.2993127556607498
    uniform_threshold: float = 64.9656377830375
    D_full: int = 6748


FROZEN = FrozenConstants()

# Frobenius re-enabled after entry aggregation fix (P0 #1 repair)
DEFAULT_BOUND_METHOD: str = "one_inf_only"
FROBENIUS_STREAMING_ENABLED: bool = True
