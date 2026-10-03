"""Frozen constants for C-0008 terminal-weighted research (spec §3)."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class FrozenConstants:
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
