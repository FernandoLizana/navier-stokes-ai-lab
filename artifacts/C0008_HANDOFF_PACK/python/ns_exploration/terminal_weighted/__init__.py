"""Terminal-weighted stretch functional for C-0008 research (N2/N5 pipeline)."""

from ns_exploration.terminal_weighted.constants import FROZEN
from ns_exploration.terminal_weighted.weights import (
    stokes_floor_hi,
    terminal_weight,
    terminal_weight_at_T_equals_r,
)

__all__ = [
    "FROZEN",
    "terminal_weight",
    "terminal_weight_at_T_equals_r",
    "stokes_floor_hi",
]
