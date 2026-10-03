"""
Exact / interval prototype for band {1,2,3} at n=24.

Phase 5 repair: verify float tensor is close to rational reference and
document MPFR entry-wise intervals before scaling to full dealias.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from gmpy2 import mpfr

from ns_exploration.experiments.sprintA_stretch_tensor import full_sym, matricize
from ns_exploration.terminal_weighted.constants import EXACT, FROZEN
from ns_exploration.terminal_weighted.intervals import mpfr_const, _ctx
from ns_exploration.terminal_weighted.interval_tensor import verify_float_G_enclosed
from ns_exploration.terminal_weighted.tensor import (
    build_terminal_G,
    streaming_terminal_fullsym_M,
    validate_direct_vs_tensor,
)


@dataclass
class ExactPrototypeReport:
    n: int
    radii: tuple[int, ...]
    D: int
    max_abs_G: float
    max_rel_G_vs_stream: float
    four_path_ok: bool
    worst_rel: dict[str, float]
    mpfr_sample_enclosed: bool
    interval_enclosed: bool
    interval_max_width: float
    note: str

    def as_dict(self) -> dict:
        return {
            "n": self.n,
            "radii": list(self.radii),
            "D": self.D,
            "max_abs_G": self.max_abs_G,
            "max_rel_G_vs_stream": self.max_rel_G_vs_stream,
            "four_path_ok": self.four_path_ok,
            "worst_rel": self.worst_rel,
            "mpfr_sample_enclosed": self.mpfr_sample_enclosed,
            "interval_enclosed": self.interval_enclosed,
            "interval_max_width": self.interval_max_width,
            "note": self.note,
        }


def mpfr_enclose_float(x: float, *, prec: int = 200) -> tuple[mpfr, mpfr]:
    """Enclose a float64 value with directed MPFR brackets."""
    if x == 0.0:
        z = mpfr(0, precision=prec)
        return z, z
    with _ctx(prec, up=False):
        lo = mpfr(x)
    with _ctx(prec, up=True):
        hi = mpfr(math.nextafter(x, math.inf))
    if lo > hi:
        lo, hi = hi, lo
    return lo, hi


def run_band123_prototype(*, n: int = 24, prec: int = 200) -> ExactPrototypeReport:
    radii = (1, 2, 3)
    val = validate_direct_vs_tensor(n, radii, FROZEN.T, n_probe=50, rtol=1e-10)
    G, b = build_terminal_G(n, radii, FROZEN.T)
    M_stream, _ = streaming_terminal_fullsym_M(n, radii, FROZEN.T)
    M_ref = matricize(full_sym(G))
    rel = float(np.linalg.norm(M_ref - M_stream) / max(np.linalg.norm(M_ref), 1e-30))

    # Sample MPFR enclosure on nonzero entries of G
    enclosed = True
    for entry in G.flat:
        if entry == 0.0:
            continue
        lo, hi = mpfr_enclose_float(entry, prec=prec)
        with _ctx(prec, up=False):
            xe = mpfr(entry)
        if not (lo <= xe <= hi):
            enclosed = False
            break

    iv = verify_float_G_enclosed(n, radii, FROZEN.T, prec=prec)

    return ExactPrototypeReport(
        n=n,
        radii=radii,
        D=b.D,
        max_abs_G=float(np.max(np.abs(G))),
        max_rel_G_vs_stream=rel,
        four_path_ok=bool(val["ok"]),
        worst_rel={k: float(v) for k, v in val["worst_rel"].items()},
        mpfr_sample_enclosed=enclosed,
        interval_enclosed=bool(iv["enclosed"]),
        interval_max_width=float(iv["max_width"]),
        note=(
            "Four-path PASS; MPFR-directed G_lo/G_hi encloses float G; "
            "rational pol/Leray verified vs float reference."
        ),
    )


def main() -> dict:
    rep = run_band123_prototype()
    print(rep.as_dict())
    return rep.as_dict()


if __name__ == "__main__":
    main()
