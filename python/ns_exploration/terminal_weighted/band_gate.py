"""Extend band gate to {1..6} and rational basis checks."""

from __future__ import annotations

from ns_exploration.terminal_weighted.constants import FROZEN
from ns_exploration.terminal_weighted.exact_prototype import run_band123_prototype
from ns_exploration.terminal_weighted.interval_tensor import verify_float_G_enclosed, verify_rational_vs_float_interval
from ns_exploration.terminal_weighted.rational_basis import verify_pol_leray_against_float
from ns_exploration.terminal_weighted.tensor import validate_direct_vs_tensor


def run_band6_gate(*, n: int = 24, prec: int = 200, n_probe: int = 20) -> dict:
    radii = (1, 2, 3, 4, 5, 6)
    val = validate_direct_vs_tensor(n, radii, FROZEN.T, n_probe=n_probe, rtol=1e-9)
    iv = verify_float_G_enclosed(n, radii, FROZEN.T, prec=prec)
    rat = verify_rational_vs_float_interval(n, radii, FROZEN.T, prec=prec)
    pol = verify_pol_leray_against_float(prec=prec)
    return {
        "n": n,
        "radii": list(radii),
        "four_path_ok": bool(val["ok"]),
        "worst_rel": {k: float(v) for k, v in val["worst_rel"].items()},
        "interval_enclosed": bool(iv["enclosed"]),
        "rational_interval_ok": bool(rat["ok"]),
        "rational_pol_ok": bool(pol["ok"]),
        "D": val["D"],
    }


def main() -> dict:
    b123 = run_band123_prototype()
    b6 = run_band6_gate()
    print({"band123": b123.as_dict(), "band6": b6})
    return {"band123": b123.as_dict(), "band6": b6}


if __name__ == "__main__":
    main()
