"""MPFR rational vs float band smoke — production path readiness."""

from __future__ import annotations

import json
from pathlib import Path

from gmpy2 import mpfr

PY = Path(__file__).resolve().parents[1]
REPO = PY.parent


def main() -> dict:
    from ns_exploration.terminal_weighted.exact_prototype import run_band123_prototype
    from ns_exploration.terminal_weighted.interval_tensor import verify_rational_vs_float_interval
    from ns_exploration.terminal_weighted.opnorm_certified import certified_shell_blocks
    from ns_exploration.terminal_weighted.rational_basis import verify_pol_leray_against_float
    from ns_exploration.terminal_weighted.constants import FROZEN

    radii = (1, 2, 3)
    proto = run_band123_prototype()
    pol = verify_pol_leray_against_float(prec=200)
    iv = verify_rational_vs_float_interval(24, radii, FROZEN.T, prec=200)
    prod = certified_shell_blocks(24, radii, prec=128, workers=1, rational_coefficients=True)
    rep = {
        "radii": list(radii),
        "four_path_ok": bool(proto.four_path_ok),
        "interval_enclosed": bool(proto.interval_enclosed),
        "pol_leray_ok": bool(pol["ok"]),
        "rational_interval_ok": bool(iv["ok"]),
        "production_opnorm_rational_wired": True,
        "production_rational_C_term_sum": str(
            sum(float(b["C_term_hi"]) for b in prod["blocks"])
        ),
        "note": "opnorm_certified rational_coefficients=True on band {1,2,3}",
        "ready_for_band_regen": bool(
            proto.four_path_ok and proto.interval_enclosed and pol["ok"] and iv["ok"]
        ),
    }
    out = REPO / "reports" / "MPFR_BAND_SMOKE.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(rep, indent=2), encoding="utf-8")
    print(json.dumps(rep, indent=2))
    return rep


if __name__ == "__main__":
    main()
