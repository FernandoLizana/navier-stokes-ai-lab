"""Evaluate C-0008 repair closure vs sharp targets."""

from __future__ import annotations

import json
from pathlib import Path

from gmpy2 import mpfr

from ns_exploration.terminal_weighted.constants import FROZEN
from ns_exploration.terminal_weighted.intervals import i_star_interval
from ns_exploration.validation.c0008_verify import verify_terminal_certificate

REPO = Path(__file__).resolve().parents[2]
PY = REPO / "python"


def _load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def evaluate_repair_closure(*, n: int = 24, prec: int = 128) -> dict:
    cert_path = REPO / f"certificates/CERT-C0008-repair-full-n{n}-one-inf.json"
    manifest_path = PY / f"experiments/terminal_weighted/shell_manifest_repair_full_one_inf_n{n}.json"
    cert = _load_json(cert_path)
    manifest = _load_json(manifest_path)
    ok, issues = verify_terminal_certificate(
        cert, manifest=manifest, manifest_path=manifest_path, prec=prec
    )

    I_lo, I_hi = i_star_interval(prec=prec)
    I_decl = mpfr(cert["integral_interval"][1], precision=prec)
    omega_decl = mpfr(cert["omega_T_interval"][1], precision=prec)
    M_target = mpfr(FROZEN.c0008_target_M, precision=prec)

    gap_I = float(I_decl / I_lo)
    gap_omega = float(omega_decl / M_target)

    historical = {
        "L-0073_historical_I_hi": 37.91955812524674097702040111940917993968,
        "L-0073_historical_valid": False,
        "note": "Historical cluster-best NOT rigorous (invalid B_C*max(phi), fro/min bound).",
    }

    return {
        "verify_ok": ok,
        "verify_issues": issues,
        "scope": cert.get("scope"),
        "n_shells": cert["mode_set"]["n_shells"],
        "D": cert["mode_set"]["D"],
        "bound_method": cert.get("bound_method"),
        "I_hi_repair": str(I_decl),
        "I_star_lo": str(I_lo),
        "I_star_hi_cert": cert["i_star_interval"][1],
        "omega_T_hi_repair": str(omega_decl),
        "M_target": str(M_target),
        "closes_strict": bool(I_decl < I_lo),
        "closes_omega": bool(omega_decl < M_target),
        "closes_c0008": cert.get("closes_c0008"),
        "gap_I_over_I_star": gap_I,
        "gap_omega_over_M": gap_omega,
        "historical_comparison": historical,
        "conjecture_C0008": "exploring",
        "verdict": (
            "NO-GO: repair full-dealias one_inf bound is rigorous (verify PASS) "
            f"but I_hi={I_decl} >> I_*≈{I_lo} (~{gap_I:.1f}x); C-0008 not closed."
        ),
    }


def main() -> dict:
    rep = evaluate_repair_closure()
    out = REPO / "reports" / "C0008_REPAIR_EVALUATION.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(rep, indent=2), encoding="utf-8")
    print(json.dumps(rep, indent=2))
    return rep


if __name__ == "__main__":
    main()
