"""Sprint 2 partial: residuals, energy cross-check, enstrophy FD gate."""

from __future__ import annotations

import json
from pathlib import Path

from ns_exploration.diagnostics.energy_crosscheck import (
    embedded_2d_operator_agreement,
    energy_decay_crosscheck,
    parseval_crosscheck,
    embedded_2d_taylor_green,
)
from ns_exploration.diagnostics.residuals import compare_residuals_two_resolutions
from ns_exploration.initial_conditions.generators import random_div_free
from ns_exploration.optimization.enstrophy_gradient import (
    one_gradient_ascent_step,
    resolution_doubling_gate,
)


def run(out_dir: str | Path = "experiments/exploratory/sprint02_partial") -> Path:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    residual_reports = [r.as_dict() for r in compare_residuals_two_resolutions((16, 32))]
    uh = embedded_2d_taylor_green(24)
    parseval = parseval_crosscheck(uh).as_dict()
    op2d = embedded_2d_operator_agreement(24)
    energy = energy_decay_crosscheck(16, nu=0.05, dt=1e-3, steps=40)

    gate = resolution_doubling_gate(
        n_low=12,
        seed=0,
        nu=0.1,
        dt=1e-3,
        t_end=0.02,
        energy_target=0.5,
        gate_tol=0.35,
    )

    # One ascent step only if gate passed; otherwise still record attempt as blocked
    ascent = None
    if gate.gate_passed:
        u0, _ = random_div_free(12, seed=0, energy_target=0.5)
        _u1, j0, j1 = one_gradient_ascent_step(
            u0, nu=0.1, dt=1e-3, t_end=0.02, energy_target=0.5, seed=1
        )
        ascent = {
            "performed": True,
            "J0": j0,
            "J1": j1,
            "delta_J": j1 - j0,
            "n": 12,
            "evidence_level": "N2",
            "notes": "Single FD ascent step after gate pass; not a maximizer.",
        }
    else:
        ascent = {
            "performed": False,
            "reason": "resolution_doubling_gate_failed",
            "relative_agreement": gate.relative_agreement,
            "evidence_level": "N2",
        }

    payload = {
        "route": "B",
        "purpose": "sprint02_partial_validation_before_adversarial",
        "residuals": residual_reports,
        "parseval": parseval,
        "embedded_2d_operator": op2d,
        "energy_decay": energy,
        "enstrophy_gradient_gate": gate.as_dict(),
        "ascent": ascent,
        "not_a_clay_claim": True,
    }
    path = out / "summary.json"
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps(payload, indent=2))
    return path


if __name__ == "__main__":
    run()
