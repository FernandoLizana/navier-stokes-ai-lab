"""Sprint 2 continued: adjoint check, gate battery, candidate, conjecture C-0001."""

from __future__ import annotations

import json
from pathlib import Path

from ns_exploration.conjectures.c0001_enstrophy_bound import evaluate_conjecture, save_conjecture
from ns_exploration.initial_conditions.generators import random_div_free
from ns_exploration.optimization.adjoint import verify_adjoint_vs_fd
from ns_exploration.optimization.candidate_pipeline import extract_one_candidate, run_gate_battery
from ns_exploration.spectral.leray import leray_project_hat


def run(out_dir: str | Path = "experiments/exploratory/sprint02_continued") -> Path:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    # Adjoint vs FD on a small field
    u, _ = random_div_free(12, seed=0, energy_target=0.5)
    d, _ = random_div_free(12, seed=17, energy_target=1.0)
    adj = verify_adjoint_vs_fd(u, d, nu=0.1, dt=1e-3, t_end=0.01, tol=0.15)

    gates = run_gate_battery(seeds=range(20), n_low=12, gate_tol=0.35)
    n_pass = sum(1 for g in gates if g["gate_passed"])
    candidate = extract_one_candidate(
        gates,
        out_dir="datasets/candidates",
    )

    conj = evaluate_conjecture(M=20.0, n=12)
    save_conjecture(conj, "conjectures/active/C-0001.json")

    payload = {
        "route": "B",
        "purpose": "adjoint_gate_battery_candidate_conjecture",
        "adjoint_vs_fd": adj.as_dict(),
        "gate_battery": {
            "n_seeds": len(gates),
            "n_passed": n_pass,
            "pass_rate": n_pass / max(len(gates), 1),
            "rows": gates,
        },
        "candidate": None if candidate is None else candidate.as_dict(),
        "conjecture_C0001": conj.as_dict(),
        "not_a_clay_claim": True,
    }
    path = out / "summary.json"
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps({k: payload[k] for k in payload if k != "gate_battery"}, indent=2))
    print(
        json.dumps(
            {
                "gate_n_passed": n_pass,
                "gate_pass_rate": n_pass / 20,
                "adjoint_agreed": adj.agreed,
                "adjoint_rel_err": adj.relative_error,
                "candidate_accepted": None if candidate is None else candidate.accepted,
                "conjecture_status": conj.status,
            },
            indent=2,
        )
    )
    return path


if __name__ == "__main__":
    run()
