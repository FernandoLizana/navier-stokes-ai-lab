"""Sprint 2.3: adjoint ascent, C-0001 adversarial pass, tighten M."""

from __future__ import annotations

import json
from pathlib import Path

from ns_exploration.conjectures.c0001_adversarial import (
    adversarial_refutation_campaign,
    update_c0001_after_campaign,
)
from ns_exploration.conjectures.c0001_enstrophy_bound import save_conjecture
from ns_exploration.initial_conditions.generators import random_div_free, save_ic
from ns_exploration.optimization.adjoint_ascent import adjoint_ascent
from ns_exploration.optimization.enstrophy_gradient import _project_fixed_energy


def run(out_dir: str | Path = "experiments/exploratory/sprint02_adjoint_ascent") -> Path:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    # Demo ascent on one seed with gates
    u0, meta = random_div_free(12, seed=1, energy_target=0.5)
    u_fin, asc = adjoint_ascent(
        u0, nu=0.1, dt=1e-3, t_end=0.02, energy_target=0.5, n_steps=6, step_size=0.08
    )
    save_ic(str(Path("datasets/candidates") / "adjoint_ascent_seed1_n12.npz"), u_fin, meta)

    # Adversarial campaign against current loose M=20, then tighten
    campaign = adversarial_refutation_campaign(
        M=20.0, n=12, n_seeds=8, ascent_steps=5, nu=0.1, dt=1e-3, t_end=0.02, energy=0.5
    )
    # Drop bulky per-step histories in printed summary but keep in file
    conj = update_c0001_after_campaign(campaign, tighten_factor=1.5)
    save_conjecture(conj, "conjectures/active/C-0001.json")
    # Archive previous loose bound note
    Path("conjectures/active/C-0001_history.md").write_text(
        f"# C-0001 history\n\nTightened from M=20 after adversarial max ETD enstrophy "
        f"{campaign['best_J_etd']:.6f}. New M={conj.proposed_bound_M}. "
        f"Status={conj.status}. Compactness lemma recorded in notes.\n",
        encoding="utf-8",
    )

    payload = {
        "route": "B",
        "purpose": "adjoint_ascent_and_C0001_adversarial",
        "demo_ascent": asc.as_dict(),
        "campaign_summary": {
            k: campaign[k]
            for k in (
                "M",
                "n",
                "best_J_etd",
                "best_seed",
                "refuted_with_gates",
                "refutation_value",
                "n_records",
                "n_gates_passed",
                "evidence_level",
                "route",
            )
        },
        "campaign_records": campaign["records"],
        "conjecture_C0001": conj.as_dict(),
        "not_a_clay_claim": True,
    }
    path = out / "summary.json"
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(
        json.dumps(
            {
                "demo_J0": asc.J_history[0],
                "demo_J_final": asc.J_final,
                "demo_gates": asc.gates_passed,
                "best_J_etd": campaign["best_J_etd"],
                "refuted": campaign["refuted_with_gates"],
                "new_M": conj.proposed_bound_M,
                "status": conj.status,
            },
            indent=2,
        )
    )
    return path


if __name__ == "__main__":
    run()
