"""Sprint: L-0055 greedy SOS feasibility scan."""

from __future__ import annotations

import json

from ns_exploration.conjectures.l0055_greedy_sos_feasibility import lemma_l0055, save_lemma_l0055


def main() -> None:
    b = lemma_l0055()
    save_lemma_l0055(b)
    out = {
        "lemma_id": b.lemma_id,
        "greedy_D": b.greedy_D,
        "C_greedy_shor": b.C_greedy_shor,
        "C_ub_sos_extrap": b.C_ub_sos_extrap,
        "extrap_closes_cr0012": b.extrap_closes_cr0012,
        "est_sdp_variables_M": round(b.est_sdp_variables / 1e6, 2),
        "est_solve_hours": round(b.est_solve_hours, 1),
        "verdict": b.verdict,
        "laptop_feasible": b.laptop_feasible,
    }
    print(json.dumps(out, indent=2), flush=True)


if __name__ == "__main__":
    main()
