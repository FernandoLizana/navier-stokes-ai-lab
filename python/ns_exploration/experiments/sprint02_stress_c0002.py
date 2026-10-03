"""Sprint 2.6: L-0004 shell bounds + stress C-0002."""

from __future__ import annotations

import json
from pathlib import Path

from ns_exploration.conjectures.c0002_stress import stress_c0002
from ns_exploration.conjectures.l0004_spectral_support import shell_a_priori_bounds, save_l0004


def run(out_dir: str | Path = "experiments/exploratory/sprint02_stress_c0002") -> Path:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    # Document shell table for low-mode families (target O(10)-O(100) caps)
    shells = {}
    for n in (12, 16):
        for k in (1, 2, 3, 4):
            b = shell_a_priori_bounds(n, k)
            shells[f"N{n}_k{k}"] = b.as_dict()

    save_l0004(shell_a_priori_bounds(12, 3), "conjectures/proved_restricted/L-0004_N12_k3.json")

    Path("conjectures/proved_restricted/L-0004_spectral_support.md").write_text(
        """# L-0004 — Spectral support restricted Galerkin bounds

**Evidence:** N7 (shell hypothesis); field-conditional diagnostics N2  
**Clay:** none

## Idea

Replace global `(K,M)` from full dealias mask with `(K_S, M_S)` for fields supported in `|k|≤k_shell`.
CMA/ES IC families with `k_max≤4` fall under these hypotheses.

## Example caps (E0=0.5, t=0.02, ν=0.1) — see experiment summary for full table

Low shells can yield **O(10)–O(10³)** ceilings vs **O(10⁴)+** for full grid.

Field-conditional `(K_eff, M_eff)` can be tighter still (diagnostic only).
""",
        encoding="utf-8",
    )

    result = stress_c0002()
    payload = {"purpose": "L0004_and_stress_C0002", "shell_table": shells, **result}
    path = out / "summary.json"
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(
        json.dumps(
            {
                "refuted": result["refuted"],
                "best_gated": result["best_gated"],
                "M_target": result["M_target"],
                "best_label": result["best_label"],
                "N12_k3_best_cap": shells["N12_k3"]["best_cap"],
                "N12_k4_best_cap": shells["N12_k4"]["best_cap"],
                "successor": result.get("successor_C0003"),
            },
            indent=2,
        )
    )
    return path


if __name__ == "__main__":
    run()
