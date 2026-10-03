"""Sprint 2.8: aggressive C-S-0001 attack + L-0006 improved shell bound."""

from __future__ import annotations

import json
from pathlib import Path

from ns_exploration.conjectures.cs0001_attack import attack_cs0001


def run(out_dir: str | Path = "experiments/exploratory/sprint02_cs0001_attack") -> Path:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    Path("conjectures/proved_restricted/L-0006_improved_shell.md").write_text(
        """# L-0006 — Improved short-time shell bound

**Evidence:** N7 (restricted)  
**Clay:** none

## Idea

Replace `||∇u||_∞ ≤ K√(3M)√(2E)` by `||∇u||_∞ ≤ √(3M)√(2Ω)` and integrate
`dΩ/dt ≤ 2√2√(3M) Ω^{3/2}` to get an algebraic short-time ceiling when it closes;
otherwise fall back to L-0001 / L-0003 and take the best finite bound.

Still shell Galerkin only — not C-S-0001 (full-grid evolution).
""",
        encoding="utf-8",
    )

    result = attack_cs0001()
    path = out / "summary.json"
    path.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(
        json.dumps(
            {
                "refuted": result["refuted"],
                "best_gated": result["best_gated"],
                "best_label": result["best_label"],
                "M_target": result["M_target"],
                "C_S_0001_M": result["conjecture"]["proposed_bound_M"],
                "C_S_0001_status": result["conjecture"]["status"],
                "L0006_N12_k3": result["L0006"]["N12_k3"]["best_cap"],
                "L0006_N12_k3_which": result["L0006"]["N12_k3"]["which_best"],
                "L0006_N12_k2": result["L0006"]["N12_k2"]["best_cap"],
            },
            indent=2,
        )
    )
    return path


if __name__ == "__main__":
    run()
