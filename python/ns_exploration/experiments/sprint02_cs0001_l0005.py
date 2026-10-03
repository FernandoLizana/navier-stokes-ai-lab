"""Sprint 2.7: C-S-0001 shell conjecture + L-0005 truncated Galerkin lemma."""

from __future__ import annotations

import json
from pathlib import Path

from ns_exploration.conjectures.cs0001_stress import stress_cs0001


def run(out_dir: str | Path = "experiments/exploratory/sprint02_cs0001_l0005") -> Path:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    Path("conjectures/proved_restricted/L-0005_shell_galerkin.md").write_text(
        """# L-0005 — Shell-truncated Galerkin enstrophy bound

**Status:** proved_restricted  
**Evidence:** N7 (via L-0001 + L-0004)  
**Clay:** none

## Statement

Let \(P_{K_0}\) be the Fourier projector onto modes with \(|k|\\le K_0\) (and dealias mask).
For the Galerkin ODE \(\\partial_t u = P_{K_0}[-P(u\\cdot\\nabla)u + \\nu\\Delta u]\),
the support remains in the shell and

\[
\\Omega(t)\\le \\min\\bigl(\\text{L-0001}(K_0,M_{K_0},t),\\;\\text{L-0003}(K_0,M_{K_0})\\bigr).
\]

Example (E0=0.5, ν=0.1, t=0.02, N=12): \(K_0=4\) gives a finite explicit ceiling
(see `L-0005.json`).

## Relation to C-S-0001

C-S-0001 allows **full** dealias evolution from shell ICs (higher modes may appear).
L-0005 does **not** prove C-S-0001; it proves the truncated companion problem.
""",
        encoding="utf-8",
    )

    result = stress_cs0001(M_target=18.0, k_shell=4)
    path = out / "summary.json"
    path.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(
        json.dumps(
            {
                "refuted": result["refuted"],
                "best_gated": result["best_gated"],
                "M_target": result["M_target"],
                "best_label": result["best_label"],
                "L0005_N12_k4": result["L0005"]["N12_k4"]["Omega_bound"],
                "L0005_N12_k3": result["L0005"]["N12_k3"]["Omega_bound"],
                "C_S_0001_status": result["conjecture"]["status"],
            },
            indent=2,
        )
    )
    return path


if __name__ == "__main__":
    run()
