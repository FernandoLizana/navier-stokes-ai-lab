"""
Bounded stress campaign against C-S-0002 (shell IC |k|_inf<=4, full-grid evolution).

Reuses the same polish/CMA machinery as the C-S-0001 attack, but targets the
successor bound M≈48.16. Records the best *gated* enstrophy and updates the
active conjecture record honestly:
  - refuted  -> archive + successor C-S-0003,
  - survived -> keep M (do NOT loosen), log best gated support.

Scope: finite-dimensional shell-IC class. NOT continuum NS. NOT Clay.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from ns_exploration.conjectures.cs0001_shell import ConjectureCS0001, save_cs0001
from ns_exploration.conjectures.cs0001_attack import _shell_ic
from ns_exploration.conjectures.l0005_shell_galerkin import project_to_shell
from ns_exploration.conjectures.l0006_improved_shell import lemma_l0006
from ns_exploration.optimization.adjoint_ascent import adjoint_ascent
from ns_exploration.optimization.cmaes_enstrophy import evolve_enstrophy_es
from ns_exploration.optimization.enstrophy_gradient import _project_fixed_energy
from ns_exploration.optimization.polish_candidate import polish_candidate


def stress_cs0002(
    active_path: str | Path = "conjectures/active/C-S-0002.json",
    k_shell: int = 4,
    nu: float = 0.1,
    dt: float = 1e-3,
    t_end: float = 0.02,
    energy: float = 0.5,
) -> dict:
    active_path = Path(active_path)
    cur = json.loads(active_path.read_text(encoding="utf-8"))
    M_target = float(cur["proposed_bound_M"])

    runs: list[dict] = []
    best_gated = -1.0
    best_label = None

    # 1) Structured shell ICs + long polish (the previous refuter was abc_n16).
    for n in (12, 16):
        for kind in ("taylor_green", "abc", "tubes"):
            label = f"{kind}_n{n}"
            u0 = _shell_ic(n, k_shell, energy, kind)
            _, pres = polish_candidate(
                u0, M_target=M_target, nu=nu, dt=dt, t_end=t_end,
                energy_target=energy, n_steps=24, step_size=0.08,
            )
            runs.append({"label": label, "type": "structured", "polish": pres.as_dict()})
            if pres.gates_passed and pres.J_etd > best_gated:
                best_gated, best_label = float(pres.J_etd), label

    # 2) Random shell ICs + polish.
    for n, seed in ((12, 0), (12, 1), (16, 0), (16, 2)):
        label = f"random_n{n}_s{seed}"
        u0 = _shell_ic(n, k_shell, energy, "random", seed=seed)
        _, pres = polish_candidate(
            u0, M_target=M_target, nu=nu, dt=dt, t_end=t_end,
            energy_target=energy, n_steps=18, step_size=0.07,
        )
        runs.append({"label": label, "type": "random", "polish": pres.as_dict()})
        if pres.gates_passed and pres.J_etd > best_gated:
            best_gated, best_label = float(pres.J_etd), label

    # 3) CMA then double polish.
    for label, kw in (
        ("cma_n16_k4", dict(n=16, k_max=4, generations=22, population=12, seed=41)),
        ("cma_n12_k4", dict(n=12, k_max=4, generations=25, population=12, seed=42)),
    ):
        u_cma, cres = evolve_enstrophy_es(
            M_target=M_target, nu=nu, dt=dt, t_end=t_end, energy_target=energy, **kw
        )
        u0 = _project_fixed_energy(project_to_shell(u_cma, k_shell), energy)
        u_mid, _ = adjoint_ascent(
            u0, nu=nu, dt=dt, t_end=t_end, energy_target=energy, n_steps=14, step_size=0.08
        )
        _, pres = polish_candidate(
            u_mid, M_target=M_target, nu=nu, dt=dt, t_end=t_end,
            energy_target=energy, n_steps=14, step_size=0.06,
        )
        runs.append({"label": label, "type": "cma_double", "cma": cres.as_dict(),
                     "polish": pres.as_dict()})
        if pres.gates_passed and pres.J_etd > best_gated:
            best_gated, best_label = float(pres.J_etd), label

    refuted = bool(best_gated > M_target)
    successor = None
    if refuted:
        c = ConjectureCS0001(
            id="C-S-0002", proposed_bound_M=M_target, k_shell=k_shell,
            numerical_support_max=best_gated, status="refuted", refuted=True,
            refutation_value=best_gated, evidence_level="N2",
            notes=f"Refuted by bounded stress ({best_label}).",
        )
        new_M = float(max(1.5 * best_gated, best_gated + 1.0))
        c3 = ConjectureCS0001(
            id="C-S-0003", proposed_bound_M=new_M, k_shell=k_shell,
            numerical_support_max=best_gated, status="exploring", parent="C-S-0002",
            notes="Successor after C-S-0002 refutation.",
        )
        save_cs0001(c, "conjectures/rejected/C-S-0002.json")
        save_cs0001(c3, "conjectures/active/C-S-0003.json")
        active_path.unlink(missing_ok=True)
        successor = c3.as_dict()
        conj = c.as_dict()
    else:
        # Survived: keep M, refresh the recorded numerical support if higher.
        support = max(float(cur.get("numerical_support_max", 0.0)), best_gated)
        cur["numerical_support_max"] = support
        cur["notes"] = (
            f"Survived bounded stress: best gated Ω={best_gated:.6f} ({best_label}) "
            f"< M={M_target:.6f}. Still unproved (N6)."
        )
        active_path.write_text(json.dumps(cur, indent=2), encoding="utf-8")
        conj = cur

    l0006 = {
        f"N{n}_k{k}": lemma_l0006(n=n, k_shell=k, t=t_end, nu=nu).as_dict()
        for n in (12, 16)
        for k in (2, 3, 4)
    }

    return {
        "M_target": M_target,
        "best_gated": best_gated,
        "best_label": best_label,
        "refuted": refuted,
        "successor": successor,
        "conjecture": conj,
        "L0006": l0006,
        "n_runs": len(runs),
        "runs": runs,
        "route": "B",
        "not_a_clay_claim": True,
    }
