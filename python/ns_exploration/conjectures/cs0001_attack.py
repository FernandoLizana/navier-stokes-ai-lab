"""
Aggressive multi-start attack on C-S-0001 (shell IC, full-grid evolution).
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from ns_exploration.conjectures.cs0001_shell import ConjectureCS0001, save_cs0001
from ns_exploration.conjectures.l0005_shell_galerkin import project_to_shell
from ns_exploration.conjectures.l0006_improved_shell import lemma_l0006, save_l0006
from ns_exploration.initial_conditions.generators import (
    ICMetadata,
    abc_flow,
    random_div_free,
    taylor_green,
    vortex_tubes_periodic,
    save_ic,
)
from ns_exploration.optimization.adjoint_ascent import adjoint_ascent, enstrophy_integrator
from ns_exploration.optimization.cmaes_enstrophy import evolve_enstrophy_es
from ns_exploration.optimization.enstrophy_gradient import _project_fixed_energy
from ns_exploration.optimization.polish_candidate import polish_candidate
from ns_exploration.spectral.fourier_conventions import kinetic_energy_from_hat
from ns_exploration.spectral.leray import divergence_l2


def _shell_ic(n: int, k_shell: int, energy: float, kind: str, seed: int = 0) -> np.ndarray:
    if kind == "random":
        uh, _ = random_div_free(n, seed=seed, k_peak=max(1, k_shell // 2), energy_target=energy)
    elif kind == "taylor_green":
        uh, _ = taylor_green(n)
    elif kind == "abc":
        uh, _ = abc_flow(n)
    elif kind == "tubes":
        uh, _ = vortex_tubes_periodic(n)
    else:
        raise ValueError(kind)
    return _project_fixed_energy(project_to_shell(uh, k_shell), energy)


def attack_cs0001(
    M_target: float = 8.894923676246325,
    k_shell: int = 4,
    nu: float = 0.1,
    dt: float = 1e-3,
    t_end: float = 0.02,
    energy: float = 0.5,
) -> dict:
    runs = []
    best_gated = -1.0
    best_u = None
    best_label = None

    # 1) Structured shell ICs + long polish
    for n in (12, 16):
        for kind in ("taylor_green", "abc", "tubes"):
            label = f"{kind}_n{n}"
            u0 = _shell_ic(n, k_shell, energy, kind)
            u1, pres = polish_candidate(
                u0,
                M_target=M_target,
                nu=nu,
                dt=dt,
                t_end=t_end,
                energy_target=energy,
                n_steps=20,
                step_size=0.08,
            )
            runs.append({"label": label, "type": "structured_polish", "polish": pres.as_dict()})
            if pres.gates_passed and pres.J_etd > best_gated:
                best_gated, best_u, best_label = float(pres.J_etd), u1, label

    # 2) Random shell ICs + long polish
    for n, seed in ((12, 0), (12, 1), (12, 2), (16, 0), (16, 3)):
        label = f"random_n{n}_s{seed}"
        u0 = _shell_ic(n, k_shell, energy, "random", seed=seed)
        u1, pres = polish_candidate(
            u0,
            M_target=M_target,
            nu=nu,
            dt=dt,
            t_end=t_end,
            energy_target=energy,
            n_steps=16,
            step_size=0.07,
        )
        runs.append({"label": label, "type": "random_polish", "polish": pres.as_dict()})
        if pres.gates_passed and pres.J_etd > best_gated:
            best_gated, best_u, best_label = float(pres.J_etd), u1, label

    # 3) CMA then long polish (shell-projected IC)
    for label, kw in (
        ("cma_n12_k4", dict(n=12, k_max=4, generations=25, population=12, seed=30)),
        ("cma_n16_k4", dict(n=16, k_max=4, generations=18, population=10, seed=31)),
        ("cma_n12_k3", dict(n=12, k_max=3, generations=20, population=12, seed=32)),
    ):
        u_cma, cres = evolve_enstrophy_es(
            M_target=M_target,
            nu=nu,
            dt=dt,
            t_end=t_end,
            energy_target=energy,
            **kw,
        )
        u0 = _project_fixed_energy(project_to_shell(u_cma, k_shell), energy)
        # Extra adjoint ascent before gated polish
        u_mid, asc = adjoint_ascent(
            u0,
            nu=nu,
            dt=dt,
            t_end=t_end,
            energy_target=energy,
            n_steps=12,
            step_size=0.08,
        )
        u1, pres = polish_candidate(
            u_mid,
            M_target=M_target,
            nu=nu,
            dt=dt,
            t_end=t_end,
            energy_target=energy,
            n_steps=12,
            step_size=0.06,
        )
        runs.append(
            {
                "label": label,
                "type": "cma_double_polish",
                "cma": cres.as_dict(),
                "ascent_J_final": asc.J_final,
                "polish": pres.as_dict(),
            }
        )
        if pres.gates_passed and pres.J_etd > best_gated:
            best_gated, best_u, best_label = float(pres.J_etd), u1, label

    # L-0006 table
    l0006 = {
        f"N{n}_k{k}": lemma_l0006(n=n, k_shell=k, t=t_end, nu=nu).as_dict()
        for n in (12, 16)
        for k in (2, 3, 4)
    }
    save_l0006(lemma_l0006(12, 3), "conjectures/proved_restricted/L-0006.json")

    refuted = bool(best_gated > M_target)
    if refuted:
        c = ConjectureCS0001(
            proposed_bound_M=M_target,
            k_shell=k_shell,
            numerical_support_max=best_gated,
            status="refuted",
            refuted=True,
            refutation_value=best_gated,
            evidence_level="N2",
            notes=f"Refuted by aggressive attack ({best_label}).",
        )
        new_M = float(max(1.5 * best_gated, best_gated + 1.0))
        c2 = ConjectureCS0001(
            id="C-S-0002",
            proposed_bound_M=new_M,
            k_shell=k_shell,
            numerical_support_max=best_gated,
            status="exploring",
            parent="C-S-0001",
            notes="Successor after C-S-0001 refutation.",
        )
        save_cs0001(c, "conjectures/rejected/C-S-0001.json")
        save_cs0001(c, "conjectures/active/C-S-0001.json")
        save_cs0001(c2, "conjectures/active/C-S-0002.json")
        successor = c2.as_dict()
    else:
        new_M = M_target
        if best_gated > 0 and M_target > 1.5 * best_gated:
            new_M = float(1.5 * best_gated)
        # If we got closer, keep max(M_target, 1.5*best) when best rose but still under M
        if best_gated > 5.92994911749755 and best_gated < M_target:
            new_M = M_target  # do not loosen; keep current M
        c = ConjectureCS0001(
            proposed_bound_M=new_M if new_M <= M_target else M_target,
            k_shell=k_shell,
            numerical_support_max=best_gated,
            status="exploring",
            notes=(
                f"Aggressive attack best gated Ω={best_gated:.6f} ({best_label}). "
                f"M={min(new_M, M_target):.6f}. L-0006 best caps logged. Still unproved."
            ),
        )
        # Prefer not loosening below previous M unless we explicitly tighten
        if best_gated < M_target / 1.5:
            c.proposed_bound_M = float(1.5 * best_gated)
        else:
            c.proposed_bound_M = M_target
        save_cs0001(c, "conjectures/active/C-S-0001.json")
        successor = None

    if best_u is not None:
        import hashlib

        payload = np.concatenate([best_u.real.ravel(), best_u.imag.ravel()])
        h = hashlib.sha256(payload.tobytes()).hexdigest()[:16]
        meta = ICMetadata(
            name=f"cs0001_attack_{best_label}",
            n=int(best_u.shape[-1]),
            seed=None,
            energy=kinetic_energy_from_hat(best_u),
            enstrophy=0.0,
            helicity=0.0,
            symmetries=["shell_attack"],
            resolution=int(best_u.shape[-1]),
            content_hash=h,
        )
        save_ic(f"datasets/candidates/cs0001_attack_{best_label}.npz", best_u, meta)

    return {
        "M_target": M_target,
        "best_gated": best_gated,
        "best_label": best_label,
        "refuted": refuted,
        "conjecture": c.as_dict(),
        "successor": successor,
        "L0006": l0006,
        "n_runs": len(runs),
        "runs": runs,
        "route": "B",
        "not_a_clay_claim": True,
    }
