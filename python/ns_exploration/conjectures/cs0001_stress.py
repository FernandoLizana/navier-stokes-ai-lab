"""
Stress C-S-0001: maximize enstrophy for shell-supported ICs under FULL evolution.
Also verify L-0005 on truncated Galerkin dynamics.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from ns_exploration.conjectures.cs0001_shell import ConjectureCS0001, save_cs0001
from ns_exploration.conjectures.l0005_shell_galerkin import (
    lemma_l0005,
    numerical_check_below_l0005,
    project_to_shell,
    save_l0005,
    evolve_galerkin_shell,
)
from ns_exploration.diagnostics.metrics import compute_diagnostics
from ns_exploration.initial_conditions.generators import ICMetadata, save_ic
from ns_exploration.optimization.adjoint_ascent import enstrophy_integrator
from ns_exploration.optimization.cmaes_enstrophy import evolve_enstrophy_es
from ns_exploration.optimization.enstrophy_gradient import _project_fixed_energy
from ns_exploration.optimization.polish_candidate import polish_candidate
from ns_exploration.spectral.fourier_conventions import kinetic_energy_from_hat
from ns_exploration.spectral.leray import divergence_l2


def stress_cs0001(
    M_target: float = 18.0,
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

    configs = [
        ("n12_k4", dict(n=12, k_max=4, generations=20, population=10, seed=20)),
        ("n12_k3", dict(n=12, k_max=3, generations=18, population=10, seed=21)),
        ("n16_k4", dict(n=16, k_max=4, generations=14, population=8, seed=22)),
    ]

    for label, kw in configs:
        # CMA already restricts DOFs to k_max; project to Euclidean shell too
        u_cma, cres = evolve_enstrophy_es(
            M_target=M_target,
            nu=nu,
            dt=dt,
            t_end=t_end,
            energy_target=energy,
            **kw,
        )
        u_cma = _project_fixed_energy(project_to_shell(u_cma, k_shell), energy)
        # Polish on full grid (may leave shell) — this is the C-S-0001 attack
        u_pol, pres = polish_candidate(
            u_cma,
            M_target=M_target,
            nu=nu,
            dt=dt,
            t_end=t_end,
            energy_target=energy,
            n_steps=8,
            step_size=0.05,
        )
        # Truncated Galerkin comparison for L-0005
        u_trunc_T = evolve_galerkin_shell(u_cma, nu, dt, t_end, k_shell)
        omega_trunc = compute_diagnostics(u_trunc_T, nu).enstrophy
        check = numerical_check_below_l0005(u_cma, k_shell, nu=nu, dt=dt, t_end=t_end)

        runs.append(
            {
                "label": label,
                "cma": cres.as_dict(),
                "polish_full": pres.as_dict(),
                "Omega_truncated_galerkin": omega_trunc,
                "L0005_check": check,
            }
        )
        if pres.gates_passed and pres.J_etd > best_gated:
            best_gated = float(pres.J_etd)
            best_u = u_pol
            best_label = label

    # L-0005 statements for k=3 and k=4
    lemmas = {
        f"N{n}_k{k}": lemma_l0005(n=n, k_shell=k, t=t_end, nu=nu).as_dict()
        for n in (12, 16)
        for k in (3, 4)
    }
    save_l0005(lemma_l0005(12, 4), "conjectures/proved_restricted/L-0005.json")

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
            notes=f"Refuted by gated full-grid polish from shell IC ({best_label}).",
        )
        new_M = float(max(1.5 * best_gated, best_gated + 1.0))
        c_next = ConjectureCS0001(
            id="C-S-0002",
            proposed_bound_M=new_M,
            k_shell=k_shell,
            numerical_support_max=best_gated,
            status="exploring",
            parent="C-S-0001",
            notes=f"Successor after C-S-0001 refutation. M=1.5× gated max.",
        )
        save_cs0001(c, "conjectures/rejected/C-S-0001.json")
        save_cs0001(c, "conjectures/active/C-S-0001.json")
        save_cs0001(c_next, "conjectures/active/C-S-0002.json")
        successor = c_next.as_dict()
    else:
        # Tighten explicit M when margin is large (still unproved)
        new_M = float(M_target)
        if best_gated > 0 and M_target > 1.5 * best_gated:
            new_M = float(1.5 * best_gated)
        c = ConjectureCS0001(
            proposed_bound_M=new_M,
            k_shell=k_shell,
            numerical_support_max=best_gated,
            status="exploring",
            notes=(
                f"Survived shell-IC stress. Best gated full-grid Ω={best_gated:.6f}. "
                f"Tightened M: {M_target} → {new_M} (=1.5× best if margin large). "
                f"L-0005 proves truncated Galerkin bounds "
                f"(N12 k3 ≈ {lemmas['N12_k3']['Omega_bound']:.4g}, "
                f"k4 ≈ {lemmas['N12_k4']['Omega_bound']:.4g}). "
                f"Full-grid C-S-0001 still open (shell not invariant)."
            ),
        )
        save_cs0001(c, "conjectures/active/C-S-0001.json")
        successor = None

    if best_u is not None:
        import hashlib

        payload = np.concatenate([best_u.real.ravel(), best_u.imag.ravel()])
        h = hashlib.sha256(payload.tobytes()).hexdigest()[:16]
        meta = ICMetadata(
            name=f"cs0001_best_{best_label}",
            n=int(best_u.shape[-1]),
            seed=None,
            energy=kinetic_energy_from_hat(best_u),
            enstrophy=0.0,
            helicity=0.0,
            symmetries=["shell_ic"],
            resolution=int(best_u.shape[-1]),
            content_hash=h,
        )
        save_ic(f"datasets/candidates/cs0001_{best_label}.npz", best_u, meta)

    return {
        "M_target": M_target,
        "k_shell": k_shell,
        "best_gated": best_gated,
        "best_label": best_label,
        "refuted": refuted,
        "conjecture": c.as_dict(),
        "successor": successor,
        "L0005": lemmas,
        "runs": runs,
        "route": "B",
        "evidence_level": "N2",
        "not_a_clay_claim": True,
    }
