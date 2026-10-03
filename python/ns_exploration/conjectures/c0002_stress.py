"""
Stress-test C-0002 (M≈20.32) with CMA + adjoint polish at N=12 and N=16.
"""

from __future__ import annotations

import json
from pathlib import Path

from ns_exploration.conjectures.c0002_enstrophy_bound import ConjectureC0002, save_conjecture_c0002
from ns_exploration.conjectures.l0004_spectral_support import (
    field_conditional_bounds,
    shell_a_priori_bounds,
)
from ns_exploration.initial_conditions.generators import ICMetadata, save_ic
from ns_exploration.optimization.cmaes_enstrophy import evolve_enstrophy_es
from ns_exploration.optimization.polish_candidate import polish_candidate
from ns_exploration.spectral.fourier_conventions import kinetic_energy_from_hat


def stress_c0002(
    M_target: float = 20.320076249551253,
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
        ("n12_k3", dict(n=12, k_max=3, generations=22, population=10, seed=10)),
        ("n12_k4", dict(n=12, k_max=4, generations=18, population=10, seed=11)),
        ("n16_k3", dict(n=16, k_max=3, generations=15, population=8, seed=12)),
        ("n16_k4", dict(n=16, k_max=4, generations=12, population=8, seed=13)),
    ]

    for label, kw in configs:
        u_cma, cres = evolve_enstrophy_es(
            M_target=M_target,
            nu=nu,
            dt=dt,
            t_end=t_end,
            energy_target=energy,
            **kw,
        )
        u_pol, pres = polish_candidate(
            u_cma,
            M_target=M_target,
            nu=nu,
            dt=dt,
            t_end=t_end,
            energy_target=energy,
            n_steps=10,
            step_size=0.05,
        )
        fcond = field_conditional_bounds(u_pol, E0=energy, t=t_end, nu=nu)
        runs.append(
            {
                "label": label,
                "cma": cres.as_dict(),
                "polish": pres.as_dict(),
                "field_L0004_best_cap": fcond.best_cap,
                "field_L0004": fcond.as_dict(),
            }
        )
        if pres.gates_passed and pres.J_etd > best_gated:
            best_gated = pres.J_etd
            best_u = u_pol
            best_label = label

    shell_bounds = {
        f"N{n}_k{k}": shell_a_priori_bounds(n, k).as_dict()
        for n in (12, 16)
        for k in (2, 3, 4)
    }

    refuted = bool(best_gated > M_target)
    if refuted:
        c2 = ConjectureC0002(
            proposed_bound_M=M_target,
            numerical_support_max=best_gated,
            status="refuted",
            refuted=True,
            refutation_value=best_gated,
            evidence_level="N2",
            notes=f"Refuted by gated polish ({best_label}).",
        )
        new_M = float(max(1.5 * best_gated, best_gated + 1.0))
        c3 = ConjectureC0002(
            id="C-0003",
            proposed_bound_M=new_M,
            numerical_support_max=best_gated,
            status="exploring",
            parent_refuted="C-0002",
            notes=f"Successor after C-0002 refutation. M=1.5× gated max {best_gated:.6f}.",
        )
        c3.statement = (
            f"For all divergence-free Fourier fields on the N^3 grid with N<={c3.n_max}, "
            f"dealiased convective term, nu={c3.nu}, kinetic energy = {c3.energy}, evolved with "
            f"ETD-RK2 (dt=1e-3) to T={c3.t_end}, the enstrophy E_ens(T) <= M "
            f"with M={c3.proposed_bound_M}. FINITE-DIMENSIONAL only. Replaces refuted C-0002."
        )
        save_conjecture_c0002(c2, "conjectures/rejected/C-0002.json")
        save_conjecture_c0002(c3, "conjectures/active/C-0003.json")
        # Keep C-0002 in active as refuted copy for audit
        save_conjecture_c0002(c2, "conjectures/active/C-0002.json")
        successor = c3.as_dict()
    else:
        margin = M_target - best_gated
        c2 = ConjectureC0002(
            proposed_bound_M=M_target,
            numerical_support_max=best_gated,
            status="exploring",
            notes=(
                f"Survived stress at N=12/16. Best gated ETD={best_gated:.6f}, "
                f"margin={margin:.6f}. Shell bounds logged in L-0004."
            ),
        )
        save_conjecture_c0002(c2, "conjectures/active/C-0002.json")
        successor = None

    if best_u is not None:
        import hashlib

        payload = __import__("numpy").concatenate([best_u.real.ravel(), best_u.imag.ravel()])
        h = hashlib.sha256(payload.tobytes()).hexdigest()[:16]
        meta = ICMetadata(
            name=f"stress_best_{best_label}",
            n=int(best_u.shape[-1]),
            seed=None,
            energy=kinetic_energy_from_hat(best_u),
            enstrophy=0.0,
            helicity=0.0,
            symmetries=["stress_c0002"],
            resolution=int(best_u.shape[-1]),
            content_hash=h,
        )
        save_ic(f"datasets/candidates/stress_c0002_{best_label}.npz", best_u, meta)

    return {
        "M_target": M_target,
        "best_gated": best_gated,
        "best_label": best_label,
        "refuted": refuted,
        "successor_C0003": successor,
        "runs": runs,
        "shell_bounds_L0004": shell_bounds,
        "route": "B",
        "evidence_level": "N2",
    }
