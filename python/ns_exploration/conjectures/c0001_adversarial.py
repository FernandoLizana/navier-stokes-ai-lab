"""
Adversarial campaign against C-0001 + optional tightening of M.

Also records the trivial finite-N existence remark (compactness): for each fixed N
the continuous objective on the compact energy sphere attains a max — so some
finite M_N exists. C-0001 is about an *explicit* uniform-looking M, not existence.
"""

from __future__ import annotations

import numpy as np

from ns_exploration.conjectures.c0001_enstrophy_bound import ConjectureC0001
from ns_exploration.initial_conditions.generators import random_div_free, taylor_green
from ns_exploration.optimization.adjoint_ascent import adjoint_ascent, enstrophy_integrator
from ns_exploration.optimization.enstrophy_gradient import _project_fixed_energy


def adversarial_refutation_campaign(
    M: float,
    n: int = 12,
    n_seeds: int = 12,
    ascent_steps: int = 6,
    nu: float = 0.1,
    dt: float = 1e-3,
    t_end: float = 0.02,
    energy: float = 0.5,
) -> dict:
    """
    Try to exceed M via adjoint ascent from random + structured ICs.
    Also evaluates objective with ETD (the conjecture's integrator).
    """
    records = []
    best = -np.inf
    best_seed = None
    refuted = False
    refutation_value = None

    # structured start
    uh, _ = taylor_green(n)
    uh = _project_fixed_energy(uh, energy)
    u_fin, asc = adjoint_ascent(
        uh, nu=nu, dt=dt, t_end=t_end, energy_target=energy, n_steps=ascent_steps, step_size=0.08
    )
    j_etd = enstrophy_integrator(u_fin, nu, dt, t_end, "etd_rk2")
    records.append({"seed": "taylor_green", **asc.as_dict(), "J_etd_conjecture_metric": j_etd})
    best = max(best, j_etd, asc.J_final)
    if j_etd > M:
        refuted = True
        refutation_value = j_etd

    for s in range(n_seeds):
        uh, _ = random_div_free(n, seed=s, energy_target=energy)
        u_fin, asc = adjoint_ascent(
            uh,
            nu=nu,
            dt=dt,
            t_end=t_end,
            energy_target=energy,
            n_steps=ascent_steps,
            step_size=0.08,
        )
        j_etd = enstrophy_integrator(u_fin, nu, dt, t_end, "etd_rk2")
        records.append({"seed": s, **asc.as_dict(), "J_etd_conjecture_metric": j_etd})
        if j_etd > best:
            best = j_etd
            best_seed = s
        if j_etd > M and asc.gates_passed:
            refuted = True
            refutation_value = j_etd
            break
        if j_etd > M and not asc.gates_passed:
            # Soft refutation flagged but not trusted without gates
            records[-1]["soft_exceed_without_gates"] = True

    return {
        "M": M,
        "n": n,
        "best_J_etd": float(best),
        "best_seed": best_seed,
        "refuted_with_gates": refuted,
        "refutation_value": refutation_value,
        "n_records": len(records),
        "n_gates_passed": sum(1 for r in records if r.get("gates_passed")),
        "records": records,
        "evidence_level": "N2",
        "route": "B",
    }


def update_c0001_after_campaign(campaign: dict, tighten_factor: float = 1.5) -> ConjectureC0001:
    """
    If refuted → status refuted.
    Else tighten M to max(best*tighten_factor, small floor) and keep exploring.
    Also attach compactness remark.
    """
    best = float(campaign["best_J_etd"])
    compactness = (
        "LEMMA (finite-N, trivial): for each fixed N the set of div-free fields with "
        "kinetic energy = E is compact in the finite-dimensional Fourier coefficient "
        "space, and u0 |-> enstrophy(T) is continuous for the discrete dynamical system; "
        "hence a finite maximum M_N(T,nu,E) exists. This does NOT yield an explicit M, "
        "NOR uniformity in N, NOR a continuum Clay statement."
    )

    if campaign["refuted_with_gates"]:
        c = ConjectureC0001(
            proposed_bound_M=float(campaign["M"]),
            numerical_support_max=best,
            status="refuted",
            refuted=True,
            refutation_value=float(campaign["refutation_value"]),
            evidence_level="N2",
            notes=f"Refuted by adjoint ascent with gates. {compactness}",
        )
        return c

    new_M = float(max(best * tighten_factor, best + 1e-6))
    c = ConjectureC0001(
        proposed_bound_M=new_M,
        numerical_support_max=best,
        status="exploring",
        refuted=False,
        evidence_level="N6",
        notes=(
            f"No gated refutation of previous M={campaign['M']}. "
            f"Tightened explicit bound to M={new_M:.6f} (= {tighten_factor}× best adversarial "
            f"ETD enstrophy {best:.6f}). Still unproved. {compactness}"
        ),
    )
    # Refresh statement with new M embedded via field
    c.statement = (
        f"For all divergence-free Fourier fields on the N^3 grid with N<={c.n_max}, "
        f"dealiased convective term, nu={c.nu}, kinetic energy = {c.energy}, evolved with "
        f"ETD-RK2 (dt=1e-3) to T={c.t_end}, the enstrophy E_ens(T) <= M "
        f"with M={c.proposed_bound_M}. FINITE-DIMENSIONAL only. "
        f"Existence of some finite max per N is trivial by compactness; explicit M is the claim."
    )
    return c
