"""Sprint 2.4: L-0002 viscous Galerkin bound + CMA-lite attack on C-0001."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

from ns_exploration.conjectures.c0001_enstrophy_bound import ConjectureC0001, save_conjecture
from ns_exploration.conjectures.l0001_galerkin_bound import bound_for_resolution, save_lemma
from ns_exploration.conjectures.l0002_viscous_galerkin import (
    bound_for_resolution_l0002,
    save_lemma_l0002,
)
from ns_exploration.initial_conditions.generators import save_ic
from ns_exploration.optimization.cmaes_enstrophy import evolve_enstrophy_es
from ns_exploration.optimization.enstrophy_gradient import _project_fixed_energy
from ns_exploration.spectral.leray import divergence_l2


def run(out_dir: str | Path = "experiments/exploratory/sprint02_l0002_cmaes") -> Path:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    M_c = 10.994055804604205
    # Analytical bounds
    l1 = {str(n): bound_for_resolution(n).as_dict() for n in (8, 12, 16)}
    l2 = {str(n): bound_for_resolution_l0002(n).as_dict() for n in (8, 12, 16)}
    save_lemma(bound_for_resolution(12), "conjectures/proved_restricted/L-0001.json")
    save_lemma_l0002(bound_for_resolution_l0002(12), "conjectures/proved_restricted/L-0002.json")

    Path("conjectures/proved_restricted/L-0002_viscous_galerkin.md").write_text(
        f"""# L-0002 — Viscous Galerkin enstrophy comparison ODE

**Status:** proved_restricted (finite Fourier–Galerkin, mean-zero)  
**Evidence:** N7  
**Clay:** none

## Statement

With \(C = 2\\sqrt{{2}}\\sqrt{{3M}}\), \(a = C/2\), \(\\Omega(0)\\le K^2 E_0\), \(z_0=\\sqrt{{\\Omega(0)}}\),

\[
\\frac{{d\\Omega}}{{dt}} \\le -2\\nu\\Omega + C\\Omega^{{3/2}}
\]

via \(\\|\\nabla u\\|_\\infty \\le \\sqrt{{3M}}\\sqrt{{2\\Omega}}\) and \(\\|\\nabla\\omega\\|^2 \\ge 2\\Omega\).
The comparison ODE for \(z=\\sqrt{{\\Omega}}\) has solution

\[
z(t) = \\frac{{\\nu}}{{a - (a - \\nu/z_0)e^{{\\nu t}}}}
\]

while the denominator stays positive; else the *estimate* fails to close (not a proof of singularity).

## Caps on NS-MRL grids (E0=0.5, t=0.02, ν=0.1)

| N | L-0001 cap | L-0002 cap / status |
|---|------------|---------------------|
| 8 | {l1['8']['Omega_t_cap']:.6e} | {l2['8']['Omega_t_cap'] if l2['8']['Omega_t_cap'] is not None else 'ESTIMATE FAILS (t*=' + str(l2['8']['t_star_estimate']) + ')'} |
| 12 | {l1['12']['Omega_t_cap']:.6e} | {l2['12']['Omega_t_cap'] if l2['12']['Omega_t_cap'] is not None else 'ESTIMATE FAILS (t*=' + str(l2['12']['t_star_estimate']) + ')'} |
| 16 | {l1['16']['Omega_t_cap']:.6e} | {l2['16']['Omega_t_cap'] if l2['16']['Omega_t_cap'] is not None else 'ESTIMATE FAILS (t*=' + str(l2['16']['t_star_estimate']) + ')'} |

Auditor note: failure of L-0002 to close at large M is expected and must not be misread as Galerkin blow-up.
""",
        encoding="utf-8",
    )

    # CMA-lite attacks
    attacks = {}
    best_overall = -1.0
    best_u = None
    for label, kwargs in (
        ("n12_k3", dict(n=12, k_max=3, generations=20, population=10, seed=0)),
        ("n12_k4", dict(n=12, k_max=4, generations=15, population=8, seed=1)),
        ("n16_k3", dict(n=16, k_max=3, generations=12, population=8, seed=2)),
    ):
        u, res = evolve_enstrophy_es(M_target=M_c, nu=0.1, dt=1e-3, t_end=0.02, **kwargs)
        attacks[label] = res.as_dict()
        if res.best_J > best_overall:
            best_overall = res.best_J
            best_u = u

    refuted = best_overall > M_c
    if refuted:
        conj = ConjectureC0001(
            proposed_bound_M=M_c,
            numerical_support_max=best_overall,
            status="refuted",
            refuted=True,
            refutation_value=best_overall,
            evidence_level="N2",
            notes="Refuted by CMA-lite coefficient ES (gated only by ETD objective).",
        )
    else:
        # Optionally tighten if best is close
        new_M = M_c
        notes = (
            f"CMA-lite attack best ETD enstrophy {best_overall:.6f} < M={M_c}. "
            f"L-0002 viscous estimate often fails to close at N≥8 (auditor: not blow-up). "
            f"L-0001 remains the constructive but loose ceiling."
        )
        if best_overall > 8.491:
            # update support only
            pass
        conj = ConjectureC0001(
            proposed_bound_M=new_M,
            numerical_support_max=max(best_overall, 8.49102383722302),
            status="exploring",
            notes=notes,
        )
        conj.statement = (
            f"For all divergence-free Fourier fields on the N^3 grid with N<={conj.n_max}, "
            f"dealiased convective term, nu={conj.nu}, kinetic energy = {conj.energy}, evolved with "
            f"ETD-RK2 (dt=1e-3) to T={conj.t_end}, the enstrophy E_ens(T) <= M "
            f"with M={conj.proposed_bound_M}. FINITE-DIMENSIONAL only."
        )

    save_conjecture(conj, "conjectures/active/C-0001.json")

    if best_u is not None:
        from ns_exploration.initial_conditions.generators import ICMetadata

        payload = np.concatenate([best_u.real.ravel(), best_u.imag.ravel()])
        h = hashlib.sha256(payload.tobytes()).hexdigest()[:16]
        meta = ICMetadata(
            name="cmaes_best",
            n=int(best_u.shape[-1]),
            seed=None,
            energy=0.5,
            enstrophy=0.0,
            helicity=0.0,
            symmetries=["es"],
            resolution=int(best_u.shape[-1]),
            content_hash=h,
        )
        save_ic("datasets/candidates/cmaes_best.npz", best_u, meta)

    summary = {
        "route": "B",
        "purpose": "L0002_and_CMA_attack_C0001",
        "L0001": l1,
        "L0002": l2,
        "cma_attacks": attacks,
        "best_overall": best_overall,
        "refuted": refuted,
        "conjecture_C0001": conj.as_dict(),
        "div_best": None if best_u is None else divergence_l2(best_u),
        "not_a_clay_claim": True,
    }
    path = out / "summary.json"
    path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(
        json.dumps(
            {
                "best_overall": best_overall,
                "refuted": refuted,
                "C0001_M": conj.proposed_bound_M,
                "C0001_status": conj.status,
                "L0002_N8_fails": l2["8"]["finite_time_blowup_of_estimate"],
                "L0002_N12_fails": l2["12"]["finite_time_blowup_of_estimate"],
                "L0001_N12": l1["12"]["Omega_t_cap"],
            },
            indent=2,
        )
    )
    return path


if __name__ == "__main__":
    run()
