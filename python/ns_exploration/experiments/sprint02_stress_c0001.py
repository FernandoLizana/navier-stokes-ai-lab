"""
Stress-test C-0001 at N=12 and N=16 with longer adjoint ascent.
Compare numerical M to L-0001 analytical Galerkin bound.
"""

from __future__ import annotations

import json
from pathlib import Path

from ns_exploration.conjectures.c0001_adversarial import (
    adversarial_refutation_campaign,
    update_c0001_after_campaign,
)
from ns_exploration.conjectures.c0001_enstrophy_bound import save_conjecture
from ns_exploration.conjectures.l0001_galerkin_bound import bound_for_resolution, save_lemma


def run(out_dir: str | Path = "experiments/exploratory/sprint02_stress_c0001") -> Path:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    # Current explicit numerical target
    current_M = 10.994055804604205

    campaigns = {}
    for n, seeds, steps in ((12, 10, 8), (16, 6, 6)):
        campaigns[str(n)] = adversarial_refutation_campaign(
            M=current_M,
            n=n,
            n_seeds=seeds,
            ascent_steps=steps,
            nu=0.1,
            dt=1e-3,
            t_end=0.02,
            energy=0.5,
        )

    best_overall = max(campaigns[k]["best_J_etd"] for k in campaigns)
    any_refute = any(campaigns[k]["refuted_with_gates"] for k in campaigns)

    # Merge into a synthetic campaign for updater (use best across N)
    merged = {
        "M": current_M,
        "n": 16,
        "best_J_etd": best_overall,
        "best_seed": None,
        "refuted_with_gates": any_refute,
        "refutation_value": best_overall if any_refute else None,
        "n_records": sum(campaigns[k]["n_records"] for k in campaigns),
        "n_gates_passed": sum(campaigns[k]["n_gates_passed"] for k in campaigns),
        "evidence_level": "N2",
        "route": "B",
    }

    if any_refute:
        conj = update_c0001_after_campaign(merged, tighten_factor=1.5)
    else:
        # Keep M if still above best*1.05; else tighten to 1.5*best
        if best_overall * 1.5 < current_M:
            conj = update_c0001_after_campaign(merged, tighten_factor=1.5)
        else:
            # Survive stress: keep M, update support max + notes
            from ns_exploration.conjectures.c0001_enstrophy_bound import ConjectureC0001

            conj = ConjectureC0001(
                proposed_bound_M=current_M,
                numerical_support_max=best_overall,
                status="exploring",
                notes=(
                    f"Survived longer adjoint stress at N=12 and N=16; best ETD enstrophy "
                    f"{best_overall:.6f} < M={current_M}. Still unproved. See L-0001 for a "
                    f"proved-but-loose Galerkin exponential bound."
                ),
            )
            conj.statement = (
                f"For all divergence-free Fourier fields on the N^3 grid with N<={conj.n_max}, "
                f"dealiased convective term, nu={conj.nu}, kinetic energy = {conj.energy}, evolved with "
                f"ETD-RK2 (dt=1e-3) to T={conj.t_end}, the enstrophy E_ens(T) <= M "
                f"with M={conj.proposed_bound_M}. FINITE-DIMENSIONAL only."
            )

    save_conjecture(conj, "conjectures/active/C-0001.json")

    # Analytical lemma for N=12 and N=16
    lemmas = {
        "N12": bound_for_resolution(12).as_dict(),
        "N16": bound_for_resolution(16).as_dict(),
    }
    save_lemma(bound_for_resolution(16), "conjectures/proved_restricted/L-0001.json")

    # Markdown dossier
    Path("conjectures/proved_restricted/L-0001_galerkin_enstrophy_exp.md").write_text(
        f"""# L-0001 — Galerkin enstrophy exponential bound

**Status:** proved_restricted (finite Fourier–Galerkin only)  
**Evidence:** N7  
**Route tag:** B (setting), **not** a Clay contribution

## Statement

Under the hypotheses in `python/ns_exploration/conjectures/l0001_galerkin_bound.py`,

$$\\Omega(t) \\le K^2 E_0 \\exp\\bigl(2 K \\sqrt{{3M}} \\sqrt{{2 E_0}}\\, t\\bigr).$$

## Constants on NS-MRL grids (dealiased)

| N | K | M | Ω(0.02) cap |
|---|-----|-----|-------------|
| 12 | {lemmas['N12']['K']:.6g} | {lemmas['N12']['M']} | {lemmas['N12']['Omega_t_cap']:.6e} |
| 16 | {lemmas['N16']['K']:.6g} | {lemmas['N16']['M']} | {lemmas['N16']['Omega_t_cap']:.6e} |

Compare to numerical conjecture C-0001 with M≈{conj.proposed_bound_M:.4f} (explicit tight target, unproved).
L-0001 is rigorous but far looser; it shows an explicit finite-N ceiling exists constructively.

## Clay

None. K,M grow with resolution.
""",
        encoding="utf-8",
    )

    payload = {
        "route": "B",
        "purpose": "stress_C0001_and_L0001",
        "current_M_tested": current_M,
        "campaigns": {
            k: {kk: campaigns[k][kk] for kk in campaigns[k] if kk != "records"}
            for k in campaigns
        },
        "campaign_records": {k: campaigns[k]["records"] for k in campaigns},
        "best_overall": best_overall,
        "refuted": any_refute,
        "conjecture_C0001": conj.as_dict(),
        "lemma_L0001": lemmas,
        "not_a_clay_claim": True,
    }
    path = out / "summary.json"
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(
        json.dumps(
            {
                "best_overall": best_overall,
                "refuted": any_refute,
                "C0001_M": conj.proposed_bound_M,
                "C0001_status": conj.status,
                "L0001_N16_cap": lemmas["N16"]["Omega_t_cap"],
                "L0001_N12_cap": lemmas["N12"]["Omega_t_cap"],
            },
            indent=2,
        )
    )
    return path


if __name__ == "__main__":
    run()
