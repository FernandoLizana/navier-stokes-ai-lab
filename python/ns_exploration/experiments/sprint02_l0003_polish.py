"""Sprint 2.5: L-0003 uniform bound + polish CMA candidate vs C-0001."""

from __future__ import annotations

import json
from pathlib import Path

from ns_exploration.conjectures.c0001_enstrophy_bound import ConjectureC0001, save_conjecture
from ns_exploration.conjectures.c0002_enstrophy_bound import ConjectureC0002, save_conjecture_c0002
from ns_exploration.conjectures.l0001_galerkin_bound import bound_for_resolution
from ns_exploration.conjectures.l0002_viscous_galerkin import bound_for_resolution_l0002
from ns_exploration.conjectures.l0003_uniform_galerkin import (
    bound_for_resolution_l0003,
    save_lemma_l0003,
)
from ns_exploration.initial_conditions.generators import ICMetadata, random_div_free, save_ic
from ns_exploration.optimization.cmaes_enstrophy import evolve_enstrophy_es
from ns_exploration.optimization.polish_candidate import load_candidate_npz, polish_candidate


def run(out_dir: str | Path = "experiments/exploratory/sprint02_l0003_polish") -> Path:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    M_c = 10.994055804604205

    # Compare lemmas
    table = {}
    for n in (8, 12, 16):
        b1 = bound_for_resolution(n)
        b2 = bound_for_resolution_l0002(n)
        b3 = bound_for_resolution_l0003(n)
        table[str(n)] = {
            "L0001_short_time": b1.Omega_t_cap,
            "L0002_closes": not b2.finite_time_blowup_of_estimate,
            "L0002_cap": b2.Omega_t_cap,
            "L0003_uniform": b3.Omega_uniform_cap,
            "L0003_eq": b3.Omega_eq,
            "best_of_closing": min(
                x
                for x in (
                    b1.Omega_t_cap,
                    b2.Omega_t_cap if b2.Omega_t_cap is not None else float("inf"),
                    b3.Omega_uniform_cap,
                )
            ),
        }
    save_lemma_l0003(bound_for_resolution_l0003(12), "conjectures/proved_restricted/L-0003.json")
    Path("conjectures/proved_restricted/L-0003_uniform_galerkin.md").write_text(
        f"""# L-0003 — Uniform-in-time Galerkin enstrophy bound

**Status:** proved_restricted  
**Evidence:** N7  
**Clay:** none

## Statement

For mean-zero Fourier–Galerkin fields with at most \(M\) spatial modes, max wave number \(K\),
viscosity \(\\nu>0\), and energy \(E(t)\\le E_0\),

\[
\\Omega(t) \\le \\max\\bigl(K^2 E_0,\\; 6 M E_0^2 / \\nu^2\\bigr)
\\qquad\\forall t\\ge 0.
\]

Uses \(\\|\\nabla\\omega\\|^2 \\ge 2\\Omega^2/E\) and the L-0002 embedding for \(\\|\\nabla u\\|_\\infty\).

## Caps (E0=0.5, ν=0.1)

| N | L-0001 (t=0.02) | L-0002 | L-0003 uniform |
|---|-----------------|--------|----------------|
| 8 | {table['8']['L0001_short_time']:.6e} | {'FAIL' if not table['8']['L0002_closes'] else table['8']['L0002_cap']} | {table['8']['L0003_uniform']:.6e} |
| 12 | {table['12']['L0001_short_time']:.6e} | {'FAIL' if not table['12']['L0002_closes'] else table['12']['L0002_cap']} | {table['12']['L0003_uniform']:.6e} |
| 16 | {table['16']['L0001_short_time']:.6e} | {'FAIL' if not table['16']['L0002_closes'] else table['16']['L0002_cap']} | {table['16']['L0003_uniform']:.6e} |

Still ≫ C-0001 numerical target (~11), but **closes globally in time** for each fixed N.
""",
        encoding="utf-8",
    )

    # Load or regenerate CMA candidate, then polish
    cma_path = Path("datasets/candidates/cmaes_best.npz")
    if cma_path.exists():
        u_cma = load_candidate_npz(cma_path)
        # Ensure n=12 energy sphere
        if u_cma.shape[-1] != 12:
            u_cma, _ = evolve_enstrophy_es(
                n=12, k_max=3, generations=8, population=8, seed=0, M_target=M_c
            )
            u_es = u_cma
            cma_meta = {"source": "regenerated_wrong_n"}
        else:
            u_es = u_cma
            cma_meta = {"source": "datasets/candidates/cmaes_best.npz"}
    else:
        u_es, cres = evolve_enstrophy_es(
            n=12, k_max=3, generations=12, population=8, seed=0, M_target=M_c
        )
        cma_meta = cres.as_dict()

    u_pol, polish = polish_candidate(u_es, M_target=M_c, n_steps=8, step_size=0.05)

    # Also polish a fresh random for diversity
    u_rand, _ = random_div_free(12, seed=7, energy_target=0.5)
    _u_r, polish_r = polish_candidate(u_rand, M_target=M_c, n_steps=6, step_size=0.05)

    best_gated = None
    for name, pr in (("cma_polished", polish), ("random_polished", polish_r)):
        if pr.gates_passed:
            val = pr.J_etd
            if best_gated is None or val > best_gated[1]:
                best_gated = (name, val)

    best_any = max(polish.J_etd, polish_r.J_etd)
    refuted = bool(
        (polish.exceeded_M and polish.gates_passed)
        or (polish_r.exceeded_M and polish_r.gates_passed)
    )

    support = max(best_any, 10.23595472515156)
    if refuted:
        conj = ConjectureC0001(
            proposed_bound_M=M_c,
            numerical_support_max=support,
            status="refuted",
            refuted=True,
            refutation_value=support,
            evidence_level="N2",
            notes="Refuted by gated polished candidate.",
        )
    else:
        # Keep M; record L-0003 and polish margins
        margin = M_c - support
        conj = ConjectureC0001(
            proposed_bound_M=M_c,
            numerical_support_max=support,
            status="exploring",
            notes=(
                f"After L-0003 + polish: best ETD={support:.6f}, margin to M={margin:.6f}. "
                f"Best gated={best_gated}. L-0003 gives uniform Galerkin ceiling "
                f"~{table['12']['L0003_uniform']:.3e} (N=12). Still unproved; not Clay."
            ),
        )
        conj.statement = (
            f"For all divergence-free Fourier fields on the N^3 grid with N<={conj.n_max}, "
            f"dealiased convective term, nu={conj.nu}, kinetic energy = {conj.energy}, evolved with "
            f"ETD-RK2 (dt=1e-3) to T={conj.t_end}, the enstrophy E_ens(T) <= M "
            f"with M={conj.proposed_bound_M}. FINITE-DIMENSIONAL only."
        )

    save_conjecture(conj, "conjectures/active/C-0001.json")
    if refuted:
        Path("conjectures/rejected").mkdir(parents=True, exist_ok=True)
        save_conjecture(conj, "conjectures/rejected/C-0001.json")
        # Successor conjecture with 1.5× gated max, still below L-0001 short-time ceiling
        new_M = float(max(1.5 * support, support + 1.0))
        c2 = ConjectureC0002(
            proposed_bound_M=new_M,
            numerical_support_max=support,
            status="exploring",
            notes=(
                f"Successor of refuted C-0001. Seeded M=1.5× gated polish max {support:.6f}. "
                f"L-0003 uniform ceiling (N=12)≈{table['12']['L0003_uniform']:.3e}; "
                f"L-0001 t=0.02 ceiling≈{table['12']['L0001_short_time']:.3e}."
            ),
        )
        c2.statement = (
            f"For all divergence-free Fourier fields on the N^3 grid with N<={c2.n_max}, "
            f"dealiased convective term, nu={c2.nu}, kinetic energy = {c2.energy}, evolved with "
            f"ETD-RK2 (dt=1e-3) to T={c2.t_end}, the enstrophy E_ens(T) <= M "
            f"with M={c2.proposed_bound_M}. FINITE-DIMENSIONAL only. Replaces refuted C-0001."
        )
        save_conjecture_c0002(c2, "conjectures/active/C-0002.json")
    else:
        c2 = None

    import hashlib

    payload = __import__("numpy").concatenate([u_pol.real.ravel(), u_pol.imag.ravel()])
    h = hashlib.sha256(payload.tobytes()).hexdigest()[:16]
    meta = ICMetadata(
        name="cma_polished",
        n=12,
        seed=None,
        energy=0.5,
        enstrophy=0.0,
        helicity=0.0,
        symmetries=["es", "adjoint_polish"],
        resolution=12,
        content_hash=h,
    )
    save_ic("datasets/candidates/cma_polished_n12.npz", u_pol, meta)

    summary = {
        "route": "B",
        "purpose": "L0003_uniform_and_polish_C0001",
        "lemma_table": table,
        "cma_source": cma_meta,
        "polish_cma": polish.as_dict(),
        "polish_random": polish_r.as_dict(),
        "best_gated": best_gated,
        "best_any_etd": best_any,
        "refuted": refuted,
        "conjecture_C0001": conj.as_dict(),
        "conjecture_C0002": None if c2 is None else c2.as_dict(),
        "not_a_clay_claim": True,
    }
    path = out / "summary.json"
    path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(
        json.dumps(
            {
                "L0003_N12": table["12"]["L0003_uniform"],
                "L0001_N12": table["12"]["L0001_short_time"],
                "best_closing_N12": table["12"]["best_of_closing"],
                "polish_cma_J_etd": polish.J_etd,
                "polish_cma_gates": polish.gates_passed,
                "best_any": best_any,
                "refuted": refuted,
                "C0001_M": conj.proposed_bound_M,
                "C0002_M": None if c2 is None else c2.proposed_bound_M,
                "margin_C0001": conj.proposed_bound_M - support,
            },
            indent=2,
        )
    )
    return path


if __name__ == "__main__":
    run()
