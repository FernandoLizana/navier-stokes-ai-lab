"""Sprint: L-0048 sparse fullsym attack on full-dealias / all-IC C-0007."""

from __future__ import annotations

import json
import os
from pathlib import Path

_default_threads = max(1, (os.cpu_count() or 4) - 2)
_ns_threads = os.environ.get("NS_THREADS", str(_default_threads))
for _var in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
):
    os.environ.setdefault(_var, _ns_threads)

from ns_exploration.conjectures.l0027_low_slab_cubic import lemma_l0027
from ns_exploration.conjectures.l0048_matrixfree_fullsym import (
    all_dealias_radii,
    sparse_C_fullsym,
    validate_sparse_vs_dense,
)


def main() -> dict:
    print(f"threads={_ns_threads}", flush=True)
    Cd = lemma_l0027(n=24, empirical=False).C_dagger
    print("validate {1..6}...", flush=True)
    val = validate_sparse_vs_dense()
    print(json.dumps(val, indent=2), flush=True)
    if not val["ok"]:
        raise RuntimeError(f"sparse≠dense: {val}")

    radii = all_dealias_radii(24)
    print(f"full dealias shells={len(radii)} -> sparse C...", flush=True)
    full = sparse_C_fullsym(24, radii, progress_every=500_000)
    print(json.dumps({k: full[k] for k in full if k != "radii"}, indent=2), flush=True)

    closes = full["C_fullsym"] <= Cd + 1e-9
    out = {
        "C_dagger": Cd,
        "validate": val,
        "full_dealias": {k: full[k] for k in full if k != "radii"},
        "closes_c0007_all_ic": closes,
        "n_shells": len(radii),
        "not_a_clay_claim": True,
        "evidence_level": "N7" if closes else "N2",
        "notes": (
            "Sparse CSR Gram of physical fullsym on all dealias shells. "
            + ("C≤C_† ⇒ all-IC via L-0027+L-0026." if closes else "If C>C_†: fullsym Shor techo for all-IC.")
        ),
    }
    Path("reports").mkdir(exist_ok=True)
    Path("experiments/exploratory/sprint02_l0048").mkdir(parents=True, exist_ok=True)
    Path("experiments/exploratory/sprint02_l0048/summary.json").write_text(
        json.dumps(out, indent=2), encoding="utf-8"
    )
    Path("reports/l0048_full_dealias_sparse.json").write_text(
        json.dumps(out, indent=2), encoding="utf-8"
    )
    Path("reports/SPRINT_02_L0048.md").write_text(
        f"""# Sprint — L-0048 sparse fullsym full-dealias

| Item | Value |
|------|-------|
| validate {{1..6}} rel err | {val['rel_err']:.3e} |
| full dealias D | {full['D']} |
| nnz | {full['nnz']} |
| C_fullsym | {full['C_fullsym']:.6f} |
| C_† | {Cd:.6f} |
| closes C-0007 all-IC | {closes} |

FINITE Galerkin only. Not continuum. Not Clay.
""",
        encoding="utf-8",
    )
    print("FINAL closes_all_ic", closes, "C", full["C_fullsym"], "Cd", Cd, flush=True)
    return out


if __name__ == "__main__":
    main()
