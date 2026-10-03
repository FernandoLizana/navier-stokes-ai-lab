"""Sprint: L-0049 rank-1 lower bound on full-dealias physical cubic."""

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

from ns_exploration.conjectures.l0049_rank1_fullsym import run_l0049, save_lemma_l0049


def main() -> dict:
    print(f"threads={_ns_threads}", flush=True)
    try:
        import psutil

        psutil.Process().nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
    except Exception:
        pass

    b = run_l0049(n=24, n_starts=64, n_iters=80, rebuild_M=False)
    save_lemma_l0049(b)

    out = b.as_dict()
    Path("reports").mkdir(exist_ok=True)
    Path("experiments/exploratory/sprint02_l0049").mkdir(parents=True, exist_ok=True)
    Path("experiments/exploratory/sprint02_l0049/summary.json").write_text(
        json.dumps(out, indent=2), encoding="utf-8"
    )
    Path("reports/l0049_rank1_full.json").write_text(
        json.dumps(out, indent=2), encoding="utf-8"
    )
    Path("reports/SPRINT_02_L0049.md").write_text(
        f"""# Sprint — L-0049 rank-1 full-dealias lower bound

| Item | Value |
|------|-------|
| C_rank1_lower | {b.C_rank1_lower:.6f} |
| C_fft at best z | {b.C_fft_at_best_z:.6f} |
| C_fullsym (L-0048) | {b.C_fullsym_full:.6f} |
| C_† | {b.C_dagger:.6f} |
| gap fullsym/rank1 | {b.gap_fullsym_over_rank1:.2f} |
| exceeds C_† (refute?) | {b.refutes_c0007} |

FINITE Galerkin only. Not continuum. Not Clay.
""",
        encoding="utf-8",
    )
    print("FINAL", json.dumps({k: out[k] for k in out if k != "notes"}, indent=2), flush=True)
    print("NOTES", b.notes, flush=True)
    return out


if __name__ == "__main__":
    main()
