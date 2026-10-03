"""Final laptop closure sprint: L-0069 + refresh maps + tests."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

from ns_exploration.conjectures.l0056_closure_map import (
    lemma_l0056,
    save_lemma_l0056,
    write_closure_map_report,
)
from ns_exploration.conjectures.l0063_greedy_cluster_export import lemma_l0063
from ns_exploration.conjectures.l0069_laptop_onepol_closure import (
    lemma_l0069,
    patch_l0064_status,
    save_lemma_l0069,
)


def main() -> dict:
    patch_l0064_status()
    l69 = lemma_l0069()
    save_lemma_l0069(l69)

    l56 = lemma_l0056()
    save_lemma_l0056(l56)
    write_closure_map_report(l56)

    l63 = lemma_l0063()
    c7 = json.loads(Path("conjectures/active/C-0007.json").read_text(encoding="utf-8"))
    rel = set(c7.get("related") or [])
    rel.update(["L-0069"])
    c7["related"] = sorted(rel)
    c7["notes"] = (
        f"Open all-IC. Laptop path closed: C-R-0002..0014 ({l69.n_proved_restricted}). "
        f"Best one-pol C-R-0014 D={l69.best_onepol_D} SOS≈{l69.best_onepol_C_ub_sos_hi:.4f}. "
        f"One-pol techo shell {l69.techo_shell}. Cluster greedy export={l63.status}."
    )
    Path("conjectures/active/C-0007.json").write_text(json.dumps(c7, indent=2), encoding="utf-8")

    tests = subprocess.run(
        [
            sys.executable,
            "-m",
            "pytest",
            "python/ns_exploration/tests/test_l0064.py",
            "python/ns_exploration/tests/test_l0065.py",
            "python/ns_exploration/tests/test_l0066.py",
            "python/ns_exploration/tests/test_l0067.py",
            "-q",
        ],
        cwd=Path.cwd(),
        env={**dict(__import__("os").environ), "PYTHONPATH": "python"},
        capture_output=True,
        text=True,
    )

    out = {
        "lemma_id": "L-0069",
        "n_proved_restricted": l69.n_proved_restricted,
        "best_onepol_cr": l69.best_onepol_cr,
        "techo_shell": l69.techo_shell,
        "cluster_greedy": l63.status,
        "pytest_ok": tests.returncode == 0,
        "all_ic_open": True,
    }
    Path("experiments/exploratory/sprint02_l0069").mkdir(parents=True, exist_ok=True)
    Path("experiments/exploratory/sprint02_l0069/summary.json").write_text(
        json.dumps(out, indent=2), encoding="utf-8"
    )
    print(json.dumps(out, indent=2))
    if tests.returncode != 0:
        print(tests.stdout, tests.stderr)
        raise RuntimeError("pytest failed")
    return out


if __name__ == "__main__":
    main()
