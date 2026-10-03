"""Sprint: Option A — all-IC C-0007 at L-0024 majorant M≈41.743 (L-0070)."""

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
from ns_exploration.conjectures.l0069_laptop_onepol_closure import (
    lemma_l0069,
    save_lemma_l0069,
)
from ns_exploration.conjectures.l0070_c0007_all_ic_l0024 import (
    lemma_l0070,
    save_lemma_l0070,
    write_c0007_proved,
    write_c0008_sharp,
)
from ns_exploration.validation.l0070_certificate import (
    build_l0070_certificate,
    save_l0070_certificate,
)


def main() -> dict:
    l70 = lemma_l0070()
    save_lemma_l0070(l70)
    write_c0007_proved(l70)
    write_c0008_sharp(l70)

    cert = build_l0070_certificate()
    save_l0070_certificate(cert)

    l69 = lemma_l0069()
    save_lemma_l0069(l69)

    l56 = lemma_l0056()
    save_lemma_l0056(l56)
    write_closure_map_report(l56)

    tests = subprocess.run(
        [
            sys.executable,
            "-m",
            "pytest",
            "python/ns_exploration/tests/test_l0070.py",
            "python/ns_exploration/tests/test_l0056.py",
            "python/ns_exploration/tests/test_l0024.py",
            "-q",
        ],
        cwd=Path.cwd(),
        env={**dict(__import__("os").environ), "PYTHONPATH": "python"},
        capture_output=True,
        text=True,
    )

    out = {
        "lemma_id": "L-0070",
        "proved_bound_M": l70.proved_bound_M,
        "sharp_target_M": l70.sharp_target_M,
        "closes_c0007_all_ic": l70.closes_c0007_all_ic,
        "closes_c0008_sharp": l70.closes_c0008_sharp,
        "c0007_status": "proved_restricted",
        "c0008_status": "exploring",
        "certificate": cert.certificate_id,
        "pytest_ok": tests.returncode == 0,
    }
    exp_dir = Path("experiments/exploratory/sprint02_l0070")
    exp_dir.mkdir(parents=True, exist_ok=True)
    (exp_dir / "summary.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))
    if tests.returncode != 0:
        print(tests.stdout, tests.stderr)
        raise RuntimeError("pytest failed")
    return out


if __name__ == "__main__":
    main()
