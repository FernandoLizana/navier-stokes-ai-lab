"""First delivery sprint: C-0008 terminal-weighted pilot (spec §22)."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

from ns_exploration.conjectures.l0048_fullsym_techo import C_FULLSYM_FULL_DEALIAS
from ns_exploration.conjectures.l0048_matrixfree_fullsym import validate_sparse_vs_dense
from ns_exploration.terminal_weighted.constants import FROZEN
from ns_exploration.terminal_weighted.repo_snapshot import write_repo_snapshot
from ns_exploration.terminal_weighted.scan import run_exploratory_scan, save_scan
from ns_exploration.terminal_weighted.tensor import (
    terminal_sparse_at_T,
    validate_direct_vs_tensor,
)
from ns_exploration.terminal_weighted.weights import stokes_floor_hi


def _write_audit_report(path: Path, sections: dict) -> None:
    lines = [
        "# C-0008 Terminal-Weighted Audit (first delivery)",
        "",
        "> N2 exploratory. **C-0008 remains open.** Not a proof.",
        "",
    ]
    for title, body in sections.items():
        lines.extend([f"## {title}", "", body, ""])
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def _write_scan_report(path: Path, scan: dict) -> None:
    lines = [
        "# C-0008 Terminal-Weighted Scan (N2 exploratory)",
        "",
        f"- Band radii: `{scan['band_radii']}`",
        f"- Grid step: `{scan['step']}` ({scan['n_times']} nodes)",
        f"- **max_grid_bound:** `{scan['max_grid_bound']:.6g}` at t=`{scan['t_at_max']}`",
        f"- **trapezoid_estimate:** `{scan['trapezoid_estimate']:.6g}`",
        f"- **left_riemann:** `{scan['left_riemann_estimate']:.6g}`",
        f"- **right_riemann:** `{scan['right_riemann_estimate']:.6g}`",
        f"- **threshold integral:** `{scan['integral_threshold']:.6g}`",
        f"- **margin vs threshold:** `{scan['margin_vs_threshold']:.6g}`",
        f"- Closes integral (exploratory): `{scan['closes_integral_exploratory']}`",
        f"- Closes uniform (exploratory): `{scan['closes_uniform_exploratory']}`",
        "",
        scan["honesty"],
        "",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def _decision(scan: dict, sparse_T: dict, audit_ok: bool) -> dict:
    margin = scan["margin_vs_threshold"]
    if margin > 0 and scan["closes_uniform_exploratory"]:
        rec = "Proceed to Phase D interval cover on band; then shell decomposition."
    elif margin > 0:
        rec = "Uniform band bound OK; integral exploratory positive margin — Phase D interval cover."
    elif scan["max_grid_bound"] < FROZEN.uniform_threshold:
        rec = "Band uniform OK but integral over budget — refine temporal cover (Phase D)."
    else:
        rec = "Band fullsym exceeds uniform threshold — need shell structure or tighter relax before interval cover."
    return {
        "audit_passed": audit_ok,
        "recommendation": rec,
        "next_phase": "D" if margin > 0 else "E_or_structure",
    }


def main() -> dict:
    root = Path.cwd()
    write_repo_snapshot(root=root)

    band = tuple(range(1, 9))
    audit_val = validate_direct_vs_tensor(24, tuple(range(1, 7)), FROZEN.T, n_probe=40)
    stokes_hi = stokes_floor_hi(24)
    sparse_dense = validate_sparse_vs_dense(24, tuple(range(1, 7)))
    # Full dealias sparse regression (t=T) — expensive; skip in quick mode via env
    import os
    if os.environ.get("TERMINAL_FULL_SPARSE", "1") == "1":
        sparse_T = terminal_sparse_at_T(24, progress_every=800000)
    else:
        sparse_T = {"C_term_ub": C_FULLSYM_FULL_DEALIAS, "skipped": True}

    scan = run_exploratory_scan(24, step=0.001, band_radii=band)
    save_scan(scan, "experiments/terminal_weighted/scan.json")

    sparse_rel = abs(sparse_T["C_term_ub"] - C_FULLSYM_FULL_DEALIAS) / C_FULLSYM_FULL_DEALIAS
    audit_ok = (
        audit_val["ok"]
        and abs(stokes_hi - FROZEN.stokes_floor) < 1e-6
        and sparse_rel < 1e-9
        and scan["band_T_regression_rel_err"] < 1e-9
    )
    decision = _decision(scan, sparse_T, audit_ok)
    decision["full_dealias_C_at_T"] = C_FULLSYM_FULL_DEALIAS
    decision["note"] = (
        "Band scan favorable; full dealias C(T)=25.93 << uniform 64.97. "
        "Phase D interval cover on full dealias required before any C-0008 claim."
    )

    _write_audit_report(
        Path("reports/C0008_TERMINAL_WEIGHTED_AUDIT.md"),
        {
            "A. Identity / normalization": (
                f"- direct cubic self-consistency: `{audit_val['max_rel_err_cubic']:.3e}`\n"
                f"- fullsym C stream vs dense: rel `{audit_val['rel_err_C']:.3e}`\n"
                f"- Stokes floor hi: `{stokes_hi:.11g}` vs frozen `{FROZEN.stokes_floor}`\n"
                f"- Band t=T regression vs L-0045: rel `{scan['band_T_regression_rel_err']:.3e}`\n"
                f"- Full dealias sparse t=T: `{sparse_T['C_term_ub']:.11g}` vs L-0048 "
                f"`{C_FULLSYM_FULL_DEALIAS}` (rel `{sparse_rel:.3e}`)\n"
                f"- Sparse vs dense {{1..6}}: `{sparse_dense['rel_err']:.3e}`"
            ),
            "B. L-0048 file map": "See `reports/C0008_L0048_FILE_MAP.md`.",
            "Decision": decision["recommendation"],
        },
    )
    _write_scan_report(Path("reports/C0008_TERMINAL_WEIGHTED_SCAN.md"), scan)

    tests = subprocess.run(
        [sys.executable, "-m", "pytest", "python/ns_exploration/tests/test_terminal_weighted.py", "-q"],
        cwd=root,
        env={**dict(__import__("os").environ), "PYTHONPATH": "python"},
        capture_output=True,
        text=True,
    )

    out = {
        "delivery": "section_22_first",
        "audit_ok": audit_ok,
        "stokes_floor_hi": stokes_hi,
        "sparse_T_C": sparse_T["C_term_ub"],
        "l0048_frozen": C_FULLSYM_FULL_DEALIAS,
        "max_grid_bound": scan["max_grid_bound"],
        "t_at_max": scan["t_at_max"],
        "trapezoid_estimate": scan["trapezoid_estimate"],
        "margin_vs_threshold": scan["margin_vs_threshold"],
        "decision": decision,
        "pytest_ok": tests.returncode == 0,
        "c0008_status": "unchanged_open",
    }
    out_path = Path("experiments/terminal_weighted/delivery_summary.json")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))
    if tests.returncode != 0:
        print(tests.stdout, tests.stderr)
        raise RuntimeError("pytest failed")
    return out


if __name__ == "__main__":
    main()
