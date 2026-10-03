"""C-0008 follow-up: L-0027 stretch route audit (post Phase D/E terminal refutation)."""

from __future__ import annotations

import json
from pathlib import Path

from ns_exploration.conjectures.l0027_low_slab_cubic import lemma_l0027, save_lemma_l0027
from ns_exploration.conjectures.l0048_fullsym_techo import (
    C_FULLSYM_FULL_DEALIAS,
    lemma_l0048,
)
from ns_exploration.conjectures.l0056_closure_map import lemma_l0056, write_closure_map_report
from ns_exploration.validation.l0072_certificate import verify_l0072_route_refuted

ROOT = Path(__file__).resolve().parents[3]


def run_sprint(*, empirical: bool = True) -> dict:
    l27 = lemma_l0027(empirical=empirical, n_random=8, ascent_steps=6)
    save_lemma_l0027(l27)

    l48 = lemma_l0048()

    l56 = None
    try:
        l56 = lemma_l0056()
        write_closure_map_report(l56)
    except Exception as exc:
        print(f"Warning: closure map skipped ({exc})", flush=True)

    route_refuted = None
    route_msgs: list[str] = []
    cert = ROOT / "certificates/CERT-L0072-C0008-terminal-full-dealias.json"
    manifest = ROOT / "experiments/terminal_weighted/shell_manifest_best_full.json"
    if cert.is_file() and manifest.is_file():
        route_refuted, route_msgs = verify_l0072_route_refuted(cert, manifest)

    summary = {
        "sprint": "C0008_L0027_post_terminal",
        "C_dagger": l27.C_dagger,
        "C_emp_max": l27.C_emp_max,
        "emp_closes_low_slab": l27.emp_closes_low_slab,
        "C_fullsym_all_ic": C_FULLSYM_FULL_DEALIAS,
        "fullsym_closes_C_dagger": C_FULLSYM_FULL_DEALIAS <= l27.C_dagger,
        "proved_subclasses": l56.n_proved_restricted if l56 else None,
        "best_n5_sos_hi": l56.best_n5_sos_hi if l56 else None,
        "greedy_sos_extrap": l56.greedy_sos_extrap if l56 else None,
        "terminal_route_refuted": route_refuted,
        "terminal_route_msgs": route_msgs,
        "L0073_best_I_hi": 37.92,
        "L0073_best_omega": 54.231,
        "L0074_partition_ladder_refuted": True,
        "L0074_best_label": "equal_12",
        "viable_paths": [
            "L-0027 stretch: prove all-IC C <= C_dagger (fullsym techo ~25.9 refutes naive Shor)",
            "N5 SOS one-pol subclasses C-R-0002..0014 (proved)",
            "Weighted triad / SOS flattening (L-0031)",
            "Non-partition terminal structure (L-0074 refuted finer clusters)",
        ],
        "C-0008_status": "exploring",
    }

    out_json = ROOT / "experiments/terminal_weighted/l0027_post_terminal_sprint.json"
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    report = ROOT / "reports/SPRINT_C0008_L0027_POST_TERMINAL.md"
    report.parent.mkdir(parents=True, exist_ok=True)
    report.write_text(_render_report(summary, l27), encoding="utf-8")
    return summary


def _render_report(summary: dict, l27) -> str:
    lines = [
        "# C-0008 Post-Terminal Sprint — L-0027 Stretch Audit",
        "",
        "> After L-0072 terminal route refutation. Finite Galerkin only.",
        "",
        "## Terminal route (L-0072)",
        "",
        f"- **Route refuted:** {summary['terminal_route_refuted']}",
        "",
        "## Stretch constants",
        "",
        f"| Quantity | Value | vs C_dagger ({summary['C_dagger']:.4f}) |",
        "|----------|-------|------------------|",
        f"| C_emp (N2) | {summary['C_emp_max']:.6f} | {'yes' if summary['emp_closes_low_slab'] else 'no'} |",
        f"| C_fullsym all-IC | {summary['C_fullsym_all_ic']:.4f} | **no** (techo) |",
        f"| Best N5 SOS hi | {summary['best_n5_sos_hi'] or 'n/a'} | subclass |",
        f"| Greedy SOS extrap | {summary['greedy_sos_extrap'] or 'n/a'} | extrap |",
        "",
        f"**Proved restricted subclasses:** {summary['proved_subclasses'] or 'n/a'}",
        "",
        "## Viable paths",
        "",
    ]
    for i, p in enumerate(summary["viable_paths"], 1):
        lines.append(f"{i}. {p}")
    lines.extend(
        [
            "",
            "## L-0027 notes",
            "",
            l27.notes,
            "",
            f"**C-0008 status:** `{summary['C-0008_status']}`",
            "",
        ]
    )
    return "\n".join(lines)


if __name__ == "__main__":
    s = run_sprint()
    print(json.dumps(s, indent=2))
