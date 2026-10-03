"""Run L-0076 N<=24 inclusion empirical audit."""

from __future__ import annotations

import json
from pathlib import Path

from ns_exploration.terminal_weighted.n_coverage import full_dealias_inclusion_audit

REPO = Path(__file__).resolve().parents[2]


def main() -> dict:
    rep = full_dealias_inclusion_audit(prec=128)
    out_json = REPO / "reports" / "L0076_N_INCLUSION_AUDIT.json"
    out_md = REPO / "reports" / "L0076_N_INCLUSION_AUDIT.md"
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(rep, indent=2), encoding="utf-8")

    r24 = rep.get("repair_full_dealias_n24", {})
    md = f"""# L-0076 N≤24 Inclusion Audit

**Evidence level:** {rep['evidence_level']}  
**C-0008:** {rep['conjecture_C0008']}

## Wavevector embedding (n=24 majorant)

- All subsets: `{rep['wavevector_domination']['all_subset']}`

## Band majorant ladders

| Band | Monotone (n=12→24) |
|------|---------------------|
| {{1,2,3}} | {rep.get('band_majorant_ladder_3', {}).get('monotone_non_decreasing')} |
| {{1..6}} | {rep['band_majorant_ladder_6']['monotone_non_decreasing']} |

## Repair full-dealias n=24 (Route A)

| Field | Value |
|-------|-------|
| Shells | {r24.get('n_shells')} |
| Complete | {r24.get('repair_complete')} |
| I_hi | {r24.get('I_hi_route_A')} |

## Conditional conclusion

{rep['conditional_conclusion']}

## Formal gaps (still open)

"""
    for g in rep["formal_gap"]:
        md += f"- {g}\n"

    out_md.write_text(md, encoding="utf-8")
    print(json.dumps({"json": str(out_json), "md": str(out_md), **rep}, indent=2))
    return rep


if __name__ == "__main__":
    main()
