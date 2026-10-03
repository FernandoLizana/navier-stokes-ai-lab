"""Sprint: L-0042 Sym-strengthening techo on {1..6}."""

from __future__ import annotations

import json
from pathlib import Path

from ns_exploration.conjectures.l0042_sym_strengthen_techo import (
    lemma_l0042,
    save_lemma_l0042,
)
from ns_exploration.validation.l0042_certificate import (
    build_l0042_certificate,
    save_l0042_certificate,
    verify_l0042_certificate,
)


def main() -> dict:
    b = lemma_l0042(n=24)
    save_lemma_l0042(b)
    cert = build_l0042_certificate(n=24)
    ok, checks = verify_l0042_certificate(cert.as_dict())
    save_l0042_certificate(cert)
    if not ok:
        raise RuntimeError(f"L-0042 certificate failed: {checks}")

    c7_path = Path("conjectures/active/C-0007.json")
    c7 = json.loads(c7_path.read_text(encoding="utf-8"))
    c7["notes"] = (
        f"Open all-IC. Closed: C-R-0002..0007. L-0042: {{1..6}} Sym={b.C_sym_123456:.4f}>C_†; "
        f"triangle/hybrid worse; rank1_lo={b.C_rank1_lo_123456:.4f}. Need SOS/rank upper bound."
    )
    related = set(c7.get("related") or [])
    related.add("L-0042")
    c7["related"] = sorted(related)
    c7_path.write_text(json.dumps(c7, indent=2), encoding="utf-8")

    out = {
        "C_sym_12345": b.C_sym_12345,
        "C_sym_123456": b.C_sym_123456,
        "C_triangle_123456": b.C_triangle_123456,
        "C_hybrid_AB": b.C_hybrid_AB,
        "C_rank1_lo_123456": b.C_rank1_lo_123456,
        "C_dagger": b.C_dagger,
        "closes_123456": False,
        "closes_c0007_all_ic": False,
        "certificate_ok": ok,
        "checks": {k: bool(v) for k, v in checks.items()},
        "not_a_clay_claim": True,
    }
    Path("reports").mkdir(exist_ok=True)
    Path("experiments/exploratory/sprint02_l0042").mkdir(parents=True, exist_ok=True)
    Path("experiments/exploratory/sprint02_l0042/summary.json").write_text(
        json.dumps(out, indent=2), encoding="utf-8"
    )
    Path("reports/SPRINT_02_L0042.md").write_text(
        f"""# Sprint — L-0042 Sym-strengthening techo

| Item | Value |
|------|-------|
| C_Sym {{1..5}} | ≈{b.C_sym_12345:.6f} ≤ C_† |
| C_Sym {{1..6}} | ≈{b.C_sym_123456:.6f} > C_† |
| Column triangle | ≈{b.C_triangle_123456:.6f} |
| A/B hybrid | ≈{b.C_hybrid_AB:.6f} |
| Rank-1 lower | ≈{b.C_rank1_lo_123456:.6f} |
| {{1..6}} / C-0007 | still open |

FINITE Galerkin only. Not continuum. Not Clay.
""",
        encoding="utf-8",
    )
    print(json.dumps(out, indent=2, ensure_ascii=True))
    return out


if __name__ == "__main__":
    main()
