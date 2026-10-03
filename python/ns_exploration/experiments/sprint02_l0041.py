"""Sprint: L-0041 shell-6 Sym + C-R-0007."""

from __future__ import annotations

import json
from pathlib import Path

from ns_exploration.conjectures.l0041_shell6_sym import (
    build_cr0007,
    lemma_l0041,
    save_cr0007,
    save_lemma_l0041,
)
from ns_exploration.validation.l0041_certificate import (
    build_l0041_certificate,
    save_l0041_certificate,
    verify_l0041_certificate,
)


def main() -> dict:
    b = lemma_l0041(n=24)
    save_lemma_l0041(b)
    cr = build_cr0007(n=24)
    save_cr0007(cr)
    cert = build_l0041_certificate(n=24)
    ok, checks = verify_l0041_certificate(cert.as_dict())
    save_l0041_certificate(cert)
    if not ok:
        raise RuntimeError(f"L-0041 certificate failed: {checks}")

    c7_path = Path("conjectures/active/C-0007.json")
    c7 = json.loads(c7_path.read_text(encoding="utf-8"))
    c7["notes"] = (
        f"Open all-IC. Closed: C-R-0002..0007. L-0041: {{1..6}} Sym C={b.band_123456_C:.4f}>C_† "
        f"(techo); C-R-0007 support with shell 6, C={b.support_cr0007_C:.4f}≤C_†."
    )
    related = set(c7.get("related") or [])
    related.update(["L-0041", "C-R-0007"])
    c7["related"] = sorted(related)
    c7_path.write_text(json.dumps(c7, indent=2), encoding="utf-8")

    out = {
        "band_123456_C": b.band_123456_C,
        "seed_12346_C": b.seed_12346_C,
        "seed_23456_C": b.seed_23456_C,
        "seed_13456_C": b.seed_13456_C,
        "support_cr0007_C": b.support_cr0007_C,
        "support_cr0007_n": len(b.support_cr0007 or []),
        "C_dagger": b.C_dagger,
        "closes_cr0007": True,
        "consecutive_6_techo": True,
        "closes_c0007_all_ic": False,
        "certificate_ok": ok,
        "checks": {k: bool(v) for k, v in checks.items()},
        "not_a_clay_claim": True,
    }
    Path("reports").mkdir(exist_ok=True)
    Path("experiments/exploratory/sprint02_l0041").mkdir(parents=True, exist_ok=True)
    Path("experiments/exploratory/sprint02_l0041/summary.json").write_text(
        json.dumps(out, indent=2), encoding="utf-8"
    )
    Path("reports/SPRINT_02_L0041.md").write_text(
        f"""# Sprint — L-0041 shell-6 Sym + C-R-0007

| Item | Value |
|------|-------|
| C_Sym {{1..6}} | ≈{b.band_123456_C:.6f} > C_† (techo) |
| Seed {{1,2,3,4,6}} | ≈{b.seed_12346_C:.6f} |
| **C-R-0007** ({len(b.support_cr0007 or [])} shells) | C≈{b.support_cr0007_C:.6f} ≤ C_† |
| C-0007 all-IC | still open |

FINITE Galerkin only. Not continuum. Not Clay.
""",
        encoding="utf-8",
    )
    print(json.dumps(out, indent=2, ensure_ascii=True))
    return out


if __name__ == "__main__":
    main()
