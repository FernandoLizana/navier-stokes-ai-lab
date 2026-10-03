"""Sprint: L-0043 one-pol Sym + C-R-0008 + shell-block/Ky-Fan techos."""

from __future__ import annotations

import json
from pathlib import Path

from ns_exploration.conjectures.l0043_onepol_shellblock import (
    lemma_l0043,
    save_cr0008,
    save_lemma_l0043,
)
from ns_exploration.validation.l0043_certificate import (
    build_l0043_certificate,
    save_l0043_certificate,
    verify_l0043_certificate,
)


def main() -> dict:
    b = lemma_l0043(n=24)
    save_lemma_l0043(b)
    save_cr0008(b)
    cert = build_l0043_certificate(n=24)
    ok, checks = verify_l0043_certificate(cert.as_dict())
    save_l0043_certificate(cert)
    if not ok:
        raise RuntimeError(f"L-0043 certificate failed: {checks}")

    c7_path = Path("conjectures/active/C-0007.json")
    c7 = json.loads(c7_path.read_text(encoding="utf-8"))
    c7["notes"] = (
        f"Open all-IC. Closed: C-R-0002..0008. L-0043: one-pol {{1,2,3,4,5,6,8}} "
        f"C={b.support_cr0008_C:.4f}≤C_† (C-R-0008); two-pol {{1..6}} still open "
        f"(shellblock_lo≳{b.C_shellblock_lo:.2f}, KyFan4≳{b.C_kyfan_partial:.2f}). "
        "Need SOS/rank UB for full two-pol consecutive band."
    )
    related = set(c7.get("related") or [])
    related.update(["L-0043", "C-R-0008"])
    c7["related"] = sorted(related)
    c7_path.write_text(json.dumps(c7, indent=2), encoding="utf-8")

    out = {
        "support_cr0008": b.support_cr0008,
        "support_cr0008_C": b.support_cr0008_C,
        "onepol_123456_C": b.onepol_123456_C,
        "onepol_techo_12345689_C": b.onepol_techo_12345689_C,
        "C_shellblock_lo": b.C_shellblock_lo,
        "C_kyfan_partial": b.C_kyfan_partial,
        "C_dagger": b.C_dagger,
        "closes_cr0008": True,
        "closes_c0007_all_ic": False,
        "certificate_ok": ok,
        "checks": {k: bool(v) for k, v in checks.items()},
        "not_a_clay_claim": True,
    }
    Path("reports").mkdir(exist_ok=True)
    Path("experiments/exploratory/sprint02_l0043").mkdir(parents=True, exist_ok=True)
    Path("experiments/exploratory/sprint02_l0043/summary.json").write_text(
        json.dumps(out, indent=2), encoding="utf-8"
    )
    Path("reports/SPRINT_02_L0043.md").write_text(
        f"""# Sprint — L-0043 one-pol Sym + C-R-0008

| Item | Value |
|------|-------|
| **C-R-0008** (one-pol shells {list(b.support_cr0008 or [])}) | C≈{b.support_cr0008_C:.6f} ≤ C_† |
| One-pol {{1..6}} | C≈{b.onepol_123456_C:.6f} ≤ C_† |
| One-pol + shell 9 (techo) | C≈{b.onepol_techo_12345689_C:.6f} > C_† |
| Shell-block lo (two-pol {{1..6}}) | ≳{b.C_shellblock_lo:.6f} > C_† |
| Ky-Fan 4-term (two-pol {{1..6}}) | ≳{b.C_kyfan_partial:.6f} > C_† |
| Two-pol {{1..6}} / C-0007 all-IC | still open |

FINITE Galerkin only. Not continuum. Not Clay.
""",
        encoding="utf-8",
    )
    print(json.dumps(out, indent=2, ensure_ascii=True))
    return out


if __name__ == "__main__":
    main()
