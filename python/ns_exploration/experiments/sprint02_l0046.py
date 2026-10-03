"""Sprint: L-0046 greedy physical fullsym + C-R-0011."""

from __future__ import annotations

import json
from pathlib import Path

from ns_exploration.conjectures.l0046_greedy_fullsym import (
    lemma_l0046,
    save_cr0011,
    save_lemma_l0046,
)
from ns_exploration.validation.l0046_certificate import (
    build_l0046_certificate,
    save_l0046_certificate,
    verify_l0046_certificate,
)


def main() -> dict:
    b = lemma_l0046(n=24)
    save_lemma_l0046(b)
    save_cr0011(b)
    cert = build_l0046_certificate(n=24, bound=b)
    ok, checks = verify_l0046_certificate(cert.as_dict())
    save_l0046_certificate(cert)
    if not ok:
        raise RuntimeError(f"L-0046 certificate failed: {checks}")

    c7_path = Path("conjectures/active/C-0007.json")
    c7 = json.loads(c7_path.read_text(encoding="utf-8"))
    c7["notes"] = (
        f"Open all-IC. Closed: C-R-0002..0011. L-0046: greedy fullsym "
        f"{len(b.support_cr0011 or [])} shells (incl. 147) C={b.C_fullsym_cr0011:.4f}≤C_†. "
        "Full dealias still open."
    )
    related = set(c7.get("related") or [])
    related.update(["L-0046", "C-R-0011"])
    c7["related"] = sorted(related)
    c7_path.write_text(json.dumps(c7, indent=2), encoding="utf-8")

    out = {
        "support_cr0011": b.support_cr0011,
        "C_fullsym_cr0011": b.C_fullsym_cr0011,
        "C_fullsym_cr0010": b.C_fullsym_cr0010,
        "C_dagger": b.C_dagger,
        "closes_cr0011": True,
        "closes_c0007_all_ic": False,
        "certificate_ok": ok,
        "checks": {k: bool(v) for k, v in checks.items()},
        "not_a_clay_claim": True,
    }
    Path("reports").mkdir(exist_ok=True)
    Path("experiments/exploratory/sprint02_l0046").mkdir(parents=True, exist_ok=True)
    Path("experiments/exploratory/sprint02_l0046/summary.json").write_text(
        json.dumps(out, indent=2), encoding="utf-8"
    )
    Path("reports/SPRINT_02_L0046.md").write_text(
        f"""# Sprint — L-0046 greedy physical fullsym + C-R-0011

| Item | Value |
|------|-------|
| **C-R-0011** ({len(b.support_cr0011 or [])} shells, incl. 147) | C≈{b.C_fullsym_cr0011:.6f} ≤ C_† |
| C-R-0010 baseline | C≈{b.C_fullsym_cr0010:.6f} |
| C-0007 all-IC | still open |

FINITE Galerkin only. Not continuum. Not Clay.
""",
        encoding="utf-8",
    )
    print(json.dumps(out, indent=2, ensure_ascii=True))
    return out


if __name__ == "__main__":
    main()
