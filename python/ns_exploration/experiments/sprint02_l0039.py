"""Sprint: L-0039 sparse support Shor + C-R-0005."""

from __future__ import annotations

import json
from pathlib import Path

from ns_exploration.conjectures.l0039_sparse_support_shor import (
    build_cr0005,
    lemma_l0039,
    save_cr0005,
    save_lemma_l0039,
)
from ns_exploration.validation.l0039_certificate import (
    build_l0039_certificate,
    save_l0039_certificate,
    verify_l0039_certificate,
)


def main() -> dict:
    b = lemma_l0039(n=24)
    save_lemma_l0039(b)
    cr = build_cr0005(n=24)
    save_cr0005(cr)
    cert = build_l0039_certificate(n=24)
    ok, checks = verify_l0039_certificate(cert.as_dict())
    save_l0039_certificate(cert)
    if not ok:
        raise RuntimeError(f"L-0039 certificate failed: {checks}")

    c7_path = Path("conjectures/active/C-0007.json")
    c7 = json.loads(c7_path.read_text(encoding="utf-8"))
    c7["notes"] = (
        f"Open all-IC. Closed: C-R-0002/0003/0004/0005. L-0039: support A "
        f"({len(b.support_A or [])} shells) C_Shor={b.support_A_C:.4f}≤C_†. "
        f"Need beyond sparse Shor for remaining shells."
    )
    related = set(c7.get("related") or [])
    related.update(["L-0039", "C-R-0005"])
    c7["related"] = sorted(related)
    c7_path.write_text(json.dumps(c7, indent=2), encoding="utf-8")

    out = {
        "quad_1346_C": b.quad_1346_C,
        "support_A_C": b.support_A_C,
        "support_A_n": len(b.support_A or []),
        "support_B_C": b.support_B_C,
        "support_B_n": len(b.support_B or []),
        "C_dagger": b.C_dagger,
        "closes_cr0005": True,
        "closes_c0007_all_ic": False,
        "certificate_ok": ok,
        "checks": {k: bool(v) for k, v in checks.items()},
        "not_a_clay_claim": True,
    }
    Path("reports").mkdir(exist_ok=True)
    Path("experiments/exploratory/sprint02_l0039").mkdir(parents=True, exist_ok=True)
    Path("experiments/exploratory/sprint02_l0039/summary.json").write_text(
        json.dumps(out, indent=2), encoding="utf-8"
    )
    Path("reports/SPRINT_02_L0039.md").write_text(
        f"""# Sprint — L-0039 sparse support Shor + C-R-0005

| Item | Value |
|------|-------|
| Quad {{1,3,4,6}} C_Shor | ≈{b.quad_1346_C:.6f} |
| Support A ({len(b.support_A or [])} shells) | C≈{b.support_A_C:.6f} ≤ C_† |
| Support B ({len(b.support_B or [])} shells) | C≈{b.support_B_C:.6f} ≤ C_† |
| **C-R-0005** | support A → Ω(T)≤M (proved) |
| C-0007 all-IC | still open |

FINITE Galerkin only. Not continuum. Not Clay.
""",
        encoding="utf-8",
    )
    print(json.dumps(out, indent=2, ensure_ascii=True))
    return out


if __name__ == "__main__":
    main()
