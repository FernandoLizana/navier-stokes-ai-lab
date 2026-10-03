"""Sprint: L-0038 band SVD Shor + C-R-0004."""

from __future__ import annotations

import json
from pathlib import Path

from ns_exploration.conjectures.l0038_band_svd_shor import (
    build_cr0004,
    lemma_l0038,
    save_cr0004,
    save_lemma_l0038,
)
from ns_exploration.validation.l0038_certificate import (
    build_l0038_certificate,
    save_l0038_certificate,
    verify_l0038_certificate,
)


def main() -> dict:
    b = lemma_l0038(n=24)
    save_lemma_l0038(b)
    cr = build_cr0004(n=24)
    save_cr0004(cr)
    cert = build_l0038_certificate(n=24)
    ok, checks = verify_l0038_certificate(cert.as_dict())
    save_l0038_certificate(cert)
    if not ok:
        raise RuntimeError(f"L-0038 certificate failed: {checks}")

    c7_path = Path("conjectures/active/C-0007.json")
    c7 = json.loads(c7_path.read_text(encoding="utf-8"))
    c7["notes"] = (
        f"Open all-IC. Closed: C-R-0002/0003/0004. L-0038: band {{1,2,3}} "
        f"C_Shor={b.band_123_C:.4f}≤C_† but {{1,2,3,4}} C={b.band_1234_C:.4f}>C_†. "
        f"Need multi-shell beyond low band."
    )
    related = set(c7.get("related") or [])
    related.update(["L-0038", "C-R-0004"])
    c7["related"] = sorted(related)
    c7_path.write_text(json.dumps(c7, indent=2), encoding="utf-8")

    out = {
        "band_12_C": b.band_12_C,
        "band_123_C": b.band_123_C,
        "band_125_C": b.band_125_C,
        "band_1234_C": b.band_1234_C,
        "C_dagger": b.C_dagger,
        "closes_cr0004": True,
        "closes_c0007_all_ic": False,
        "certificate_ok": ok,
        "checks": {k: bool(v) for k, v in checks.items()},
        "not_a_clay_claim": True,
    }
    Path("reports").mkdir(exist_ok=True)
    Path("experiments/exploratory/sprint02_l0038").mkdir(parents=True, exist_ok=True)
    Path("experiments/exploratory/sprint02_l0038/summary.json").write_text(
        json.dumps(out, indent=2), encoding="utf-8"
    )
    Path("reports/SPRINT_02_L0038.md").write_text(
        f"""# Sprint — L-0038 band SVD Shor + C-R-0004

| Item | Value |
|------|-------|
| C_Shor {{1,2}} | ≈{b.band_12_C:.6f} |
| C_Shor {{1,2,3}} | ≈{b.band_123_C:.6f} ≤ C_† |
| C_Shor {{1,2,5}} | ≈{b.band_125_C:.6f} ≤ C_† |
| C_Shor {{1,2,3,4}} | ≈{b.band_1234_C:.6f} > C_† |
| **C-R-0004** | shells {{1,2,3}} → Ω(T)≤M (proved) |
| C-0007 all-IC | still open |

FINITE Galerkin only. Not continuum. Not Clay.
""",
        encoding="utf-8",
    )
    print(json.dumps(out, indent=2, ensure_ascii=True))
    return out


if __name__ == "__main__":
    main()
