"""Sprint: L-0040 Sym Shor + C-R-0006."""

from __future__ import annotations

import json
from pathlib import Path

from ns_exploration.conjectures.l0040_sym_shor import (
    build_cr0006,
    lemma_l0040,
    save_cr0006,
    save_lemma_l0040,
)
from ns_exploration.validation.l0040_certificate import (
    build_l0040_certificate,
    save_l0040_certificate,
    verify_l0040_certificate,
)


def main() -> dict:
    b = lemma_l0040(n=24)
    save_lemma_l0040(b)
    cr = build_cr0006(n=24)
    save_cr0006(cr)
    cert = build_l0040_certificate(n=24)
    ok, checks = verify_l0040_certificate(cert.as_dict())
    save_l0040_certificate(cert)
    if not ok:
        raise RuntimeError(f"L-0040 certificate failed: {checks}")

    c7_path = Path("conjectures/active/C-0007.json")
    c7 = json.loads(c7_path.read_text(encoding="utf-8"))
    c7["notes"] = (
        f"Open all-IC. Closed: C-R-0002..0006. L-0040 Sym Shor: {{1..5}} "
        f"C={b.band_12345_C:.4f}≤C_† but {{1..6}} C={b.band_123456_C:.4f}>C_†."
    )
    related = set(c7.get("related") or [])
    related.update(["L-0040", "C-R-0006"])
    c7["related"] = sorted(related)
    c7_path.write_text(json.dumps(c7, indent=2), encoding="utf-8")

    out = {
        "band_1234_C_sym": b.band_1234_C_sym,
        "band_12345_C": b.band_12345_C,
        "band_123456_C": b.band_123456_C,
        "greedy9_C": b.greedy9_C,
        "C_dagger": b.C_dagger,
        "closes_cr0006": True,
        "closes_c0007_all_ic": False,
        "certificate_ok": ok,
        "checks": {k: bool(v) for k, v in checks.items()},
        "not_a_clay_claim": True,
    }
    Path("reports").mkdir(exist_ok=True)
    Path("experiments/exploratory/sprint02_l0040").mkdir(parents=True, exist_ok=True)
    Path("experiments/exploratory/sprint02_l0040/summary.json").write_text(
        json.dumps(out, indent=2), encoding="utf-8"
    )
    Path("reports/SPRINT_02_L0040.md").write_text(
        f"""# Sprint — L-0040 Sym Shor + C-R-0006

| Item | Value |
|------|-------|
| C_Sym {{1..4}} | ≈{b.band_1234_C_sym:.6f} |
| C_Sym {{1..5}} | ≈{b.band_12345_C:.6f} ≤ C_† |
| C_Sym {{1..6}} | ≈{b.band_123456_C:.6f} > C_† |
| Greedy 9-shell | ≈{b.greedy9_C:.6f} |
| **C-R-0006** | shells {{1..5}} → Ω(T)≤M (proved) |
| C-0007 all-IC | still open |

FINITE Galerkin only. Not continuum. Not Clay.
""",
        encoding="utf-8",
    )
    print(json.dumps(out, indent=2, ensure_ascii=True))
    return out


if __name__ == "__main__":
    main()
