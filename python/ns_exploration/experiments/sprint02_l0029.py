"""Sprint: L-0029 Frobenius stretch embedding techo."""

from __future__ import annotations

import json
from pathlib import Path

from ns_exploration.conjectures.l0029_frobenius_stretch import (
    lemma_l0029,
    save_lemma_l0029,
)
from ns_exploration.validation.l0029_certificate import (
    build_l0029_certificate,
    save_l0029_certificate,
    verify_l0029_certificate,
)


def main() -> dict:
    b = lemma_l0029(n=24)
    save_lemma_l0029(b)
    cert = build_l0029_certificate(n=24)
    ok, checks = verify_l0029_certificate(cert.as_dict())
    save_l0029_certificate(cert)
    if not ok:
        raise RuntimeError(f"L-0029 certificate failed: {checks}")
    if b.frobenius_closes_c0007:
        raise RuntimeError("unexpected closure")

    c7_path = Path("conjectures/active/C-0007.json")
    c7 = json.loads(c7_path.read_text(encoding="utf-8"))
    c7["notes"] = (
        f"Open all-IC. C-R-0002: high slab. L-0027: need C≤C_†≈{b.C_dagger:.4f}. "
        f"L-0028/29: embeddings fail (L-0029 Frobenius C_F≈{b.C_F_full:.2f} on full mask; "
        f"r=1 meets C_† but λ_max=1≪λ★≈{b.lambda_star:.2f}). Next: triad-tensor stretch."
    )
    related = set(c7.get("related") or [])
    related.update(["C-R-0002", "L-0026", "L-0027", "L-0028", "L-0029"])
    c7["related"] = sorted(related)
    c7_path.write_text(json.dumps(c7, indent=2), encoding="utf-8")

    out = {
        "C_dagger": b.C_dagger,
        "C_F_full": b.C_F_full,
        "C_L0002_full": b.C_L0002_full,
        "M_star_F": b.M_star_F,
        "C_F_shell_r1": b.C_F_shell_r1,
        "shell_r1_meets_Cdagger": b.shell_r1_meets_Cdagger,
        "shell_r1_reaches_low_slab_edge": b.shell_r1_reaches_low_slab_edge,
        "frobenius_closes_c0007": False,
        "certificate_ok": ok,
        "checks": {k: bool(v) for k, v in checks.items()},
        "not_a_clay_claim": True,
    }
    Path("reports").mkdir(exist_ok=True)
    Path("experiments/exploratory/sprint02_l0029").mkdir(parents=True, exist_ok=True)
    Path("experiments/exploratory/sprint02_l0029/summary.json").write_text(
        json.dumps(out, indent=2), encoding="utf-8"
    )
    Path("reports/SPRINT_02_L0029.md").write_text(
        f"""# Sprint — L-0029 Frobenius stretch embedding

| Item | Value |
|------|-------|
| C_F (full mask) | ≈{b.C_F_full:.4f} |
| C L-0002 (full) | ≈{b.C_L0002_full:.4f} |
| Improvement | √3 ≈ {b.improvement_factor:.4f} |
| C_† | ≈{b.C_dagger:.4f} |
| M_†^F | ≈{b.M_star_F:.4f} |
| r=1 C_F | ≈{b.C_F_shell_r1:.4f} (meets C_†, cannot reach λ★) |
| C-0007 | still open |

FINITE Galerkin only. Not continuum. Not Clay.
""",
        encoding="utf-8",
    )
    print(json.dumps(out, indent=2, ensure_ascii=True))
    return out


if __name__ == "__main__":
    main()
