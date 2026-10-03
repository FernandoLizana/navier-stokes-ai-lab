"""Sprint: L-0028 embedding techo for C_† / C-0007 low slab."""

from __future__ import annotations

import json
from pathlib import Path

from ns_exploration.conjectures.l0028_embedding_techo import lemma_l0028, save_lemma_l0028
from ns_exploration.validation.l0028_certificate import (
    build_l0028_certificate,
    save_l0028_certificate,
    verify_l0028_certificate,
)


def main() -> dict:
    b = lemma_l0028(n=24, empirical=True, n_random=12, ascent_steps=8)
    save_lemma_l0028(b)
    cert = build_l0028_certificate(n=24)
    ok, checks = verify_l0028_certificate(cert.as_dict())
    save_l0028_certificate(cert)
    if not ok:
        raise RuntimeError(f"L-0028 certificate failed: {checks}")
    if b.embedding_can_meet_Cdagger:
        raise RuntimeError("unexpected: embedding meets C_dagger")

    c7_path = Path("conjectures/active/C-0007.json")
    c7 = json.loads(c7_path.read_text(encoding="utf-8"))
    c7["notes"] = (
        f"Open all-IC. C-R-0002 covers Ω0≥Ω★. L-0027: need stretch C≤C_†≈{b.C_dagger:.4f}. "
        f"L-0028: L∞/√M embeddings cannot meet C_† (need M≤{b.M_star:.2f}; r=1 already "
        f"M=6,C≈{b.C_shell_r1:.2f}). Next: triad-tensor stretch bound."
    )
    related = set(c7.get("related") or [])
    related.update(["C-R-0002", "L-0026", "L-0027", "L-0028"])
    c7["related"] = sorted(related)
    c7_path.write_text(json.dumps(c7, indent=2), encoding="utf-8")

    out = {
        "C_dagger": b.C_dagger,
        "M_star": b.M_star,
        "C_full_L0002": b.C_full_L0002,
        "C_shell_r1": b.C_shell_r1,
        "C_hard_R2": b.C_hard_R2,
        "embedding_can_meet_Cdagger": b.embedding_can_meet_Cdagger,
        "C_emp_max": b.C_emp_max,
        "closes_c0007": False,
        "certificate_ok": ok,
        "checks": {k: bool(v) for k, v in checks.items()},
        "not_a_clay_claim": True,
    }
    Path("reports").mkdir(exist_ok=True)
    Path("experiments/exploratory/sprint02_l0028").mkdir(parents=True, exist_ok=True)
    Path("experiments/exploratory/sprint02_l0028/summary.json").write_text(
        json.dumps(out, indent=2), encoding="utf-8"
    )
    Path("reports/SPRINT_02_L0028.md").write_text(
        f"""# Sprint — L-0028 embedding techo

| Item | Value |
|------|-------|
| C_† (L-0027) | ≈{b.C_dagger:.6f} |
| M_★ for L-0002 embedding | ≈{b.M_star:.4f} |
| Full-mask C (L-0002) | ≈{b.C_full_L0002:.4f} |
| r=1 shell C | ≈{b.C_shell_r1:.4f} |
| \|k\|²≤2 C | ≈{b.C_hard_R2:.4f} |
| Embedding meets C_†? | **No** |
| C_emp | ≈{b.C_emp_max:.6g} |
| C-0007 | still open |

Next: triad-aware cubic-tensor bound on stretch (not √M).

FINITE Galerkin only. Not continuum. Not Clay.
""",
        encoding="utf-8",
    )
    print(json.dumps(out, indent=2, ensure_ascii=True))
    return out


if __name__ == "__main__":
    main()
