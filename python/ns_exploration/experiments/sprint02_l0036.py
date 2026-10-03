"""Sprint: L-0036 signed Shor stretch techo."""

from __future__ import annotations

import json
from pathlib import Path

from ns_exploration.conjectures.l0036_signed_shor_stretch import (
    lemma_l0036,
    save_lemma_l0036,
    signed_flattening_C_lower,
)
from ns_exploration.validation.l0036_certificate import (
    build_l0036_certificate,
    save_l0036_certificate,
    verify_l0036_certificate,
)


def main() -> dict:
    b = lemma_l0036(n=12, iters=16, include_n24_probe=False)
    d24 = signed_flattening_C_lower(n=24, iters=3)
    b.C_shor_lower_N24 = float(d24["C_shor_lower"])
    b.notes = (
        f"Signed Shor stretch: N=12 C_Shor≳{b.C_shor_lower:.4g}; "
        f"N=24 C_Shor≳{b.C_shor_lower_N24:.4g} (both ≫ C_†={b.C_dagger:.4g}). "
        f"Closes C-0007: False."
    )
    save_lemma_l0036(b)
    cert = build_l0036_certificate(n=12, iters=12)
    ok, checks = verify_l0036_certificate(cert.as_dict())
    save_l0036_certificate(cert)
    if not ok:
        raise RuntimeError(f"L-0036 certificate failed: {checks}")

    c7_path = Path("conjectures/active/C-0007.json")
    c7 = json.loads(c7_path.read_text(encoding="utf-8"))
    c7["notes"] = (
        f"Open all-IC. Subclasses C-R-0002/0003. L-0036: signed Shor C≳{b.C_shor_lower:.1f} "
        f"(N12) / ≳{b.C_shor_lower_N24:.0f} (N24) ≫ C_†. Need SOS/rank-constrained or "
        f"shell-diff blocks."
    )
    related = set(c7.get("related") or [])
    related.add("L-0036")
    c7["related"] = sorted(related)
    c7_path.write_text(json.dumps(c7, indent=2), encoding="utf-8")

    out = {
        "C_shor_lower_N12": b.C_shor_lower,
        "C_shor_lower_N24": b.C_shor_lower_N24,
        "C_dagger": b.C_dagger,
        "L_op_lower_N12": b.L_op_lower,
        "closes_c0007": False,
        "certificate_ok": ok,
        "checks": {k: bool(v) for k, v in checks.items()},
        "not_a_clay_claim": True,
    }
    Path("reports").mkdir(exist_ok=True)
    Path("experiments/exploratory/sprint02_l0036").mkdir(parents=True, exist_ok=True)
    Path("experiments/exploratory/sprint02_l0036/summary.json").write_text(
        json.dumps(out, indent=2), encoding="utf-8"
    )
    Path("reports/SPRINT_02_L0036.md").write_text(
        f"""# Sprint — L-0036 signed Shor stretch

| Item | Value |
|------|-------|
| C_Shor ≳ (N=12) | ≈{b.C_shor_lower:.4f} |
| C_Shor ≳ (N=24) | ≈{b.C_shor_lower_N24:.4f} |
| C_† | ≈{b.C_dagger:.4f} |
| C-0007 | still open |

FINITE Galerkin only. Not continuum. Not Clay.
""",
        encoding="utf-8",
    )
    print(json.dumps(out, indent=2, ensure_ascii=True))
    return out


if __name__ == "__main__":
    main()
