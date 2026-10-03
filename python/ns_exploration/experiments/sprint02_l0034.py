"""Sprint: L-0034 shell-factorized ‖N‖ bound."""

from __future__ import annotations

import json
from pathlib import Path

from ns_exploration.conjectures.l0034_shell_N_bound import lemma_l0034, save_lemma_l0034
from ns_exploration.validation.l0034_certificate import (
    build_l0034_certificate,
    save_l0034_certificate,
    verify_l0034_certificate,
)


def main() -> dict:
    b = lemma_l0034(n=24, n_grid=11)
    save_lemma_l0034(b)
    cert = build_l0034_certificate(n=24)
    ok, checks = verify_l0034_certificate(cert.as_dict())
    save_l0034_certificate(cert)
    if not ok:
        raise RuntimeError(f"L-0034 certificate failed: {checks}")

    c7_path = Path("conjectures/active/C-0007.json")
    c7 = json.loads(c7_path.read_text(encoding="utf-8"))
    c7["notes"] = (
        f"Open all-IC. Subclasses C-R-0002/0003. L-0034: shell N≤2.45 at Ω=E "
        f"(≪ Young {b.N_young:.1f}) but C_eff≳{b.C_eff_at_equipartition:.0f}≫C_† "
        f"and low-slab ODE≈{b.low_slab_ode_worst:.1f}>M. Need signed triad SOS."
    )
    related = set(c7.get("related") or [])
    related.add("L-0034")
    c7["related"] = sorted(related)
    c7_path.write_text(json.dumps(c7, indent=2), encoding="utf-8")

    out = {
        "N_at_equipartition": b.N_at_equipartition,
        "N_young": b.N_young,
        "C_eff_at_equipartition": b.C_eff_at_equipartition,
        "C_eff_at_Omega_star": b.C_eff_at_Omega_star,
        "low_slab_ode_worst": b.low_slab_ode_worst,
        "C_dagger": b.C_dagger,
        "closes_c0007": False,
        "certificate_ok": ok,
        "checks": {k: bool(v) for k, v in checks.items()},
        "not_a_clay_claim": True,
    }
    Path("reports").mkdir(exist_ok=True)
    Path("experiments/exploratory/sprint02_l0034").mkdir(parents=True, exist_ok=True)
    Path("experiments/exploratory/sprint02_l0034/summary.json").write_text(
        json.dumps(out, indent=2), encoding="utf-8"
    )
    Path("reports/SPRINT_02_L0034.md").write_text(
        f"""# Sprint — L-0034 shell-factorized ‖N‖

| Item | Value |
|------|-------|
| N @ Ω=E | ≈{b.N_at_equipartition:.4f} (Young ≈{b.N_young:.4f}) |
| C_eff @ Ω=E | ≈{b.C_eff_at_equipartition:.2f} |
| C_† | ≈{b.C_dagger:.4f} |
| Low-slab ODE | ≈{b.low_slab_ode_worst:.4f} |
| C-0007 | still open |

FINITE Galerkin only. Not continuum. Not Clay.
""",
        encoding="utf-8",
    )
    print(json.dumps(out, indent=2, ensure_ascii=True))
    return out


if __name__ == "__main__":
    main()
