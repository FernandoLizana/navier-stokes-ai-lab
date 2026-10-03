"""Sprint: L-0035 shell two-point ‖∇u‖_∞ stretch techo."""

from __future__ import annotations

import json
from pathlib import Path

from ns_exploration.conjectures.l0035_ginf_stretch import lemma_l0035, save_lemma_l0035
from ns_exploration.validation.l0035_certificate import (
    build_l0035_certificate,
    save_l0035_certificate,
    verify_l0035_certificate,
)


def main() -> dict:
    b = lemma_l0035(n=24, n_grid=9)
    save_lemma_l0035(b)
    cert = build_l0035_certificate(n=24)
    ok, checks = verify_l0035_certificate(cert.as_dict())
    save_l0035_certificate(cert)
    if not ok:
        raise RuntimeError(f"L-0035 certificate failed: {checks}")

    c7_path = Path("conjectures/active/C-0007.json")
    c7 = json.loads(c7_path.read_text(encoding="utf-8"))
    c7["notes"] = (
        f"Open all-IC. Subclasses C-R-0002/0003. L-0035: C_ginf(E)≈{b.C_at_equipartition:.2f}≤C_† "
        f"but C_max(low)≈{b.C_max_low_slab:.1f} and min-prod ODE≈{b.minprod_ode_worst:.1f}>M. "
        f"Need signed triad SOS."
    )
    related = set(c7.get("related") or [])
    related.add("L-0035")
    c7["related"] = sorted(related)
    c7_path.write_text(json.dumps(c7, indent=2), encoding="utf-8")

    out = {
        "C_at_equipartition": b.C_at_equipartition,
        "C_max_low_slab": b.C_max_low_slab,
        "Omega_c": b.Omega_c,
        "C_dagger": b.C_dagger,
        "ginf_ode_worst": b.ginf_ode_worst,
        "minprod_ode_worst": b.minprod_ode_worst,
        "closes_c0007": False,
        "certificate_ok": ok,
        "checks": {k: bool(v) for k, v in checks.items()},
        "not_a_clay_claim": True,
    }
    Path("reports").mkdir(exist_ok=True)
    Path("experiments/exploratory/sprint02_l0035").mkdir(parents=True, exist_ok=True)
    Path("experiments/exploratory/sprint02_l0035/summary.json").write_text(
        json.dumps(out, indent=2), encoding="utf-8"
    )
    Path("reports/SPRINT_02_L0035.md").write_text(
        f"""# Sprint — L-0035 shell two-point ‖∇u‖_∞

| Item | Value |
|------|-------|
| C @ Ω=E | ≈{b.C_at_equipartition:.4f} ≤ C_†≈{b.C_dagger:.4f} |
| C_max (low slab) | ≈{b.C_max_low_slab:.4f} |
| Ω_c (C≤C_†) | ≈{b.Omega_c:.4f} |
| ginf ODE | ≈{b.ginf_ode_worst:.4f} |
| min-prod ODE | ≈{b.minprod_ode_worst:.4f} |
| C-0007 | still open |

FINITE Galerkin only. Not continuum. Not Clay.
""",
        encoding="utf-8",
    )
    print(json.dumps(out, indent=2, ensure_ascii=True))
    return out


if __name__ == "__main__":
    main()
