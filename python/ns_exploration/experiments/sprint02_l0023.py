"""Sprint: L-0023 cubic-dissipation C_★ toward C-0005."""

from __future__ import annotations

import json
from pathlib import Path

from ns_exploration.conjectures.l0023_cubic_dissipation import lemma_l0023, save_lemma_l0023
from ns_exploration.validation.l0023_certificate import (
    build_l0023_certificate,
    save_l0023_certificate,
    verify_l0023_certificate,
)


def main() -> dict:
    b = lemma_l0023(empirical=True, n_random=20, ascent_steps=18)
    save_lemma_l0023(b)
    cert = build_l0023_certificate(empirical=False)
    cert.C_emp_max = b.C_emp_max
    cert.emp_closes = b.emp_closes
    cert.notes = b.notes
    ok, checks = verify_l0023_certificate(cert.as_dict())
    save_l0023_certificate(cert)

    active = Path("conjectures/active/C-0005.json")
    if active.exists():
        c5 = json.loads(active.read_text(encoding="utf-8"))
        c5["notes"] = (
            f"Open all-IC: M={c5['proposed_bound_M']:.6f}. "
            f"L-0022: uniform N_*/γ ruled out. "
            f"L-0023: closes if stretch C≤C_★≈{b.C_star:.4f} in "
            f"|⟨ω,curl N⟩|≤C Ω^{{3/2}} (cubic dissipation ODE); "
            f"C_emp≈{b.C_emp_max:.4f}≪C_★; proved C≤C_★ still open."
        )
        c5["L0023_C_star"] = b.C_star
        c5["L0023_C_emp"] = b.C_emp_max
        active.write_text(json.dumps(c5, indent=2), encoding="utf-8")

    out = {
        "C_star": b.C_star,
        "Omega_at_Cstar": b.Omega_at_Cstar,
        "C_emp_max": b.C_emp_max,
        "emp_closes": b.emp_closes,
        "proved_closes_c0005": False,
        "certificate_ok": ok,
        "checks": {k: bool(v) for k, v in checks.items()},
        "not_a_clay_claim": True,
    }
    Path("reports").mkdir(exist_ok=True)
    Path("experiments/exploratory/sprint02_l0023").mkdir(parents=True, exist_ok=True)
    Path("experiments/exploratory/sprint02_l0023/summary.json").write_text(
        json.dumps(out, indent=2), encoding="utf-8"
    )
    print(json.dumps(out, indent=2, ensure_ascii=True))
    return out


if __name__ == "__main__":
    main()
