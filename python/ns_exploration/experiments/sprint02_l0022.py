"""Sprint: L-0022 bootstrap arithmetic + adversarial N_* reality check."""

from __future__ import annotations

import json
from pathlib import Path

from ns_exploration.conjectures.l0022_gamma_bootstrap import lemma_l0022, save_lemma_l0022
from ns_exploration.validation.l0022_certificate import (
    build_l0022_certificate,
    save_l0022_certificate,
    verify_l0022_certificate,
)


def main() -> dict:
    b = lemma_l0022(empirical=True, n_random=10, ascent_steps=18)
    save_lemma_l0022(b)
    cert = build_l0022_certificate(empirical=False)
    cert.N_emp_max = b.N_emp_max
    cert.gamma_emp_max = b.gamma_emp_max
    cert.Omega_T_emp_max = b.Omega_T_emp_max
    cert.uniform_hypothesis_holds_empirically = b.uniform_hypothesis_holds_empirically
    cert.notes = b.notes
    ok, checks = verify_l0022_certificate(cert.as_dict())
    save_l0022_certificate(cert)

    active = Path("conjectures/active/C-0005.json")
    if active.exists():
        c5 = json.loads(active.read_text(encoding="utf-8"))
        c5["notes"] = (
            f"Open all-IC: M={c5['proposed_bound_M']:.6f}. "
            f"L-0020 N_crit≈{b.Ncrit_c0005:.4f}; L-0022 γ_crit≈{b.gamma_crit:.4f}. "
            f"Adversarial ‖N‖_emp≈{b.N_emp_max:.2f}≫N_crit and γ_emp≈{b.gamma_emp_max:.2f}≫γ_crit "
            f"(uniform Duhamel bootstrap cannot close). "
            f"Ω(T)_emp≲{b.Omega_T_emp_max:.2f}≪M so C-0005 not refuted; need a sharper majorant."
        )
        c5["L0022_gamma_crit"] = b.gamma_crit
        c5["L0022_gamma_emp"] = b.gamma_emp_max
        c5["L0022_N_emp"] = b.N_emp_max
        c5["L0022_Omega_T_emp"] = b.Omega_T_emp_max
        active.write_text(json.dumps(c5, indent=2), encoding="utf-8")

    out = {
        "gamma_crit": b.gamma_crit,
        "Ncrit": b.Ncrit_c0005,
        "N_emp_max": b.N_emp_max,
        "gamma_emp_max": b.gamma_emp_max,
        "Omega_T_emp_max": b.Omega_T_emp_max,
        "uniform_hypothesis_holds_empirically": b.uniform_hypothesis_holds_empirically,
        "proved_closes_c0005": False,
        "certificate_ok": ok,
        "checks": {k: bool(v) for k, v in checks.items()},
        "not_a_clay_claim": True,
    }
    Path("reports").mkdir(exist_ok=True)
    Path("experiments/exploratory/sprint02_l0022").mkdir(parents=True, exist_ok=True)
    Path("experiments/exploratory/sprint02_l0022/summary.json").write_text(
        json.dumps(out, indent=2), encoding="utf-8"
    )
    print(json.dumps(out, indent=2, ensure_ascii=True))
    return out


if __name__ == "__main__":
    main()
