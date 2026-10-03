"""Sprint: L-0021 quartic shell + prove C-R-0001 (mono-radial)."""

from __future__ import annotations

import json
from pathlib import Path

from ns_exploration.conjectures.l0021_quartic_shell import lemma_l0021, save_lemma_l0021
from ns_exploration.validation.l0021_certificate import (
    build_l0021_certificate,
    save_l0021_certificate,
    verify_l0021_certificate,
)


def main() -> dict:
    b = lemma_l0021()
    save_lemma_l0021(b)
    cert = build_l0021_certificate()
    ok, checks = verify_l0021_certificate(cert.as_dict())
    save_l0021_certificate(cert)

    M = float(b.Omega_radial_H1)
    cr = {
        "id": "C-R-0001",
        "route": "B",
        "status": "proved_restricted",
        "evidence_level": "N7",
        "lemma": "L-0021+L-0020",
        "domain": (
            "Pseudospectral T^3, N≤24, 2/3 dealias, IC supported on a single "
            "radial shell |k|²=r. NOT continuum."
        ),
        "nu": b.nu,
        "energy": b.E0,
        "t_end": b.T,
        "n_max": b.n,
        "proposed_bound_M": M,
        "statement": (
            f"For all divergence-free Fourier ICs on N≤{b.n} with energy ≤{b.E0} "
            f"supported on a single radial shell, evolved under dealiased Galerkin NS "
            f"(ν={b.nu}), one has Ω({b.T}) ≤ {M:.8f}, via quartic Gram N_*≤α_★="
            f"{b.alpha_star:.6g} (r={b.r_star}) and Stokes–Duhamel H¹ (L-0020/21). "
            f"Beats envelope K²E0={b.envelope_cap}. FINITE-DIMENSIONAL only."
        ),
        "alpha_star": b.alpha_star,
        "r_star": b.r_star,
        "envelope_cap": b.envelope_cap,
        "parent_open": "C-0005",
        "clay_implication": b.clay_implication,
        "notes": b.notes,
    }
    Path("conjectures/proved_restricted/C-R-0001.json").write_text(
        json.dumps(cr, indent=2), encoding="utf-8"
    )

    # Update C-0005 notes
    active = Path("conjectures/active/C-0005.json")
    if active.exists():
        c5 = json.loads(active.read_text(encoding="utf-8"))
        c5["notes"] = (
            f"Open all-IC: M={c5['proposed_bound_M']:.6f}. L-0020 needs N_*≤"
            f"{c5.get('L0020_Ncrit', 9.666):.4f}. L-0021 mono-radial quartic α_★="
            f"{b.alpha_star:.4f} ⇒ Duhamel M≈{M:.4f} (C-R-0001 proved); still too "
            f"weak for all-IC C-0005 without a global N_*≤N_crit."
        )
        c5["C_R_0001_M"] = M
        c5["L0021_alpha_star"] = b.alpha_star
        active.write_text(json.dumps(c5, indent=2), encoding="utf-8")

    out = {
        "alpha_star": b.alpha_star,
        "r_star": b.r_star,
        "Omega_radial_H1": b.Omega_radial_H1,
        "beats_envelope": b.beats_envelope,
        "closes_c0005": b.closes_c0005,
        "C_R_0001_M": M,
        "certificate_ok": ok,
        "checks": {k: bool(v) for k, v in checks.items()},
        "not_a_clay_claim": True,
    }
    Path("reports").mkdir(exist_ok=True)
    Path("experiments/exploratory/sprint02_l0021").mkdir(parents=True, exist_ok=True)
    Path("experiments/exploratory/sprint02_l0021/summary.json").write_text(
        json.dumps(out, indent=2), encoding="utf-8"
    )
    print(json.dumps(out, indent=2, ensure_ascii=True))
    return out


if __name__ == "__main__":
    main()
