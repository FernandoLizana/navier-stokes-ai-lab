"""Sprint: L-0020 Duhamel H¹ → N_crit for C-0005 (Young too weak)."""

from __future__ import annotations

import json
from pathlib import Path

from ns_exploration.conjectures.l0020_duhamel_h1 import lemma_l0020, save_lemma_l0020
from ns_exploration.validation.l0020_certificate import (
    build_l0020_certificate,
    save_l0020_certificate,
    verify_l0020_certificate,
)


def main() -> dict:
    b = lemma_l0020()
    save_lemma_l0020(b)
    cert = build_l0020_certificate()
    ok, checks = verify_l0020_certificate(cert.as_dict())
    save_l0020_certificate(cert)

    # Update active C-0005 notes with N_crit target
    active_path = Path("conjectures/active/C-0005.json")
    if active_path.exists():
        c5 = json.loads(active_path.read_text(encoding="utf-8"))
        c5["L0020_Ncrit"] = b.Ncrit_c0005
        c5["L0020_Nstar_young"] = b.Nstar_young
        c5["L0020_Omega_young_H1"] = b.Omega_young_H1
        c5["notes"] = (
            f"Open: M={c5['proposed_bound_M']:.6f} = 1.5×Stokes_floor({c5['stokes_floor']:.6f}). "
            f"L-0020: closes if ‖N‖_ℓ₂ ≤ N_crit≈{b.Ncrit_c0005:.4f} on [0,T] "
            f"(empirical ‖N‖≲2; Young N_*≈{b.Nstar_young:.2f} too weak, Ω_H1≈{b.Omega_young_H1:.1f})."
        )
        active_path.write_text(json.dumps(c5, indent=2), encoding="utf-8")

    out = {
        "Stokes_floor": b.Stokes_floor,
        "I_sigma": b.I_sigma,
        "Nstar_young": b.Nstar_young,
        "Omega_young_H1": b.Omega_young_H1,
        "envelope_cap": b.envelope_cap,
        "c0005_M": b.c0005_M,
        "Ncrit_c0005": b.Ncrit_c0005,
        "young_closes_c0005": b.young_closes_c0005,
        "certificate_ok": ok,
        "checks": {k: bool(v) for k, v in checks.items()},
        "not_a_clay_claim": True,
        "next": "L-0023: prove stretch C <= C_star (~2.15); uniform N_* ruled out by L-0022.",
    }
    Path("reports").mkdir(exist_ok=True)
    Path("experiments/exploratory/sprint02_l0020").mkdir(parents=True, exist_ok=True)
    Path("experiments/exploratory/sprint02_l0020/summary.json").write_text(
        json.dumps(out, indent=2), encoding="utf-8"
    )
    print(json.dumps(out, indent=2, ensure_ascii=True))
    return out


if __name__ == "__main__":
    main()
