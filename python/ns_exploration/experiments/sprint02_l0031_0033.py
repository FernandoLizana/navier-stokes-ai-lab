"""Sprint: L-0031/32/33 techos + C-R-0003 (α=0 mono-radial)."""

from __future__ import annotations

import json
from pathlib import Path

from ns_exploration.conjectures.l0031_0033_techos import (
    lemma_l0031_33,
    save_lemma_l0031_33,
)
from ns_exploration.validation.l0031_0033_certificate import (
    build_l0031_33_certificate,
    save_l0031_33_certificate,
    verify_l0031_33_certificate,
)


def main() -> dict:
    b = lemma_l0031_33(n=24, n_grid=21)
    save_lemma_l0031_33(b)
    cert = build_l0031_33_certificate(n=24)
    ok, checks = verify_l0031_33_certificate(cert.as_dict())
    save_l0031_33_certificate(cert)
    if not ok:
        raise RuntimeError(f"certificate failed: {checks}")

    cr = {
        "id": "C-R-0003",
        "route": "B",
        "status": "proved_restricted",
        "evidence_level": "N7",
        "lemma": "L-0020+L-0021(α=0)+L-0031/33",
        "domain": (
            "Pseudospectral T^3, N<=24, 2/3 dealias, mono-radial IC on a shell "
            "with quartic alpha_r=0 (no self-triads). NOT continuum."
        ),
        "nu": 0.1,
        "energy": 0.5,
        "t_end": 0.02,
        "n_max": 24,
        "proposed_bound_M": b.c0007_M,
        "zero_alpha_shells": b.zero_alpha_shells,
        "proved_bound_Omega_T": b.zero_alpha_duhamel_OmT,
        "parent_open": "C-0007",
        "statement": (
            f"For mono-radial ICs on shells r∈{b.zero_alpha_shells} (α_r=0), "
            f"dealiased Galerkin NS (N≤24, ν=0.1, E≤0.5) satisfies "
            f"Ω(0.02)≤{b.zero_alpha_duhamel_OmT:.8f}≤{b.c0007_M:.8f} "
            f"(Stokes–Duhamel with N_*=0). FINITE-DIMENSIONAL only."
        ),
        "clay_implication": (
            "None. Mono-radial α=0 subclass only; not all-IC C-0007; not Clay."
        ),
        "notes": (
            "Shells with vanishing L-0021 quartic α (no dealias self-triads). "
            "Duhamel collapses to Stokes floor."
        ),
    }
    Path("conjectures/proved_restricted/C-R-0003.json").write_text(
        json.dumps(cr, indent=2), encoding="utf-8"
    )

    c7_path = Path("conjectures/active/C-0007.json")
    c7 = json.loads(c7_path.read_text(encoding="utf-8"))
    c7["notes"] = (
        f"Open all-IC. Closed subclasses: C-R-0002 (Ω0≥Ω★), C-R-0003 (mono-radial "
        f"α=0 shells {b.zero_alpha_shells}). Techos L-0026…L-0033 block embeddings, "
        f"raw R_★, compatible-S (Ω(T)≈{b.compatible_worst_OmT:.3f}), full S_max, "
        f"absolute cubic F (C≲{b.C_from_cubic_F:.0f}). Need signed/SOS triad stretch."
    )
    related = set(c7.get("related") or [])
    related.update(
        [
            "C-R-0002",
            "C-R-0003",
            "L-0026",
            "L-0027",
            "L-0028",
            "L-0029",
            "L-0030",
            "L-0031",
            "L-0032",
            "L-0033",
        ]
    )
    c7["related"] = sorted(related)
    c7_path.write_text(json.dumps(c7, indent=2), encoding="utf-8")

    out = {
        "compatible_worst_OmT": b.compatible_worst_OmT,
        "L0024_majorant": b.L0024_majorant,
        "S_max_full": b.S_max_full,
        "N_from_Smax": b.N_from_Smax,
        "N_young": b.N_young,
        "C_from_cubic_F": b.C_from_cubic_F,
        "C_dagger": b.C_dagger,
        "zero_alpha_shells": b.zero_alpha_shells,
        "C_R_0003_OmT": b.zero_alpha_duhamel_OmT,
        "closes_c0007": False,
        "certificate_ok": ok,
        "checks": {k: bool(v) for k, v in checks.items()},
        "not_a_clay_claim": True,
    }
    Path("reports").mkdir(exist_ok=True)
    Path("experiments/exploratory/sprint02_l0031_0033").mkdir(parents=True, exist_ok=True)
    Path("experiments/exploratory/sprint02_l0031_0033/summary.json").write_text(
        json.dumps(out, indent=2), encoding="utf-8"
    )
    Path("reports/SPRINT_02_L0031_0033.md").write_text(
        f"""# Sprint — L-0031/32/33 techos + C-R-0003

| Item | Result |
|------|--------|
| L-0031 compatible-S Ω(T) | ≈{b.compatible_worst_OmT:.6f} (worse than L-0024) |
| L-0032 full S_max N | ≈{b.N_from_Smax:.2f} ≫ Young ≈{b.N_young:.2f} |
| L-0033 absolute cubic C | ≲{b.C_from_cubic_F:.1f} ≫ C_†≈{b.C_dagger:.2f} |
| **C-R-0003** α=0 shells | {b.zero_alpha_shells} → Ω≤{b.zero_alpha_duhamel_OmT:.6f}≤M |
| C-0007 all-IC | **still open** |

FINITE Galerkin only. Not continuum. Not Clay.
""",
        encoding="utf-8",
    )
    print(json.dumps(out, indent=2, ensure_ascii=True))
    return out


if __name__ == "__main__":
    main()
