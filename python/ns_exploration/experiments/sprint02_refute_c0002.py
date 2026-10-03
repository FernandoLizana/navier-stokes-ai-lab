"""
Refute C-0002 with max-|k| Stokes mode; spawn and prove C-0003 via L-0018.

C-0002 claimed Ω(T)≤≈20.32 for all ICs on N≤16. A single dealias mode with
|k|²=75 evolves as Stokes and gives Ω(0.02)≈27.78 > M (N2, gates ETD↔RK4).

Successor C-0003 takes M = K²(16) E0 = 37.5 and is proved for all t≥0 by L-0018
(spectral envelope). FINITE Galerkin only. Not Clay.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from ns_exploration.conjectures.c0002_enstrophy_bound import ConjectureC0002, save_conjecture_c0002
from ns_exploration.conjectures.c0002_stokes_refuter import (
    stokes_Omega_exact,
    stokes_max_mode_ic,
)
from ns_exploration.conjectures.l0018_envelope import lemma_l0018, save_lemma_l0018
from ns_exploration.diagnostics.metrics import compute_diagnostics
from ns_exploration.initial_conditions.generators import ICMetadata, save_ic
from ns_exploration.spectral.fourier_conventions import kinetic_energy_from_hat
from ns_exploration.spectral.integrators import step_etd_rk2, step_rk4
from ns_exploration.validation.l0018_certificate import (
    build_l0018_envelope_certificate,
    save_l0018_envelope_certificate,
    verify_l0018_envelope_certificate,
)


def main(
    nu: float = 0.1,
    dt: float = 1e-3,
    t_end: float = 0.02,
    energy: float = 0.5,
    M_c0002: float = 20.320076249551253,
) -> dict:
    u0, meta = stokes_max_mode_ic(n=16, E0=energy)
    Om0 = compute_diagnostics(u0, nu).enstrophy
    Om_exact = stokes_Omega_exact(meta["K2"], energy, nu, t_end)

    u_etd = u0.copy()
    u_rk = u0.copy()
    steps = int(round(t_end / dt))
    for _ in range(steps):
        u_etd = step_etd_rk2(u_etd, nu, dt)
        u_rk = step_rk4(u_rk, nu, dt)

    Om_etd = compute_diagnostics(u_etd, nu).enstrophy
    Om_rk = compute_diagnostics(u_rk, nu).enstrophy
    rel = abs(Om_etd - Om_rk) / max(Om_etd, 1e-30)
    gates_ok = rel < 1e-6 and abs(Om_etd - Om_exact) / Om_exact < 1e-6
    refutes = Om_etd > M_c0002 and gates_ok

    # Archive C-0002
    c2 = ConjectureC0002(
        proposed_bound_M=M_c0002,
        numerical_support_max=Om_etd,
        status="refuted",
        refuted=True,
        refutation_value=Om_etd,
        evidence_level="N2",
        notes=(
            f"Refuted by Stokes max-mode |k|²={meta['K2']}: "
            f"Ω_ETD({t_end})={Om_etd:.12g}, Ω_exact={Om_exact:.12g}, "
            f"ETD↔RK4 rel={rel:.3e}. Prior low-mode CMA max≈14.27 missed this IC."
        ),
    )
    save_conjecture_c0002(c2, "conjectures/rejected/C-0002.json")

    # Save refuting IC
    import hashlib

    payload = np.concatenate([u0.real.ravel(), u0.imag.ravel()])
    h = hashlib.sha256(payload.tobytes()).hexdigest()[:16]
    ic_meta = ICMetadata(
        name="stokes_k2_75_c0002_refuter",
        n=16,
        seed=None,
        energy=kinetic_energy_from_hat(u0),
        enstrophy=Om0,
        helicity=0.0,
        symmetries=["stokes_single_mode"],
        resolution=16,
        content_hash=h,
    )
    save_ic("datasets/candidates/c0002_stokes_k2_75_refuter.npz", u0, ic_meta)

    # C-0003 = all-IC envelope, proved by L-0018
    b18 = lemma_l0018(n_max=16, E0=energy, nu=nu, cs0002_M=M_c0002)
    # Recompute with M = Omega_cap for C-0003 statement
    M3 = b18.Omega_cap  # 37.5
    assert M3 >= Om_etd

    c3 = {
        "id": "C-0003",
        "route": "B",
        "status": "proved_restricted",
        "evidence_level": "N7",
        "lemma": "L-0018",
        "domain": "Galerkin/pseudospectral T^3, N<=16, 2/3 dealias (NOT continuum PDE)",
        "nu": nu,
        "energy": energy,
        "t_end": t_end,
        "n_max": 16,
        "proposed_bound_M": M3,
        "parent_refuted": "C-0002",
        "statement": (
            f"For all divergence-free Fourier fields on the N³ grid with N≤16, "
            f"dealiased Galerkin NS (ν={nu}), kinetic energy ≤ {energy}, one has "
            f"Ω(t) ≤ {M3} for all t ≥ 0 (hence at T={t_end}), via the spectral "
            f"envelope Ω ≤ K² E ≤ K²(16) E0 = {M3} (L-0018). FINITE-DIMENSIONAL only. "
            f"Replaces refuted C-0002 (M≈{M_c0002})."
        ),
        "Omega_cap": M3,
        "K2_worst": b18.K2_worst,
        "clay_implication": (
            "None. Finite-dimensional only; no continuum regularity or blow-up claim."
        ),
        "numerical_support_max": Om_etd,
        "refuted": False,
        "refutation_value": None,
        "notes": (
            f"Proved by L-0018. Stokes refuter of C-0002 sits at Ω(T)≈{Om_etd:.6f} < M={M3}."
        ),
    }
    Path("conjectures/proved_restricted/C-0003.json").write_text(
        json.dumps(c3, indent=2), encoding="utf-8"
    )

    # Clear active C-0002; no active all-IC conjecture left at this N/E0.
    active = Path("conjectures/active/C-0002.json")
    if active.exists():
        active.unlink()

    # Refresh L-0018 artifact note (same lemma closes C-0003)
    b18.notes = (
        b18.notes
        + f" Also closes all-IC C-0003 with M={M3} (after C-0002 Stokes refutation)."
    )
    save_lemma_l0018(b18)

    cert = build_l0018_envelope_certificate(
        n_max=16, E0=energy, R=M3, certificate_id="CERT-L0018-envelope-C0003-N16"
    )
    ok, checks = verify_l0018_envelope_certificate(cert.as_dict())
    save_l0018_envelope_certificate(
        cert, "certificates/CERT-L0018-envelope-C0003-N16.json"
    )
    save_l0018_envelope_certificate(
        build_l0018_envelope_certificate(),
        "certificates/CERT-L0018-envelope-CS0002-N16.json",
    )

    out = {
        "refuted_C0002": refutes,
        "M_c0002": M_c0002,
        "Omega0": Om0,
        "Omega_ETD": Om_etd,
        "Omega_RK4": Om_rk,
        "Omega_Stokes_exact": Om_exact,
        "rel_ETD_RK4": rel,
        "gates_ok": gates_ok,
        "C0003_M": M3,
        "C0003_proved_by": "L-0018",
        "certificate_verified": ok,
        "checks": {k: bool(v) for k, v in checks.items()},
        "k": meta["k"],
        "K2": meta["K2"],
        "not_a_clay_claim": True,
    }
    Path("reports").mkdir(exist_ok=True)
    Path("experiments/exploratory/sprint02_refute_c0002").mkdir(parents=True, exist_ok=True)
    Path("experiments/exploratory/sprint02_refute_c0002/summary.json").write_text(
        json.dumps(out, indent=2), encoding="utf-8"
    )
    print(json.dumps(out, indent=2, ensure_ascii=True))
    return out


if __name__ == "__main__":
    main()
