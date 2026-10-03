"""Sprint: prove C-0004 (N≤24 envelope); open C-0005 (time-T below envelope)."""

from __future__ import annotations

import json
from pathlib import Path

from ns_exploration.conjectures.c0002_stokes_refuter import stokes_max_mode_ic
from ns_exploration.conjectures.c0004_n24_line import build_c0004_and_c0005
from ns_exploration.conjectures.l0019_stokes_majorant import lemma_l0019, save_lemma_l0019
from ns_exploration.diagnostics.metrics import compute_diagnostics
from ns_exploration.spectral.integrators import step_etd_rk2, step_rk4
from ns_exploration.validation.l0019_certificate import (
    build_l0019_certificate,
    save_l0019_certificate,
    verify_l0019_certificate,
)


def _gate_stokes_n24(nu: float = 0.1, dt: float = 1e-3, t_end: float = 0.02) -> dict:
    u0, meta = stokes_max_mode_ic(n=24, E0=0.5)
    u_etd, u_rk = u0.copy(), u0.copy()
    for _ in range(int(round(t_end / dt))):
        u_etd = step_etd_rk2(u_etd, nu, dt)
        u_rk = step_rk4(u_rk, nu, dt)
    om_e = compute_diagnostics(u_etd, nu).enstrophy
    om_r = compute_diagnostics(u_rk, nu).enstrophy
    rel = abs(om_e - om_r) / max(om_e, 1e-30)
    return {
        "K2": meta["K2"],
        "Omega_ETD": om_e,
        "Omega_RK4": om_r,
        "rel_ETD_RK4": rel,
        "gates_ok": rel < 1e-6,
    }


def main() -> dict:
    out = build_c0004_and_c0005()
    gates = _gate_stokes_n24()
    out["stokes_gates"] = gates
    assert gates["gates_ok"]
    assert gates["Omega_ETD"] < out["C0005_M"]
    assert gates["Omega_ETD"] > 40.0

    b19 = lemma_l0019()
    save_lemma_l0019(b19)
    cert19 = build_l0019_certificate()
    ok19, checks19 = verify_l0019_certificate(cert19.as_dict())
    save_l0019_certificate(cert19)
    out["L0019_Stokes_floor"] = b19.Stokes_floor
    out["L0019_cert_ok"] = ok19
    out["L0019_checks"] = {k: bool(v) for k, v in checks19.items()}
    assert ok19 and b19.below_c0005

    Path("reports").mkdir(exist_ok=True)
    Path("experiments/exploratory/sprint02_c0004_n24").mkdir(parents=True, exist_ok=True)
    Path("experiments/exploratory/sprint02_c0004_n24/summary.json").write_text(
        json.dumps(out, indent=2), encoding="utf-8"
    )
    print(json.dumps(out, indent=2, ensure_ascii=True))
    return out


if __name__ == "__main__":
    main()
