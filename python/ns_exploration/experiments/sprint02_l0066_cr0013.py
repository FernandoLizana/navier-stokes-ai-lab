"""Sprint: close C-R-0013 (L-0066) + refresh closure map and Lean bridge."""

from __future__ import annotations

import json
from pathlib import Path

from ns_exploration.conjectures.l0056_closure_map import (
    lemma_l0056,
    save_lemma_l0056,
    write_closure_map_report,
)
from ns_exploration.conjectures.l0059_lean_sos_cert_bridge import lemma_l0059, save_lemma_l0059
from ns_exploration.conjectures.l0066_cr0013_onepol_greedy import (
    lemma_l0066,
    save_cr0013,
    save_lemma_l0066,
)
from ns_exploration.validation.l0066_certificate import (
    build_l0066_certificate,
    save_l0066_certificate,
    verify_l0066_certificate,
)


def _patch_scale_results() -> None:
    path = Path("tools/sos_julia/data/scale_results.json")
    if not path.is_file():
        return
    data = json.loads(path.read_text(encoding="utf-8"))
    entry = {
        "support": [1, 2, 3, 4, 5, 6, 8, 10, 12, 16, 19],
        "branch": "second",
        "D": 154,
        "C_shor_sym": 9.189036845789548,
        "C_fullsym": 9.060826005460791,
        "C_ub": 4.01493826157398,
        "C_ub_hi": 4.014938261573984,
        "solver": "COSMO",
        "beats_shor": True,
        "closes_vs_Cdagger": True,
        "solve_s_approx": 217,
        "cert": "CERT-L0065-sos-onepol-greedy-second-N24",
        "restricted": "C-R-0013",
    }
    rows = data.setdefault("onepol_greedy_second_n5", entry)
    if isinstance(rows, dict):
        data["onepol_greedy_second_n5"] = entry
    path.write_text(json.dumps(data, indent=2), encoding="utf-8")


def main() -> dict:
    b = lemma_l0066(n=24)
    save_lemma_l0066(b)
    save_cr0013(b)
    cert = build_l0066_certificate(n=24)
    ok, checks = verify_l0066_certificate(cert.as_dict())
    save_l0066_certificate(cert)
    if not ok:
        raise RuntimeError(f"L-0066 certificate failed: {checks}")

    l56 = lemma_l0056(n=24)
    save_lemma_l0056(l56)
    write_closure_map_report(l56)

    l59 = lemma_l0059(n=24, run_lean=False)
    save_lemma_l0059(l59)

    _patch_scale_results()

    c7_path = Path("conjectures/active/C-0007.json")
    c7 = json.loads(c7_path.read_text(encoding="utf-8"))
    related = set(c7.get("related") or [])
    related.update(["L-0066", "C-R-0013", "L-0065", "L-0064"])
    c7["related"] = sorted(related)
    c7["notes"] = (
        f"Open all-IC. Closed: C-R-0002..0013 ({l56.n_proved_restricted} restricted). "
        f"C-R-0013 one-pol 11 shells D={b.support_cr0013_D} "
        f"C_Shor≈{b.C_shor_sym_cr0013:.4f} SOS≈{b.C_ub_sos_hi:.4f}."
    )
    c7_path.write_text(json.dumps(c7, indent=2), encoding="utf-8")

    out = {
        "C-R-0013": True,
        "C_shor_sym": b.C_shor_sym_cr0013,
        "C_ub_sos_hi": b.C_ub_sos_hi,
        "D": b.support_cr0013_D,
        "n_proved_restricted": l56.n_proved_restricted,
        "certificate_ok": ok,
        "checks": {k: bool(v) for k, v in checks.items()},
        "not_a_clay_claim": True,
    }
    Path("experiments/exploratory/sprint02_l0066").mkdir(parents=True, exist_ok=True)
    Path("experiments/exploratory/sprint02_l0066/summary.json").write_text(
        json.dumps(out, indent=2), encoding="utf-8"
    )
    print(json.dumps(out, indent=2, ensure_ascii=True))
    return out


if __name__ == "__main__":
    main()
