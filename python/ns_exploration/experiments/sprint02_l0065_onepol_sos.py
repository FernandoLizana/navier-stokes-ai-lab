"""Sprint: L-0065 one-pol greedy-second COSMO SOS + N5 certificate."""

from __future__ import annotations

import json
from pathlib import Path

from ns_exploration.conjectures.l0065_onepol_greedy_second_sos import lemma_l0065, save_lemma_l0065
from ns_exploration.validation.l0065_certificate import (
    build_l0065_certificate,
    save_l0065_certificate,
    verify_l0065_certificate,
)


def main() -> dict:
    b = lemma_l0065(n=24)
    save_lemma_l0065(b)
    cert = build_l0065_certificate(n=24)
    ok, checks = verify_l0065_certificate(cert.as_dict())
    save_l0065_certificate(cert)
    if not ok:
        raise RuntimeError(f"L-0065 certificate failed: {checks}")

    c7_path = Path("conjectures/active/C-0007.json")
    if c7_path.is_file():
        c7 = json.loads(c7_path.read_text(encoding="utf-8"))
        related = set(c7.get("related") or [])
        related.update(["L-0065", "C-R-0013"])
        c7["related"] = sorted(related)
        c7["notes"] = (
            f"Open all-IC. L-0065: one-pol greedy-second SOS C_ub_hi≈{b.C_ub_hi:.4f} "
            f"({len(b.support_radii or [])} shells, D={b.D}); proposed C-R-0013."
        )
        c7_path.write_text(json.dumps(c7, indent=2), encoding="utf-8")

    out = {
        "lemma_id": b.lemma_id,
        "D": b.D,
        "C_ub_float": b.C_ub_float,
        "C_ub_hi": b.C_ub_hi,
        "C_shor_sym": b.C_shor_sym,
        "closes_vs_Cdagger": b.closes_vs_Cdagger,
        "certificate_ok": ok,
        "checks": {k: bool(v) for k, v in checks.items()},
        "proposed_cr": b.proposed_cr,
        "not_a_clay_claim": True,
    }
    print(json.dumps(out, indent=2, ensure_ascii=True))
    return out


if __name__ == "__main__":
    main()
