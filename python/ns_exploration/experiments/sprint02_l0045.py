"""Sprint: L-0045 streaming physical fullsym + C-R-0010."""

from __future__ import annotations

import json
from pathlib import Path

from ns_exploration.conjectures.l0045_physical_band import (
    lemma_l0045,
    save_cr0010,
    save_lemma_l0045,
)
from ns_exploration.validation.l0045_certificate import (
    build_l0045_certificate,
    save_l0045_certificate,
    verify_l0045_certificate,
)


def main() -> dict:
    b = lemma_l0045(n=24)
    save_lemma_l0045(b)
    save_cr0010(b)
    cert = build_l0045_certificate(n=24, bound=b)
    ok, checks = verify_l0045_certificate(cert.as_dict())
    save_l0045_certificate(cert)
    if not ok:
        raise RuntimeError(f"L-0045 certificate failed: {checks}")

    c7_path = Path("conjectures/active/C-0007.json")
    c7 = json.loads(c7_path.read_text(encoding="utf-8"))
    c7["notes"] = (
        f"Open all-IC. Closed: C-R-0002..0010. L-0045: physical fullsym "
        f"r≤25 C={b.C_fullsym_rmax25:.4f}≤C_† (C-R-0010); r≤26 C={b.C_fullsym_rmax26:.4f}>C_†. "
        "Full dealias (D≈6748) still open."
    )
    related = set(c7.get("related") or [])
    related.update(["L-0045", "C-R-0010"])
    c7["related"] = sorted(related)
    c7_path.write_text(json.dumps(c7, indent=2), encoding="utf-8")

    out = {
        "support_cr0010": b.support_cr0010,
        "C_fullsym_rmax25": b.C_fullsym_rmax25,
        "C_fullsym_rmax26": b.C_fullsym_rmax26,
        "C_dagger": b.C_dagger,
        "closes_cr0010": True,
        "closes_c0007_all_ic": False,
        "certificate_ok": ok,
        "checks": {k: bool(v) for k, v in checks.items()},
        "not_a_clay_claim": True,
    }
    Path("reports").mkdir(exist_ok=True)
    Path("experiments/exploratory/sprint02_l0045").mkdir(parents=True, exist_ok=True)
    Path("experiments/exploratory/sprint02_l0045/summary.json").write_text(
        json.dumps(out, indent=2), encoding="utf-8"
    )
    Path("reports/SPRINT_02_L0045.md").write_text(
        f"""# Sprint — L-0045 streaming physical fullsym + C-R-0010

| Item | Value |
|------|-------|
| **C-R-0010** (Hermitian \\|k\\|²≤25, {len(b.support_cr0010 or [])} shells) | C≈{b.C_fullsym_rmax25:.6f} ≤ C_† |
| Techo r≤26 | C≈{b.C_fullsym_rmax26:.6f} > C_† |
| C-0007 all-IC | still open (full dealias D≈6748) |

FINITE Galerkin only. Not continuum. Not Clay.
""",
        encoding="utf-8",
    )
    print(json.dumps(out, indent=2, ensure_ascii=True))
    return out


if __name__ == "__main__":
    main()
