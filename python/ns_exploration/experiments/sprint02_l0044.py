"""Sprint: L-0044 physical Hermitian fullsym + C-R-0009."""

from __future__ import annotations

import json
from pathlib import Path

from ns_exploration.conjectures.l0044_physical_fullsym import (
    lemma_l0044,
    save_cr0009,
    save_lemma_l0044,
)
from ns_exploration.validation.l0044_certificate import (
    build_l0044_certificate,
    save_l0044_certificate,
    verify_l0044_certificate,
)


def main() -> dict:
    b = lemma_l0044(n=24)
    save_lemma_l0044(b)
    save_cr0009(b)
    cert = build_l0044_certificate(n=24)
    ok, checks = verify_l0044_certificate(cert.as_dict())
    save_l0044_certificate(cert)
    if not ok:
        raise RuntimeError(f"L-0044 certificate failed: {checks}")

    c7_path = Path("conjectures/active/C-0007.json")
    c7 = json.loads(c7_path.read_text(encoding="utf-8"))
    c7["notes"] = (
        f"Open all-IC. Closed: C-R-0002..0009. L-0044: physical Hermitian "
        f"fullsym on {{1..6}} C={b.C_fullsym_123456:.4f}≤C_† (C-R-0009); "
        f"(i,j)-Sym still {b.C_sym_123456:.4f}>C_†. FFT↔tensor err={b.fft_tensor_err:.2e}. "
        "All-IC needs full-mask physical fullsym (memory)."
    )
    related = set(c7.get("related") or [])
    related.update(["L-0044", "C-R-0009"])
    c7["related"] = sorted(related)
    c7_path.write_text(json.dumps(c7, indent=2), encoding="utf-8")

    out = {
        "support_cr0009": b.support_cr0009,
        "C_fullsym_123456": b.C_fullsym_123456,
        "C_sym_123456": b.C_sym_123456,
        "fft_tensor_err": b.fft_tensor_err,
        "C_dagger": b.C_dagger,
        "closes_cr0009": True,
        "closes_c0007_all_ic": False,
        "certificate_ok": ok,
        "checks": {k: bool(v) for k, v in checks.items()},
        "not_a_clay_claim": True,
    }
    Path("reports").mkdir(exist_ok=True)
    Path("experiments/exploratory/sprint02_l0044").mkdir(parents=True, exist_ok=True)
    Path("experiments/exploratory/sprint02_l0044/summary.json").write_text(
        json.dumps(out, indent=2), encoding="utf-8"
    )
    Path("reports/SPRINT_02_L0044.md").write_text(
        f"""# Sprint — L-0044 physical Hermitian fullsym + C-R-0009

| Item | Value |
|------|-------|
| **C-R-0009** (Hermitian shells {{1..6}}) | C_fullsym≈{b.C_fullsym_123456:.6f} ≤ C_† |
| (i,j)-only Sym (same tensor) | ≈{b.C_sym_123456:.6f} > C_† |
| FFT ↔ tensor max rel err | {b.fft_tensor_err:.2e} |
| C-0007 all-IC | still open |

Key audit: physical triad uses factor `i`; full symmetrization then closes {{1..6}}.

FINITE Galerkin only. Not continuum. Not Clay.
""",
        encoding="utf-8",
    )
    print(json.dumps(out, indent=2, ensure_ascii=True))
    return out


if __name__ == "__main__":
    main()
