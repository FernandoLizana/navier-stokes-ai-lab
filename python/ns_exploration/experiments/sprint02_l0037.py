"""Sprint: L-0037 shell-diff / bi-shell Shor techo."""

from __future__ import annotations

import json
from pathlib import Path

from ns_exploration.conjectures.l0037_shell_diff_stretch import lemma_l0037, save_lemma_l0037
from ns_exploration.validation.l0037_certificate import (
    build_l0037_certificate,
    save_l0037_certificate,
    verify_l0037_certificate,
)


def main() -> dict:
    b = lemma_l0037(n=24, iters=16, max_pairs=100)
    save_lemma_l0037(b)
    cert = build_l0037_certificate(n=24)
    ok, checks = verify_l0037_certificate(cert.as_dict())
    save_l0037_certificate(cert)
    if not ok:
        raise RuntimeError(f"L-0037 certificate failed: {checks}")

    c7_path = Path("conjectures/active/C-0007.json")
    c7 = json.loads(c7_path.read_text(encoding="utf-8"))
    c7["notes"] = (
        f"Open all-IC. Subclasses C-R-0002/0003. L-0037: mono stretch=0; "
        f"bi-shell Shor near{b.near_pair} C≳{b.C_near_pair:.2f}≤C_† but far "
        f"max≳{b.C_bishell_shor_max_sample:.1f}>C_†. Need multi-shell SOS/rank."
    )
    related = set(c7.get("related") or [])
    related.add("L-0037")
    c7["related"] = sorted(related)
    c7_path.write_text(json.dumps(c7, indent=2), encoding="utf-8")

    out = {
        "mono_stretch_zero": b.mono_stretch_zero,
        "C_near_pair": b.C_near_pair,
        "near_pair": b.near_pair,
        "C_far_pair": b.C_far_pair,
        "far_pair": b.far_pair,
        "C_bishell_shor_max_sample": b.C_bishell_shor_max_sample,
        "C_dagger": b.C_dagger,
        "n_pairs_scanned": b.n_pairs_scanned,
        "closes_c0007": False,
        "certificate_ok": ok,
        "checks": {k: bool(v) for k, v in checks.items()},
        "not_a_clay_claim": True,
    }
    Path("reports").mkdir(exist_ok=True)
    Path("experiments/exploratory/sprint02_l0037").mkdir(parents=True, exist_ok=True)
    Path("experiments/exploratory/sprint02_l0037/summary.json").write_text(
        json.dumps(out, indent=2), encoding="utf-8"
    )
    Path("reports/SPRINT_02_L0037.md").write_text(
        f"""# Sprint — L-0037 shell-diff / bi-shell Shor

| Item | Value |
|------|-------|
| Mono-radial stretch | 0 (exact) |
| Near pair {b.near_pair} C_Shor ≳ | ≈{b.C_near_pair:.4f} |
| Far/sample max C_Shor ≳ | ≈{b.C_bishell_shor_max_sample:.4f} |
| C_† | ≈{b.C_dagger:.4f} |
| C-0007 | still open |

FINITE Galerkin only. Not continuum. Not Clay.
""",
        encoding="utf-8",
    )
    print(json.dumps(out, indent=2, ensure_ascii=True))
    return out


if __name__ == "__main__":
    main()
