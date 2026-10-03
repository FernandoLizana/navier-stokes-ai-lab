"""Sprint: hygiene + L-0007 bridge (shell IC → full grid) + N5 exp certificate."""

from __future__ import annotations

import json
from pathlib import Path

from ns_exploration.conjectures.l0007_shell_ic_fullgrid import lemma_l0007, save_lemma_l0007
from ns_exploration.validation.l0007_certificate import (
    build_l0007_exp_certificate,
    save_l0007_exp_certificate,
    verify_l0007_exp_certificate,
)


def main() -> dict:
    table = []
    for n in (12, 16):
        for kind, k in (("euclidean", 2), ("euclidean", 3), ("linf", 4)):
            b = lemma_l0007(n=n, k_ic=k, ic_kind=kind)
            table.append(
                {
                    "n": n,
                    "ic_kind": kind,
                    "k_ic": k,
                    "K_ic2": b.K_ic_squared,
                    "M_full": b.M_full,
                    "which": b.which_best,
                    "best_cap": b.best_cap,
                    "L0003_cap": b.L0003_cap,
                    "proves_cs0002": b.proves_cs0002,
                }
            )

    # Canonical lemma record for C-S-0002 class
    canon = lemma_l0007(n=16, k_ic=4, ic_kind="linf")
    save_lemma_l0007(canon, "conjectures/proved_restricted/L-0007.json")

    # N5 certificate for the sharpest short-time exp case at N=12
    cert = build_l0007_exp_certificate(n=12, k_ic=2, ic_kind="euclidean")
    ok, checks = verify_l0007_exp_certificate(cert.as_dict())
    save_l0007_exp_certificate(cert, "certificates/CERT-L0007-exp-euclidean-k2-N12.json")

    out = {
        "table": table,
        "canonical_cs0002_class": {
            "best_cap": canon.best_cap,
            "which": canon.which_best,
            "proves_cs0002": canon.proves_cs0002,
            "cs0002_M": canon.cs0002_M,
            "gap_factor": canon.best_cap / float(canon.cs0002_M),
        },
        "certificate": cert.as_dict(),
        "certificate_verified": ok,
        "certificate_checks": {k: bool(v) for k, v in checks.items()},
        "hygiene": {
            "active_should_be": ["C-0002", "C-S-0002"],
            "roadmap_updated": True,
        },
        "not_a_clay_claim": True,
    }
    Path("reports").mkdir(exist_ok=True)
    print(json.dumps(out, indent=2, ensure_ascii=True))
    return out


if __name__ == "__main__":
    main()
