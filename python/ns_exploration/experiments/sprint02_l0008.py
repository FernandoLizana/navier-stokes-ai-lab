"""Sprint: L-0008 cascade (E_H=0) + N5 cosh certificate."""

from __future__ import annotations

import json
from pathlib import Path

from ns_exploration.conjectures.l0008_cascade import lemma_l0008, save_lemma_l0008
from ns_exploration.validation.l0008_certificate import (
    build_l0008_cosh_certificate,
    save_l0008_cosh_certificate,
    verify_l0008_cosh_certificate,
)


def main() -> dict:
    table = []
    for n, kind, k in (
        (12, "euclidean", 2),
        (12, "euclidean", 3),
        (12, "linf", 4),
        (16, "linf", 4),
    ):
        b = lemma_l0008(n=n, k_ic=k, ic_kind=kind)
        table.append(
            {
                "n": n,
                "ic_kind": kind,
                "k_ic": k,
                "cosh_cap": b.cosh_cap,
                "L0007_cap": b.L0007_cap,
                "L0003_cap": b.L0003_cap,
                "best": b.best_cap,
                "which": b.which_best,
                "improvement_vs_L0007": b.L0007_cap / b.cosh_cap if b.cosh_cap > 0 else None,
                "proves_cs0002": b.proves_cs0002,
            }
        )

    canon12 = lemma_l0008(n=12, k_ic=2, ic_kind="euclidean", cs0002_M=None)
    save_lemma_l0008(canon12, "conjectures/proved_restricted/L-0008.json")

    cert = build_l0008_cosh_certificate(n=12, k_ic=2, ic_kind="euclidean")
    ok, checks = verify_l0008_cosh_certificate(cert.as_dict())
    save_l0008_cosh_certificate(cert, "certificates/CERT-L0008-cosh-euclidean-k2-N12.json")

    cs = lemma_l0008(n=16, k_ic=4, ic_kind="linf")
    out = {
        "table": table,
        "n12_eucl_k2": {
            "cosh_cap": canon12.cosh_cap,
            "L0007_cap": canon12.L0007_cap,
            "ratio_improvement": canon12.L0007_cap / canon12.cosh_cap,
        },
        "cs0002_class_n16": {
            "best_cap": cs.best_cap,
            "which": cs.which_best,
            "cosh_cap": cs.cosh_cap,
            "proves": cs.proves_cs0002,
            "gap_vs_M": cs.best_cap / float(cs.cs0002_M),
        },
        "certificate_verified": ok,
        "certificate_bound_hi": cert.Omega_bound_hi,
        "certificate_checks": {k: bool(v) for k, v in checks.items()},
        "not_a_clay_claim": True,
    }
    Path("reports").mkdir(exist_ok=True)
    print(json.dumps(out, indent=2, ensure_ascii=True))
    return out


if __name__ == "__main__":
    main()
