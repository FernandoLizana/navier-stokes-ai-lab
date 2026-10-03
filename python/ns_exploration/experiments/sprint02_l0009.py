"""Sprint: L-0009 two-scale embedding + N5 certificate."""

from __future__ import annotations

import json
from pathlib import Path

from ns_exploration.conjectures.l0009_twoscale import lemma_l0009, save_lemma_l0009
from ns_exploration.validation.l0009_certificate import (
    build_l0009_certificate,
    save_l0009_certificate,
    verify_l0009_certificate,
)


def main() -> dict:
    table = []
    for n, kind, k in (
        (12, "euclidean", 2),
        (12, "euclidean", 3),
        (12, "linf", 4),
        (16, "linf", 4),
    ):
        b = lemma_l0009(n=n, k_ic=k, ic_kind=kind)
        table.append(
            {
                "n": n,
                "ic_kind": kind,
                "k_ic": k,
                "M_L": b.M_L,
                "M_H": b.M_H,
                "closes": b.two_scale_closes,
                "cap9": b.two_scale_cap,
                "L0008": b.L0008_cap,
                "best": b.best_cap,
                "which": b.which_best,
                "proves_cs0002": b.proves_cs0002,
            }
        )

    canon = lemma_l0009(n=12, k_ic=2, ic_kind="euclidean", cs0002_M=None)
    save_lemma_l0009(canon, "conjectures/proved_restricted/L-0009.json")

    cert = build_l0009_certificate()
    ok, checks = verify_l0009_certificate(cert.as_dict())
    save_l0009_certificate(cert, "certificates/CERT-L0009-twoscale-euclidean-k2-N12.json")

    cs = lemma_l0009(n=16, k_ic=4, ic_kind="linf")
    out = {
        "table": table,
        "n12_eucl_k2": {
            "cap9": canon.two_scale_cap,
            "L0008": canon.L0008_cap,
            "improvement_vs_L0008": canon.L0008_cap / float(canon.two_scale_cap),
        },
        "cs0002_class": {
            "closes": cs.two_scale_closes,
            "best": cs.best_cap,
            "which": cs.which_best,
            "gap": cs.best_cap / float(cs.cs0002_M),
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
