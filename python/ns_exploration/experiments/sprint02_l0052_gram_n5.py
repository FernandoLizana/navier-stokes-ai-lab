"""Sprint: L-0052 N5 Gram rationalization on band {1,2,3}."""

from __future__ import annotations

import json

from ns_exploration.conjectures.l0052_sos_gram_n5 import lemma_l0052, save_lemma_l0052
from ns_exploration.validation.l0052_certificate import (
    build_l0052_certificate,
    save_l0052_certificate,
    verify_l0052_certificate,
)


def main() -> None:
    b = lemma_l0052()
    save_lemma_l0052(b)
    cert = build_l0052_certificate()
    save_l0052_certificate(cert)
    ok, checks = verify_l0052_certificate(cert.as_dict())
    out = {
        "lemma_id": b.lemma_id,
        "cert_id": cert.certificate_id,
        "C_ub_hi": b.C_ub_hi,
        "C_dagger": b.C_dagger,
        "n_gram_blocks": b.n_gram_blocks,
        "verify_ok": ok,
        "checks": checks,
    }
    print(json.dumps(out, indent=2), flush=True)
    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
