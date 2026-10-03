"""Sprint: L-0053 N5 Gram on one-pol C-R-0008 (D=92)."""

from __future__ import annotations

import json

from ns_exploration.conjectures.l0053_onepol_sos_gram_n5 import lemma_l0053, save_lemma_l0053
from ns_exploration.validation.l0053_certificate import (
    build_l0053_certificate,
    save_l0053_certificate,
    verify_l0053_certificate,
)


def main() -> None:
    b = lemma_l0053()
    save_lemma_l0053(b)
    cert = build_l0053_certificate()
    save_l0053_certificate(cert)
    ok, checks = verify_l0053_certificate(cert.as_dict())
    out = {
        "lemma_id": b.lemma_id,
        "cert_id": cert.certificate_id,
        "C_ub_hi": b.C_ub_hi,
        "C_dagger": b.C_dagger,
        "C_fullsym": b.C_fullsym,
        "n_gram_blocks": b.n_gram_blocks,
        "verify_ok": ok,
        "checks": checks,
    }
    print(json.dumps(out, indent=2), flush=True)
    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
