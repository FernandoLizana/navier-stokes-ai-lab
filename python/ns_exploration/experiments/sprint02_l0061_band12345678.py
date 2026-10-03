"""Sprint 02: L-0061 band {1..8} N5 + certificate."""

from __future__ import annotations

from ns_exploration.conjectures.l0061_band12345678_sos_n5 import GRAM_JSON, lemma_l0061, save_lemma_l0061
from ns_exploration.validation.l0061_certificate import build_l0061_certificate, save_l0061_certificate


def main() -> None:
    l61 = lemma_l0061(require_gram=GRAM_JSON.is_file())
    save_lemma_l0061(l61)
    cert = build_l0061_certificate()
    save_l0061_certificate(cert)
    print(f"L-0061: status={l61.status} C_ub_hi={l61.C_ub_hi:.6f} solver={l61.solver}")
    print(f"  cert={cert.certificate_id}")


if __name__ == "__main__":
    main()
