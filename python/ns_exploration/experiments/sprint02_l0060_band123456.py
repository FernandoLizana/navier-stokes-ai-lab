"""Sprint 02: L-0060 band {1..6} N5 Gram + certificate."""

from __future__ import annotations

from ns_exploration.conjectures.l0060_band123456_sos_gram_n5 import (
    GRAM_JSON,
    lemma_l0060,
    save_lemma_l0060,
)
from ns_exploration.validation.l0060_certificate import build_l0060_certificate, save_l0060_certificate


def main() -> None:
    gram_ok = GRAM_JSON.is_file()
    l60 = lemma_l0060(require_gram=gram_ok)
    save_lemma_l0060(l60)
    cert = build_l0060_certificate()
    save_l0060_certificate(cert)
    print(f"L-0060: status={l60.status} C_ub_hi={l60.C_ub_hi:.6f}")
    print(f"  cert={cert.certificate_id}")


if __name__ == "__main__":
    main()
