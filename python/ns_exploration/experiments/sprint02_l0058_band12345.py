"""Sprint 02: L-0058 band {1..5} N5 Gram + certificate."""

from __future__ import annotations

from ns_exploration.conjectures.l0058_band12345_sos_gram_n5 import lemma_l0058, save_lemma_l0058
from ns_exploration.validation.l0058_certificate import build_l0058_certificate, save_l0058_certificate


def main() -> None:
    l58 = lemma_l0058()
    save_lemma_l0058(l58)
    cert = build_l0058_certificate()
    save_l0058_certificate(cert)
    print(f"L-0058: status={l58.status} C_ub_hi={l58.C_ub_hi:.6f}")
    print(f"  cert={cert.certificate_id}")


if __name__ == "__main__":
    main()
