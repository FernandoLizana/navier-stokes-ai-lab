"""Build repair full/partial certificate from repair manifest."""

from __future__ import annotations

import json
from pathlib import Path

from ns_exploration.terminal_weighted.portable_paths import cert_root, data_root, rational_tag
from ns_exploration.validation.c0008_repair_band_cert import (
    build_repair_full_certificate,
    build_repair_full_partial_certificate,
)
from ns_exploration.validation.c0008_verify import verify_terminal_certificate

PY = Path(__file__).resolve().parents[1]
REPO = PY.parent


def _manifest_path(n: int, *, rational: bool = False) -> Path:
    tag = rational_tag(rational)
    for base in (data_root(), PY / "experiments/terminal_weighted", REPO / "experiments/terminal_weighted"):
        full = base / f"shell_manifest_repair_full_one_inf_n{n}{tag}.json"
        if full.is_file():
            return full
        partial = base / f"shell_manifest_repair_full_one_inf_n{n}{tag}_partial.json"
        if partial.is_file():
            return partial
    raise FileNotFoundError(f"repair manifest for n={n} rational={rational} not found")


def main() -> dict:
    import argparse

    p = argparse.ArgumentParser()
    p.add_argument("--n", type=int, default=24)
    p.add_argument("--prec", type=int, default=128)
    p.add_argument("--rational", action="store_true")
    args = p.parse_args()

    manifest_path = _manifest_path(args.n, rational=args.rational)
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    tag = rational_tag(args.rational)
    cert_path = cert_root() / f"CERT-C0008-repair-full-n{args.n}-one-inf{tag}.json"
    if manifest.get("repair_complete") and len(manifest.get("blocks", [])) >= manifest.get(
        "n_shells", 0
    ):
        cert = build_repair_full_certificate(manifest, prec=args.prec)
    else:
        cert = build_repair_full_partial_certificate(manifest, prec=args.prec)
    cert_path.parent.mkdir(parents=True, exist_ok=True)
    cert_path.write_text(json.dumps(cert, indent=2), encoding="utf-8")
    ok, issues = verify_terminal_certificate(
        cert, manifest=manifest, manifest_path=manifest_path, prec=args.prec
    )
    rep = {
        "manifest": str(manifest_path),
        "cert": str(cert_path),
        "verify_ok": ok,
        "issues": issues,
        "integral_hi": cert["integral_interval"][1],
        "omega_T_hi": cert["omega_T_interval"][1],
        "closes_strict": cert.get("closes_strict"),
    }
    print(rep)
    return rep


if __name__ == "__main__":
    main()
