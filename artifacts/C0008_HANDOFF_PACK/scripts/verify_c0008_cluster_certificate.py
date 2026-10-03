#!/usr/bin/env python3
"""Independent verifier for C-0008 Route A' cluster certificate."""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

from gmpy2 import mpfr

from ns_exploration.terminal_weighted.cluster_bounds import cluster_integral_certificate
from ns_exploration.terminal_weighted.constants import FROZEN
from ns_exploration.terminal_weighted.intervals import i_star_interval
from ns_exploration.terminal_weighted.recompute_integral import omega_T_hi_from_integral


def _manifest_hash(manifest: dict) -> str:
    m = dict(manifest)
    m.pop("manifest_sha256", None)
    return hashlib.sha256(json.dumps(m, sort_keys=True, default=str).encode()).hexdigest()


def _recompute_integral_hi(manifest: dict, *, prec: int) -> mpfr:
    total = mpfr(0, precision=prec)
    for blk in manifest.get("blocks", []):
        total += mpfr(blk["contrib_hi"], precision=prec)
    return total


def verify_certificate(
    cert_path: str | Path, *, manifest_path: str | Path | None = None
) -> tuple[bool, list[str]]:
    cert_path = Path(cert_path)
    msgs: list[str] = []
    if not cert_path.is_file():
        return False, [f"Missing certificate: {cert_path}"]

    cert = json.loads(cert_path.read_text(encoding="utf-8"))
    ok = True

    if cert.get("claim") != "C-0008":
        ok = False
        msgs.append("claim != C-0008")

    if cert.get("method") != "route_a_prime_cluster_best":
        ok = False
        msgs.append("method != route_a_prime_cluster_best")

    params = cert.get("parameters", {})
    for key, val in [("nu", FROZEN.nu), ("E0", FROZEN.E0), ("T", FROZEN.T)]:
        if abs(float(params.get(key, 0)) - val) > 1e-15:
            ok = False
            msgs.append(f"parameter mismatch {key}")

    manifest_sha = cert.get("source_hashes", {}).get("cluster_manifest")
    manifest: dict | None = None
    if manifest_path:
        mp = Path(manifest_path)
        if mp.is_file():
            manifest = json.loads(mp.read_text(encoding="utf-8"))
    elif manifest_sha:
        default = Path("experiments/terminal_weighted/shell_manifest_cluster_best_full.json")
        if default.is_file():
            manifest = json.loads(default.read_text(encoding="utf-8"))

    prec = cert.get("arithmetic_backend", {}).get("precision_bits", 128)
    if manifest_sha and manifest:
        if _manifest_hash(manifest) != manifest_sha:
            ok = False
            msgs.append("cluster_manifest hash mismatch")
        I_re = _recompute_integral_hi(manifest, prec=int(prec))
        cert_I = mpfr(cert.get("integral_interval", ["0", "0"])[1])
        if cert_I > I_re * mpfr("1.0000000000000001"):
            ok = False
            msgs.append("integral_hi exceeds recomputation from cluster blocks")

    I_lo, _ = i_star_interval(prec=int(prec))
    intv = cert.get("integral_interval", [])
    if len(intv) != 2:
        ok = False
        msgs.append("missing integral_interval")
    else:
        I_term_hi = mpfr(intv[1])
        if I_term_hi >= I_lo:
            ok = False
            msgs.append(f"FAIL: I_term_hi={I_term_hi} >= I_star_lo={I_lo}")

    omega = cert.get("omega_T_interval", [])
    if len(omega) == 2:
        omega_hi = mpfr(omega[1])
        if omega_hi >= mpfr(FROZEN.c0008_target_M):
            ok = False
            msgs.append("omega_T_hi >= M_target")
        if len(intv) == 2:
            omega_re = omega_T_hi_from_integral(mpfr(intv[1]), prec=int(prec))
            if omega_hi > omega_re * mpfr("1.0000000000000001"):
                ok = False
                msgs.append("omega_T_hi exceeds recomputation from integral")

    if cert.get("closes_strict"):
        ok = False
        msgs.append("closes_strict must be false for open certificate")

    if ok:
        msgs.append("PASS")
    else:
        msgs.append("FAIL (expected for open C-0008)")

    return ok, msgs


def main(argv: list[str] | None = None) -> int:
    argv = argv or sys.argv[1:]
    if not argv:
        print(
            "Usage: python verify_c0008_cluster_certificate.py <certificate.json> [manifest.json]"
        )
        return 2
    manifest_path = argv[1] if len(argv) > 1 else None
    ok, msgs = verify_certificate(argv[0], manifest_path=manifest_path)
    for m in msgs:
        print(m)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
