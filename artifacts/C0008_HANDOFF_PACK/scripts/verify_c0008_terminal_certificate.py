#!/usr/bin/env python3
"""Independent verifier for C-0008 terminal-weighted certificate (Phase D)."""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

from gmpy2 import mpfr

from ns_exploration.terminal_weighted.constants import FROZEN
from ns_exploration.terminal_weighted.intervals import i_star_interval
from ns_exploration.terminal_weighted.recompute_integral import (
    omega_T_hi_from_integral,
    shell_integral_hi_from_blocks,
)
from ns_exploration.terminal_weighted.weights import stokes_floor_hi


def _manifest_hash(manifest: dict) -> str:
    m = dict(manifest)
    m.pop("manifest_sha256", None)
    return hashlib.sha256(json.dumps(m, sort_keys=True, default=str).encode()).hexdigest()


def verify_certificate(path: str | Path, *, manifest_path: str | Path | None = None) -> tuple[bool, list[str]]:
    path = Path(path)
    msgs: list[str] = []
    if not path.is_file():
        return False, [f"Missing certificate: {path}"]

    cert = json.loads(path.read_text(encoding="utf-8"))
    ok = True

    if cert.get("claim") != "C-0008":
        ok = False
        msgs.append("claim != C-0008")

    if cert.get("status") not in ("candidate_certificate", "verified"):
        ok = False
        msgs.append("invalid status")

    params = cert.get("parameters", {})
    for key, val in [("nu", FROZEN.nu), ("E0", FROZEN.E0), ("T", FROZEN.T)]:
        if abs(float(params.get(key, 0)) - val) > 1e-15:
            ok = False
            msgs.append(f"parameter mismatch {key}")

    if int(params.get("N_max", 0)) != 24:
        ok = False
        msgs.append("N_max != 24")

    I_lo, I_hi = i_star_interval()
    istar = cert.get("i_star_interval", [])
    if len(istar) == 2:
        cert_I_lo = mpfr(istar[0])
        cert_I_hi = mpfr(istar[1])
        if cert_I_lo > I_lo * mpfr("1.0000000000000001") or cert_I_hi < I_hi:
            ok = False
            msgs.append("i_star_interval inconsistent with recomputation")

    floor = cert.get("stokes_floor_interval", [])
    sf = stokes_floor_hi(24)
    if len(floor) == 2 and float(floor[1]) < sf - 1e-9:
        ok = False
        msgs.append("stokes_floor_interval too low")

    mode_set = cert.get("mode_set", {})
    if mode_set.get("radii_mode") == "band":
        ok = False
        msgs.append("band radii pilot rejected for all-IC full-dealias claim")

    if mode_set.get("radii_mode") == "full":
        if mode_set.get("D") != FROZEN.D_full:
            ok = False
            msgs.append(f"D={mode_set.get('D')} != D_full={FROZEN.D_full}")
        if mode_set.get("n_shells", 0) != 87:
            ok = False
            msgs.append(f"n_shells={mode_set.get('n_shells')} != 87")

    src = cert.get("source_hashes", {})
    manifest_sha = src.get("shell_manifest")
    manifest: dict | None = None
    if manifest_path:
        mp = Path(manifest_path)
        if mp.is_file():
            manifest = json.loads(mp.read_text(encoding="utf-8"))
    elif manifest_sha:
        default_mp = Path("experiments/terminal_weighted/shell_manifest.json")
        if default_mp.is_file():
            manifest = json.loads(default_mp.read_text(encoding="utf-8"))

    if manifest_sha and manifest:
        recomputed = _manifest_hash(manifest)
        if recomputed != manifest_sha:
            ok = False
            msgs.append("shell_manifest hash mismatch")
        blocks = manifest.get("blocks", [])
        if blocks:
            prec = cert.get("arithmetic_backend", {}).get("precision_bits", 200)
            I_re = shell_integral_hi_from_blocks(blocks, prec=prec)
            route_a = cert.get("routes", {}).get("A", {})
            if route_a.get("integral_hi"):
                if mpfr(route_a["integral_hi"]) > I_re * mpfr("1.0000000000000001"):
                    ok = False
                    msgs.append("Route A integral_hi exceeds recomputation from blocks")

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
            omega_re = omega_T_hi_from_integral(mpfr(intv[1]))
            if omega_hi > omega_re * mpfr("1.0000000000000001"):
                ok = False
                msgs.append("omega_T_hi exceeds recomputation from integral")

    if cert.get("method") == "band_only":
        ok = False
        msgs.append("band-only certificate rejected for all-IC")

    backend = cert.get("arithmetic_backend", {})
    if backend.get("name") != "gmpy2_mpfr":
        ok = False
        msgs.append("arithmetic backend must be gmpy2_mpfr")

    if ok:
        msgs.append("PASS")
    else:
        msgs.append("FAIL")

    return ok, msgs


def main(argv: list[str] | None = None) -> int:
    argv = argv or sys.argv[1:]
    if not argv:
        print("Usage: python verify_c0008_terminal_certificate.py <certificate.json> [manifest.json]")
        return 2
    manifest_path = argv[1] if len(argv) > 1 else None
    ok, msgs = verify_certificate(argv[0], manifest_path=manifest_path)
    for m in msgs:
        print(m)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
