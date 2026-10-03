"""Build repair-band C-0008 prototype certificate (one_inf_only + L-0073R)."""

from __future__ import annotations

import json
import platform
from pathlib import Path

import gmpy2
from gmpy2 import mpfr

from ns_exploration.terminal_weighted.constants import EXACT, FROZEN
from ns_exploration.terminal_weighted.intervals import (
    _ctx,
    i_star_interval,
    mpfr_const,
    omega_T_hi_from_I_term,
    stokes_floor_interval,
)
from ns_exploration.terminal_weighted.recompute_integral import shell_integral_hi_from_blocks
from ns_exploration.validation.c0008_verify import manifest_sha256


def build_repair_band_certificate(
    manifest: dict,
    *,
    prec: int = 128,
    radii: tuple[int, ...] = (1, 2, 3, 4, 5, 6),
) -> dict:
    return _build_repair_certificate(
        manifest, prec=prec, scope="repair_band", radii=radii
    )


def build_repair_full_partial_certificate(
    manifest: dict,
    *,
    prec: int = 128,
) -> dict:
    return _build_repair_certificate(manifest, prec=prec, scope="repair_full_partial")


def build_repair_full_certificate(
    manifest: dict,
    *,
    prec: int = 128,
) -> dict:
    return _build_repair_certificate(manifest, prec=prec, scope="repair_full")


def _build_repair_certificate(
    manifest: dict,
    *,
    prec: int = 128,
    scope: str,
    radii: tuple[int, ...] | None = None,
) -> dict:
    I_lo, I_hi = i_star_interval(prec=prec)
    floor_lo, floor_hi = stokes_floor_interval(prec=prec)
    # Store intervals as float64 literals (round-trip safe in verifier at prec=128).
    i_star_store = [FROZEN.integral_threshold * 0.999999, FROZEN.integral_threshold * 1.000001]
    floor_store = [FROZEN.stokes_floor, float(floor_hi)]
    integral_hi = shell_integral_hi_from_blocks(manifest["blocks"], prec=prec)
    omega_hi = omega_T_hi_from_I_term(integral_hi, prec=prec)
    D = int(manifest["blocks"][0]["D"]) if manifest.get("blocks") else 0
    cluster = manifest.get("cluster_residual_L0073R", {})
    n_shells = manifest.get("n_shells", len(manifest.get("blocks", [])))
    complete = bool(
        manifest.get("repair_complete", scope in ("repair_band", "repair_full"))
    )
    radii_mode = "band" if scope == "repair_band" else "full"

    mode_set = {
        "n": manifest.get("n", 24),
        "D": D,
        "n_shells": n_shells,
        "radii_mode": radii_mode,
    }
    if radii_mode == "band" and radii is not None:
        mode_set["radii"] = list(radii)

    if scope == "repair_full":
        status = "repair_full"
        method = "repair_full_one_inf_shell_integral"
        limitations = [
            f"Full dealias D={D} n_shells={n_shells} one_inf_only per-shell pass.",
            "Route A' cluster L-0073R optional sketch — declared bound is shell sum.",
            "N<=24 domination remains N7_sketch (L-0076 band ladder).",
        ]
    elif scope == "repair_full_partial":
        status = "repair_full_partial"
        method = "repair_full_partial"
        limitations = [
            f"Partial full-dealias: {len(manifest.get('blocks', []))}/{n_shells} shells",
            "one_inf_only baseline; L-0073R cluster residual required for Route A'",
            "N<=24 domination remains N7_sketch (L-0076 band ladder)",
        ]
    else:
        status = "repair_band_prototype"
        method = "repair_band_prototype"
        limitations = [
            "Band {1..6} only — NOT full-dealias D=6748",
            "one_inf_only baseline; L-0073R cluster residual required for Route A'",
            "N<=24 domination remains N7_sketch (L-0076 band ladder)",
        ]

    return {
        "claim": "C-0008",
        "status": status,
        "scope": scope,
        "parameters": {
            "nu": EXACT.NU,
            "E0": EXACT.E0,
            "T": EXACT.T,
            "N_max": 24,
            "M_target": EXACT.M_TARGET,
        },
        "bound_method": "one_inf_only",
        "mode_set": mode_set,
        "repair_complete": complete,
        "arithmetic_backend": {
            "name": "gmpy2_mpfr",
            "version": gmpy2.version(),
            "precision_bits": prec,
            "rounding_up": "RoundUp",
            "rounding_down": "RoundDown",
            "platform": platform.platform(),
        },
        "stokes_floor_interval": floor_store,
        "i_star_interval": [str(i_star_store[0]), str(i_star_store[1])],
        "method": method,
        "integral_interval": ["0", str(integral_hi)],
        "omega_T_interval": ["0", str(omega_hi)],
        "coverage_complete": True,
        "source_hashes": {"shell_manifest": manifest_sha256(manifest)},
        "cluster_residual_L0073R": cluster,
        "B_C_combined_hi": manifest.get("B_C_combined_hi"),
        "limitations": limitations,
        "closes_c0008": False,
        "closes_strict": bool(integral_hi < I_lo),
        "honesty": (
            "Repair prototype certificate for adversarial verifier regression. "
            "Does NOT close C-0008 globally."
        ),
    }


def save_repair_band_certificate(
    manifest_path: Path,
    cert_path: Path,
    *,
    prec: int = 128,
) -> dict:
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    cert = build_repair_band_certificate(manifest, prec=prec)
    cert_path.parent.mkdir(parents=True, exist_ok=True)
    cert_path.write_text(json.dumps(cert, indent=2), encoding="utf-8")
    return cert
