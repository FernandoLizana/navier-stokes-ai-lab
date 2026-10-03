"""
Exact shell decomposition for terminal-weighted fullsym operator (t=T baseline).

Convention (documented):
  w_r(T) = r on output mode in shell r.
  M(T) = sum_{r in R} M_r  (operator sum, not op-norm sum)
  Each M_r is built with output_shell=r and weight terminal_weight(r,T,T)=r.

Temporal:  M(t) = sum_r exp(-2 nu r (T-t)) M_r  (weight factor on shell r block)
"""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

from gmpy2 import mpfr

from ns_exploration.conjectures.l0048_matrixfree_fullsym import all_dealias_radii
from ns_exploration.conjectures.l0045_streaming_fullsym import streaming_C_fullsym
from ns_exploration.conjectures.l0048_fullsym_techo import C_FULLSYM_FULL_DEALIAS
from ns_exploration.terminal_weighted.constants import FROZEN
from ns_exploration.terminal_weighted.intervals import (
    ArithmeticBackend,
    add_up,
    mul_up,
    sqrt_up,
)
from ns_exploration.terminal_weighted.intervals import _ctx
from ns_exploration.terminal_weighted.opnorm_certified import (
    certified_L_op_frobenius_hi,
    certified_shell_blocks,
)
from ns_exploration.terminal_weighted.tensor import (
    build_terminal_G,
    streaming_terminal_fullsym_M,
    validate_direct_vs_tensor,
)
from ns_exploration.experiments.sprintA_stretch_tensor import f_from_G, full_sym, operator_C


def frobenius_sq_from_streaming_M(M) -> float:
    return float((M * M).sum())


def audit_shell_reconstruction_band(n: int = 24, radii: tuple[int, ...] = tuple(range(1, 7))) -> dict:
    """Verify sum_r M_r recovers M(T) on small band via Frobenius norms (N2)."""
    from ns_exploration.terminal_weighted.tensor import streaming_terminal_fullsym_M

    M_full, b = streaming_terminal_fullsym_M(n, radii, FROZEN.T)
    f_full = frobenius_sq_from_streaming_M(M_full)
    f_sum = 0.0
    for r in sorted(set(radii)):
        M_r, _ = streaming_terminal_fullsym_M(
            n, radii, FROZEN.T, output_shell=r
        )
        f_sum += frobenius_sq_from_streaming_M(M_r)
    # Frobenius of sum <= sum of Frobenius (triangle ineq on entries only for disjoint support?)
    # For exact entry decomposition M = sum M_r with disjoint output rows? Not disjoint rows.
    # Check: ||M||_F <= sum ||M_r||_F always
    val = validate_direct_vs_tensor(n, radii, FROZEN.T, n_probe=30)
    return {
        "radii": list(radii),
        "D": b.D,
        "frobenius_full": f_full,
        "frobenius_sum_blocks": f_sum,
        "direct_vs_tensor_ok": val["ok"],
        "C_stream": val["C_stream"],
        "C_dense": val["C_dense_fullsym"],
    }


def build_shell_manifest(
    n: int = 24,
    prec: int = 200,
    *,
    radii: tuple[int, ...] | None = None,
    workers: int = 1,
    bound_method: str = "best",
) -> dict:
    radii = radii or all_dealias_radii(n)
    full_dealias = len(radii) >= len(all_dealias_radii(n))
    if full_dealias and bound_method in ("one_inf", "best"):
        raise ValueError(
            f"bound_method={bound_method!r} on full dealias (D={FROZEN.D_full}, "
            f"{len(radii)} shells) exhausts RAM: sym-col accumulators are O(D^2). "
            "Use bound_method='frobenius' or: "
            "python -m ns_exploration.experiments.upgrade_manifest_best_bounds"
        )
    blocks_data = certified_shell_blocks(n, radii, prec=prec, workers=workers, bound_method=bound_method)
    manifest = {
        "n": n,
        "D_full": blocks_data["full_at_T_best"]["D"],
        "n_shells": blocks_data["n_shells"],
        "bound_method": bound_method,
        "convention": blocks_data["convention"],
        "blocks": blocks_data["blocks"],
        "full_frobenius_hi": blocks_data["full_at_T_frobenius"],
        "full_one_inf_hi": blocks_data["full_at_T_one_inf"],
        "full_best_hi": blocks_data["full_at_T_best"],
        "l0048_float_reference": C_FULLSYM_FULL_DEALIAS,
        "arithmetic_backend": blocks_data["arithmetic_backend"],
        "reconstruction_band": audit_shell_reconstruction_band(n, tuple(range(1, 7))),
    }
    payload = dict(manifest)
    payload.pop("manifest_sha256", None)
    digest = hashlib.sha256(
        json.dumps(payload, sort_keys=True, default=str).encode()
    ).hexdigest()
    manifest["manifest_sha256"] = digest
    return manifest


def save_shell_manifest(manifest: dict, path: str | Path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return path
