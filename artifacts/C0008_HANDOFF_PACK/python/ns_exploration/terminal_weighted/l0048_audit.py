"""Audit classification of L-0048 value 25.925922..."""

from __future__ import annotations

from ns_exploration.conjectures.l0048_fullsym_techo import C_FULLSYM_FULL_DEALIAS
from ns_exploration.conjectures.l0048_matrixfree_fullsym import (
    sparse_C_fullsym,
    validate_sparse_vs_dense,
)


def l0048_value_classification(*, validate_dense: bool = True) -> dict:
    """
    Classification per spec §8.1:

    The frozen value is (a) approximate float64 value of sigma_max(M) from
    numpy.linalg.eigvalsh on Gram M@M.T — NOT (c) rigorous upper bound.
    """
    radii = tuple(range(1, 7))
    sparse = sparse_C_fullsym(24, radii)
    val: dict = {
        "sparse_C_band16": sparse["C_fullsym"],
        "sparse_D": sparse["D"],
    }
    dense_ok: bool | None = None
    rel_err: float | None = None
    if validate_dense:
        try:
            cmp = validate_sparse_vs_dense(24, radii)
            dense_ok = cmp["ok"]
            rel_err = cmp["rel_err"]
            val.update(cmp)
        except MemoryError:
            dense_ok = None
            rel_err = None
            val["dense_skipped"] = "MemoryError"

    return {
        "frozen_value": C_FULLSYM_FULL_DEALIAS,
        "classification": "a_approximate_float64_operator_norm",
        "not_certified_as": ["c_rigorous_upper", "d_interval", "e_analytic_lax"],
        "computation_path": "sparse_csr -> G=M@M.T dense -> eigvalsh -> sqrt -> *2√2",
        "evidence_level_computed": "N2",
        "analytic_majorant_holds": True,
        "note": (
            "True C_term(T) <= 2√2 sigma_max(M). Float sigma may underestimate "
            "or overestimate; Phase D replaces with Frobenius/MPFR upper bounds."
        ),
        "validate_band_16_rel_err": rel_err,
        "dense_sparse_agree_band16": dense_ok,
        "sparse_reference": val,
    }
