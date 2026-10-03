"""Audit: does N=24 full-dealias majorant cover all N<=24?"""

from __future__ import annotations

import json
from pathlib import Path

from gmpy2 import mpfr

from ns_exploration.conjectures.l0021_quartic_shell import shell_modes_by_r
from ns_exploration.terminal_weighted.opnorm_certified import certified_shell_blocks
from ns_exploration.terminal_weighted.recompute_integral import shell_integral_hi_from_blocks
from ns_exploration.validation.l0003_certificate import full_dealias_exact_stats


def wavevector_shells_subset(n_small: int, n_large: int) -> dict:
    by_s = shell_modes_by_r(n_small)
    by_l = shell_modes_by_r(n_large)
    ks_s = {k for modes in by_s.values() for k in modes}
    ks_l = {k for modes in by_l.values() for k in modes}
    return {
        "n_small": n_small,
        "n_large": n_large,
        "modes_small": len(ks_s),
        "modes_large": len(ks_l),
        "subset": ks_s.issubset(ks_l),
        "only_large": len(ks_l - ks_s),
    }


def embed_basis_index_check(n_small: int, n_large: int) -> dict:
    """
    L-0076: each wavevector k in n_small appears in n_large shell lists with same k.
    Index alignment in Hermitian basis still requires explicit embedding map (N7 open).
    """
    by_s = shell_modes_by_r(n_small)
    by_l = shell_modes_by_r(n_large)
    ks_s = {k for modes in by_s.values() for k in modes}
    ks_l = {k for modes in by_l.values() for k in modes}
    missing = sorted(ks_s - ks_l)
    return {
        "n_small": n_small,
        "n_large": n_large,
        "wavevectors_small": len(ks_s),
        "missing_in_large": len(missing),
        "all_wavevectors_embedded": len(missing) == 0,
        "sample_missing": [list(m) for m in missing[:5]],
        "basis_index_map": "NOT_IMPLEMENTED — requires Hermitian basis sprintB index lookup",
        "evidence_level": "N7_sketch_empirical",
    }


def band_majorant_ladder(
    radii: tuple[int, ...] = (1, 2, 3),
    *,
    prec: int = 128,
) -> dict:
    """
    Empirical monotonicity: on fixed shell band, C_term sum at n=24 should
    majorize smaller grids (mode-set inclusion on band).
    """
    ladder = []
    for n in (12, 16, 20, 24):
        data = certified_shell_blocks(
            n, radii, prec=prec, workers=1, bound_method="one_inf_only"
        )
        total = mpfr(0, precision=prec)
        for blk in data["blocks"]:
            total += mpfr(blk["C_term_hi"], precision=prec)
        ladder.append(
            {
                "n": n,
                "D": data["blocks"][0]["D"] if data["blocks"] else 0,
                "C_term_sum_hi": str(total),
            }
        )
    vals = [float(r["C_term_sum_hi"]) for r in ladder]
    monotone = all(vals[i] <= vals[i + 1] + 1e-6 for i in range(len(vals) - 1))
    return {
        "radii": list(radii),
        "prec": prec,
        "ladder": ladder,
        "monotone_non_decreasing": monotone,
        "evidence_level": "N7_sketch_empirical",
    }


def n_domination_analysis(n_max: int = 24, *, run_band_ladder: bool = False) -> dict:
    """
    Option 1 (conditional): If majorant C_term(N) is computed as sup over all
    dealias modes at grid n, then for n' <= n the embedded subspace is contained
    in the larger mode set (wavevectors subset verified). Any field on n' embeds
    in n by zero-padding; stretch on embedded field equals n' restriction.
    The sup over z on R^D_n contains sup over z supported on embedded subspace,
    hence C_term(N_max) >= C_term|_{embedded n'} (majorant monotone in mode inclusion).

    Requires: same Fourier normalization, same dealias 2/3 convention per n.
    """
    subs = [wavevector_shells_subset(n, n_max) for n in [12, 16, 20, 24] if n <= n_max]
    K2, M = full_dealias_exact_stats(n_max)
    out = {
        "n_max": n_max,
        "K2": K2,
        "M_modes": M,
        "wavevector_subsets": subs,
        "all_subset": all(s["subset"] for s in subs if s["n_small"] < n_max),
        "conclusion": (
            "Wavevectors at n=12,16 embed in n=24 dealias shell lists. "
            "Certifying full-dealias at N=24 majorant dominates embedded truncations "
            "n<=24 IF the terminal functional L-0071 is applied on the same n grid "
            "and majorant sup increases with mode-set inclusion. "
            "NOT automatic for different n indexing without embedding proof (N7 sketch)."
        ),
        "evidence_level": "N7_sketch",
        "requires_per_n_cert": False,
    }
    if run_band_ladder:
        out["band_majorant_ladder"] = band_majorant_ladder((1, 2, 3), prec=128)
    return out


def repair_manifest_integral_hi(manifest: dict, *, prec: int = 128) -> mpfr:
    """Route A integral from repair per-shell manifest blocks."""
    return shell_integral_hi_from_blocks(manifest["blocks"], prec=prec)


def load_repair_manifest(n: int = 24) -> dict | None:
    bases = (
        Path("experiments/terminal_weighted"),
        Path("python/experiments/terminal_weighted"),
    )
    for base in bases:
        full = base / f"shell_manifest_repair_full_one_inf_n{n}.json"
        if full.is_file():
            return json.loads(full.read_text(encoding="utf-8"))
        partial = base / f"shell_manifest_repair_full_one_inf_n{n}_partial.json"
        if partial.is_file():
            return json.loads(partial.read_text(encoding="utf-8"))
    return None


def full_dealias_inclusion_audit(*, prec: int = 128) -> dict:
    """
    L-0076 empirical audit: wavevector embedding + band ladders + n=24 repair I_hi.
    Formal N7 proof still open.
    """
    n_max = 24
    dom = n_domination_analysis(n_max, run_band_ladder=True)
    band6 = band_majorant_ladder((1, 2, 3, 4, 5, 6), prec=prec)
    embed_checks = [embed_basis_index_check(n, n_max) for n in (12, 16, 20) if n < n_max]
    manifest = load_repair_manifest(n_max)
    repair = {}
    if manifest is not None:
        I_hi = repair_manifest_integral_hi(manifest, prec=prec)
        repair = {
            "n": manifest.get("n"),
            "n_shells": len(manifest.get("blocks", [])),
            "repair_complete": manifest.get("repair_complete"),
            "D_full": manifest.get("D_full"),
            "I_hi_route_A": str(I_hi),
            "bound_method": manifest.get("bound_method"),
        }
    return {
        "lemma": "L-0076",
        "evidence_level": "N7_sketch_empirical",
        "conjecture_C0008": "exploring",
        "wavevector_domination": dom,
        "embed_basis_checks": embed_checks,
        "all_wavevectors_embedded": all(c["all_wavevectors_embedded"] for c in embed_checks),
        "band_majorant_ladder_3": dom.get("band_majorant_ladder"),
        "band_majorant_ladder_6": band6,
        "repair_full_dealias_n24": repair,
        "formal_gap": [
            "Embedding z_{n'} into z_n preserving terminal functional (N7).",
            "Per-n full-dealias certificates for all n<=24 (only n=24 complete).",
            "Monotonicity of sup f_n under mode-set inclusion not proved.",
        ],
        "conditional_conclusion": (
            "IF mode embedding + majorant monotonicity hold, then n=24 repair "
            "full-dealias I_hi majorizes embedded n'<=24 truncations."
        ),
    }
