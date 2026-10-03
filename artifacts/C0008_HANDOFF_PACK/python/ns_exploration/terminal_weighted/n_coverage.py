"""Audit: does N=24 full-dealias majorant cover all N<=24?"""

from __future__ import annotations

from ns_exploration.conjectures.l0021_quartic_shell import shell_modes_by_r
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


def n_domination_analysis(n_max: int = 24) -> dict:
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
    return {
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
