"""Sprint: L-0051 structured ladder synthesis (lightweight, no Julia SDP)."""

from __future__ import annotations

import json

from ns_exploration.conjectures.l0051_structured_ladder import lemma_l0051, save_lemma_l0051


def main() -> None:
    b = lemma_l0051()
    save_lemma_l0051(b)
    out = {
        "lemma_id": b.lemma_id,
        "C_dagger": b.C_dagger,
        "ratio_median_large": b.ratio_median_large,
        "C_ub_sos_extrapolated_full": b.C_ub_sos_extrapolated_full,
        "brute_sos_closes_all_ic": b.brute_sos_closes_all_ic,
        "n_scale_rows": len(b.sos_scale_rows or []),
        "planned_next": "refresh-closure-map-l0060",
    }
    print(json.dumps(out, indent=2), flush=True)
    if b.structured_routes:
        print("\nStructured routes (ranked):", flush=True)
        for r in b.structured_routes:
            print(
                f"  [{r.status}] {r.id}: {r.target} - {r.resource}",
                flush=True,
            )


if __name__ == "__main__":
    main()
