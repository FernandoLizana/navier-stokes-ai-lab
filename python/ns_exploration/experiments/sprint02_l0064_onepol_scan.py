"""Sprint: L-0064 one-pol Shor scanner (laptop-friendly structured search)."""

from __future__ import annotations

import json
from pathlib import Path

from ns_exploration.conjectures.l0064_onepol_scanner import (
    lemma_l0064,
    save_lemma_l0064,
    save_report,
)


def main() -> dict:
    bound = lemma_l0064(n=24)
    save_lemma_l0064(bound)
    save_report(bound)

    c7_path = Path("conjectures/active/C-0007.json")
    if c7_path.is_file():
        c7 = json.loads(c7_path.read_text(encoding="utf-8"))
        related = set(c7.get("related") or [])
        related.add("L-0064")
        c7["related"] = sorted(related)
        c7["notes"] = (
            f"Open all-IC. L-0064 one-pol scan: {bound.n_closing}/{bound.n_candidates} "
            f"close C_†; greedy first {len(bound.greedy_growth_first or [])} shells; "
            f"SOS queue {len(bound.sos_queue or [])}. Cluster greedy SOS deferred."
        )
        c7_path.write_text(json.dumps(c7, indent=2), encoding="utf-8")

    out = {
        "lemma_id": bound.lemma_id,
        "n_candidates": bound.n_candidates,
        "n_closing": bound.n_closing,
        "best_D": bound.best_closing.D if bound.best_closing else None,
        "best_C": bound.best_closing.C_shor_sym if bound.best_closing else None,
        "greedy_first_n": len(bound.greedy_growth_first or []),
        "greedy_second_n": len(bound.greedy_growth_second or []),
        "sos_queue_len": len(bound.sos_queue or []),
        "not_a_clay_claim": True,
    }
    Path("experiments/exploratory/sprint02_l0064").mkdir(parents=True, exist_ok=True)
    Path("experiments/exploratory/sprint02_l0064/summary.json").write_text(
        json.dumps(out, indent=2), encoding="utf-8"
    )
    print(json.dumps(out, indent=2, ensure_ascii=True))
    if bound.sos_queue:
        print("\nSOS laptop queue (top 5):", flush=True)
        for r in bound.sos_queue[:5]:
            print(
                f"  {r.scan_id}: D={r.D} C={r.C_shor_sym:.4f} radii={r.radii}",
                flush=True,
            )
    return out


if __name__ == "__main__":
    main()
