"""Sprint: L-0054 structured SOS synthesis."""

from __future__ import annotations

import json

from ns_exploration.conjectures.l0054_structured_sos_synthesis import lemma_l0054, save_lemma_l0054


def main() -> None:
    b = lemma_l0054()
    save_lemma_l0054(b)
    out = {
        "lemma_id": b.lemma_id,
        "C_dagger": b.C_dagger,
        "best_sos_n5_hi": b.best_sos_n5_hi,
        "best_sos_subclass": b.best_sos_subclass,
        "n_rows": len(b.rows or []),
        "all_ic_open": b.all_ic_open,
    }
    print(json.dumps(out, indent=2), flush=True)


if __name__ == "__main__":
    main()
