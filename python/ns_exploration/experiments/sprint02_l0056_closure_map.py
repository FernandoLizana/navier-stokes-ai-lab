"""Sprint 02: L-0056 C-0007 structured closure map."""

from __future__ import annotations

from ns_exploration.conjectures.l0056_closure_map import (
    lemma_l0056,
    save_lemma_l0056,
    write_closure_map_report,
)


def main() -> None:
    l56 = lemma_l0056()
    save_lemma_l0056(l56)
    report = write_closure_map_report(l56)
    print(f"L-0056: {l56.lemma_id} status={l56.status}")
    print(f"  proved={l56.n_proved_subclasses} open={l56.n_open_subclasses}")
    print(f"  best_n5={l56.best_n5_sos_hi} structured={l56.best_structured_sos_hi}")
    print(f"  report={report}")


if __name__ == "__main__":
    main()
