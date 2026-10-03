"""Sprint 02: L-0059 Lean N8 SOS cert bridge."""

from __future__ import annotations

from ns_exploration.conjectures.l0059_lean_sos_cert_bridge import lemma_l0059, save_lemma_l0059


def main() -> None:
    l59 = lemma_l0059(run_lean=True)
    save_lemma_l0059(l59)
    print(f"L-0059: status={l59.status} lean_build={l59.lean_build_ok}")
    for p in l59.pins or []:
        print(f"  {p.cert_id} -> {p.lean_theorem}")


if __name__ == "__main__":
    main()
