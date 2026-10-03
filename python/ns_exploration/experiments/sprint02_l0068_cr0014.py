"""Sprint: L-0067/68 close C-R-0014 (12-shell one-pol + shell 24)."""

from __future__ import annotations

import json
from pathlib import Path

from ns_exploration.conjectures.l0056_closure_map import (
    lemma_l0056,
    save_lemma_l0056,
    write_closure_map_report,
)
from ns_exploration.conjectures.l0067_onepol_greedy_r24_sos import lemma_l0067, save_lemma_l0067
from ns_exploration.conjectures.l0068_cr0014_onepol_greedy import (
    lemma_l0068,
    save_cr0014,
    save_lemma_l0068,
)
from ns_exploration.validation.l0067_certificate import (
    build_l0067_certificate,
    save_l0067_certificate,
    verify_l0067_certificate,
)
from ns_exploration.validation.l0068_certificate import (
    build_l0068_certificate,
    save_l0068_certificate,
    verify_l0068_certificate,
)


def main() -> dict:
    l67 = lemma_l0067()
    save_lemma_l0067(l67)
    c67 = build_l0067_certificate()
    ok67, chk67 = verify_l0067_certificate(c67.as_dict())
    save_l0067_certificate(c67)
    if not ok67:
        raise RuntimeError(chk67)

    l68 = lemma_l0068()
    save_lemma_l0068(l68)
    save_cr0014(l68)
    c68 = build_l0068_certificate()
    ok68, chk68 = verify_l0068_certificate(c68.as_dict())
    save_l0068_certificate(c68)
    if not ok68:
        raise RuntimeError(chk68)

    l56 = lemma_l0056()
    save_lemma_l0056(l56)
    write_closure_map_report(l56)

    c7 = json.loads(Path("conjectures/active/C-0007.json").read_text(encoding="utf-8"))
    rel = set(c7.get("related") or [])
    rel.update(["L-0067", "L-0068", "C-R-0014"])
    c7["related"] = sorted(rel)
    c7["notes"] = (
        f"Open all-IC. Closed C-R-0002..0014 ({l56.n_proved_restricted}). "
        f"C-R-0014: 12 shells D={l68.support_cr0014_D} C_Shor≈{l68.C_shor_sym_cr0014:.4f} "
        f"SOS≈{l68.C_ub_sos_hi:.4f}."
    )
    Path("conjectures/active/C-0007.json").write_text(json.dumps(c7, indent=2), encoding="utf-8")

    out = {
        "C-R-0014": True,
        "D": l68.support_cr0014_D,
        "C_shor": l68.C_shor_sym_cr0014,
        "C_ub_sos_hi": l68.C_ub_sos_hi,
        "n_proved": l56.n_proved_restricted,
        "ok67": ok67,
        "ok68": ok68,
    }
    print(json.dumps(out, indent=2))
    return out


if __name__ == "__main__":
    main()
