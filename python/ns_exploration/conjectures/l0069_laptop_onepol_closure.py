"""
L-0069: Laptop one-pol path closure synthesis (L-0064…L-0068).

Documents techo at shell 25 for greedy-second, marks cluster greedy SOS
deferred, and packages final structured subclass ladder through C-R-0014.

FINITE Galerkin only. Not all-IC. Not Clay.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path

from ns_exploration.conjectures.l0026_c0007_techo import lemma_l0026
from ns_exploration.conjectures.l0027_low_slab_cubic import lemma_l0027
from ns_exploration.conjectures.l0043_onepol_shellblock import sym_onepol_C_shor
from ns_exploration.conjectures.l0067_onepol_greedy_r24_sos import SUPPORT_CR0014, POL_BRANCH
from ns_exploration.conjectures.l0065_onepol_greedy_second_sos import SUPPORT_GREEDY_SECOND

SUPPORT_TECHO_25 = tuple(sorted(SUPPORT_CR0014 + [25]))


@dataclass
class LaptopClosureRow:
    id: str
    shells: list[int]
    D: int
    C_shor: float
    C_ub_sos_hi: float | None
    status: str


@dataclass
class GalerkinBoundL0069:
    lemma_id: str = "L-0069"
    route: str = "B"
    status: str = "validated"
    evidence_level: str = "N7"
    n: int = 24
    C_dagger: float = 0.0
    n_proved_restricted: int = 0
    best_onepol_cr: str = "C-R-0014"
    best_onepol_D: int = 178
    best_onepol_C_shor: float = 0.0
    best_onepol_C_ub_sos_hi: float | None = None
    techo_shell: int = 25
    techo_C_shor: float = 0.0
    cluster_greedy_sos: str = "deferred"
    all_ic_open: bool = False
    c0007_proved_M: float = 0.0
    c0008_sharp_open: bool = True
    ladder: list[LaptopClosureRow] | None = None
    clay_implication: str = "None. Laptop structured path closure map; not Clay."
    notes: str = ""

    def as_dict(self) -> dict:
        d = asdict(self)
        if self.ladder:
            d["ladder"] = [asdict(r) for r in self.ladder]
        return d


def lemma_l0069(n: int = 24) -> GalerkinBoundL0069:
    Cd = lemma_l0027(n=n, empirical=False).C_dagger
    cr_dir = Path("conjectures/proved_restricted")
    n_cr = len(list(cr_dir.glob("C-R-*.json")))

    l65 = json.loads(Path("conjectures/active/L-0065.json").read_text(encoding="utf-8"))
    l67 = json.loads(Path("conjectures/active/L-0067.json").read_text(encoding="utf-8"))
    cr14 = json.loads((cr_dir / "C-R-0014.json").read_text(encoding="utf-8"))

    techo = sym_onepol_C_shor(n, SUPPORT_TECHO_25, branch=POL_BRANCH)
    assert float(techo["C_shor_sym"]) > Cd + 1e-6

    ladder = [
        LaptopClosureRow(
            id="C-R-0008",
            shells=[1, 2, 3, 4, 5, 6, 8],
            D=92,
            C_shor=8.4764,
            C_ub_sos_hi=3.356,
            status="proved_N5",
        ),
        LaptopClosureRow(
            id="C-R-0013",
            shells=list(SUPPORT_GREEDY_SECOND),
            D=154,
            C_shor=float(l65["C_shor_sym"]),
            C_ub_sos_hi=float(l65["C_ub_hi"]),
            status="proved_N5",
        ),
        LaptopClosureRow(
            id="C-R-0014",
            shells=list(SUPPORT_CR0014),
            D=int(cr14["D"]),
            C_shor=float(cr14["C_shor_sym"]),
            C_ub_sos_hi=float(l67["C_ub_hi"]),
            status="proved_N5",
        ),
    ]

    b26 = lemma_l0026(n=n)
    M7 = b26.L0024_majorant
    notes = (
        f"Laptop one-pol path closed at C-R-0014 ({len(SUPPORT_CR0014)} shells, D={cr14['D']}). "
        f"Adding shell 25 gives C_Shor≈{float(techo['C_shor_sym']):.4f}>C_†≈{Cd:.4f} (techo). "
        f"{n_cr} restricted subclasses proved. Greedy cluster SOS (L-0057) deferred. "
        f"C-0007 all-IC closed at M≈{M7:.4f} (L-0070); C-0008 sharp still open."
    )
    return GalerkinBoundL0069(
        n=n,
        C_dagger=Cd,
        n_proved_restricted=n_cr,
        best_onepol_C_shor=float(cr14["C_shor_sym"]),
        best_onepol_C_ub_sos_hi=float(l67["C_ub_hi"]),
        techo_shell=25,
        techo_C_shor=float(techo["C_shor_sym"]),
        all_ic_open=False,
        c0007_proved_M=float(M7),
        c0008_sharp_open=True,
        ladder=ladder,
        notes=notes,
    )


def save_lemma_l0069(
    bound: GalerkinBoundL0069,
    path: str | Path = "conjectures/active/L-0069.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(bound.as_dict(), indent=2), encoding="utf-8")


def patch_l0064_status() -> None:
    path = Path("conjectures/active/L-0064.json")
    if not path.is_file():
        return
    d = json.loads(path.read_text(encoding="utf-8"))
    d["status"] = "validated"
    d["notes"] = (
        "One-pol scan complete. Greedy-second closes through shell 24 (C-R-0014). "
        "Shell 25 techo. Cluster greedy SOS deferred."
    )
    path.write_text(json.dumps(d, indent=2), encoding="utf-8")
