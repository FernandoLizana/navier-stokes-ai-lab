"""
Sprint B (continuación): cross-check FFT del stretch + tensor Hermitiano físico.

Sprint A verificó que el tensor de l0040 coincide consigo mismo (matricización Sym),
pero ese objeto usa amplitudes reales independientes en cada modo complejo y
OMITE el factor i de ∂ → ik. El stretch físico (campos reales Hermitianos +
pseudospectral) SÍ lleva ese factor.

Este módulo construye el tensor cúbico correcto para campos reales:
  û_k = Σ (c_α - i s_α)/(|k|√2) · p_α ,   û_{-k} = conj(û_k),
  stretch = ⟨ω, curl N⟩ = f(z),   ‖z‖² = 2Ω,
y lo valida contra stretch_inner (FFT) a ~1e-13.

Luego reporta C_sym / C_fullsym / C_rank1_lo sobre ESE objeto físico.
NO modifica lemmas. NO Clay. Finito Galerkin.
"""

from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from ns_exploration.conjectures.l0021_quartic_shell import (
    _leray_vec,
    _pol_basis,
    shell_modes_by_r,
)
from ns_exploration.conjectures.l0023_cubic_dissipation import stretch_inner
from ns_exploration.conjectures.l0027_low_slab_cubic import lemma_l0027
from ns_exploration.diagnostics.metrics import compute_diagnostics
from ns_exploration.experiments.sprintA_stretch_tensor import (
    BAND_123456,
    f_from_G,
    full_sym,
    operator_C,
    rank1_lower_C,
    sym_ij,
)
from ns_exploration.spectral.dealias import dealias_mask
from ns_exploration.spectral.fourier_conventions import wave_numbers_1d


def in_half(k: tuple[int, int, int]) -> bool:
    if k[0] > 0:
        return True
    if k[0] < 0:
        return False
    if k[1] > 0:
        return True
    if k[1] < 0:
        return False
    return k[2] > 0


@dataclass
class HermitianBasis:
    """Real (c,s) coords for Hermitian div-free fields on a shell band."""

    entries: list[tuple[tuple[int, int, int], np.ndarray, float]]
    psis: list[list[tuple[tuple[int, int, int], np.ndarray]]]
    D: int
    n: int
    idx: dict[int, int]


def build_hermitian_basis(n: int, radii: tuple[int, ...]) -> HermitianBasis:
    by = shell_modes_by_r(n)
    kx1 = wave_numbers_1d(n)
    idx = {int(round(v)): i for i, v in enumerate(kx1)}
    entries: list[tuple[tuple[int, int, int], np.ndarray, float]] = []
    for r in radii:
        if r not in by:
            continue
        for k in by[r]:
            if not in_half(k):
                continue
            lam = float(k[0] * k[0] + k[1] * k[1] + k[2] * k[2])
            for p in _pol_basis(np.array(k, float)):
                entries.append((k, np.asarray(p, float), lam))
    psis: list[list[tuple[tuple[int, int, int], np.ndarray]]] = []
    for k, p, lam in entries:
        scale = 1.0 / (math.sqrt(lam) * math.sqrt(2.0))
        mk = (-k[0], -k[1], -k[2])
        # c-coord: û_k += scale*p, û_{-k} += scale*p
        psis.append([(k, scale * p), (mk, scale * p)])
        # s-coord: û_k += -i*scale*p, û_{-k} += i*scale*p
        psis.append([(k, (-1j * scale) * p), (mk, (1j * scale) * p)])
    return HermitianBasis(entries=entries, psis=psis, D=len(psis), n=n, idx=idx)


def build_physical_G(n: int, radii: tuple[int, ...]) -> tuple[np.ndarray, HermitianBasis]:
    """Cubic G with i-factor: stretch_inner(uh(z)) = einsum G z z z."""
    b = build_hermitian_basis(n, radii)
    by = shell_modes_by_r(n)
    out_set: set[tuple[int, int, int]] = set()
    for r in radii:
        if r not in by:
            continue
        out_set.update(by[r])
    D = b.D
    G = np.zeros((D, D, D), dtype=np.float64)
    for a in range(D):
        for bb in range(D):
            Nloc: dict[tuple[int, int, int], np.ndarray] = {}
            for kp, va in b.psis[a]:
                for kq, vb in b.psis[bb]:
                    s = (kp[0] + kq[0], kp[1] + kq[1], kp[2] + kq[2])
                    if s not in out_set:
                        continue
                    coef = 1j * np.dot(va, np.array(kq, float))
                    if s not in Nloc:
                        Nloc[s] = np.zeros(3, dtype=np.complex128)
                    Nloc[s] += coef * vb
            for s, vec in list(Nloc.items()):
                Nloc[s] = -_leray_vec(np.array(s, float), vec)
            for m in range(D):
                acc = 0.0
                for km, vm in b.psis[m]:
                    if km not in Nloc:
                        continue
                    lam = float(km[0] * km[0] + km[1] * km[1] + km[2] * km[2])
                    acc += lam * float(np.real(np.vdot(vm, Nloc[km])))
                G[m, a, bb] = acc
    return G, b


def z_to_uhat(z: np.ndarray, b: HermitianBasis) -> np.ndarray:
    uh = np.zeros((3, b.n, b.n, b.n), dtype=np.complex128)
    for a in range(b.D):
        for k, v in b.psis[a]:
            uh[:, b.idx[k[0]], b.idx[k[1]], b.idx[k[2]]] += z[a] * v
    return uh


@dataclass
class SprintBResult:
    n: int
    radii: list[int]
    D: int
    C_dagger: float
    max_fft_vs_tensor_err: float
    max_omega_rel_err: float
    C_sym_physical: float
    C_fullsym_physical: float
    C_rank1_lower_physical: float
    C_sym_l0040_galerkin: float
    physical_closes_band: bool
    notes: str

    def as_dict(self) -> dict:
        return asdict(self)


def run(
    n: int = 24,
    radii: tuple[int, ...] = BAND_123456,
    n_probe: int = 200,
) -> SprintBResult:
    Cd = lemma_l0027(n=n, empirical=False).C_dagger
    G, basis = build_physical_G(n, radii)
    D = basis.D
    Gsym = sym_ij(G)
    A = full_sym(G)

    rng = np.random.default_rng(20260718)
    err_fft = 0.0
    err_om = 0.0
    for _ in range(n_probe):
        z = rng.standard_normal(D)
        uh = z_to_uhat(z, basis)
        st = stretch_inner(uh)
        ft = f_from_G(G, z)
        scale = max(1.0, abs(st), abs(ft))
        err_fft = max(err_fft, abs(st - ft) / scale)
        om = compute_diagnostics(uh, 0.1).enstrophy
        omz = 0.5 * float(z @ z)
        err_om = max(err_om, abs(om - omz) / max(1.0, om, omz))

    C_sym = operator_C(Gsym)
    C_full = operator_C(A)
    C_r1 = rank1_lower_C(A, n_restart=80, iters=80)
    # reference Galerkin (l0040) value already known
    C_gal = 11.344112520990302
    notes = (
        f"FFT vs physical tensor err={err_fft:.2e}; Ω err={err_om:.2e}; "
        f"C_sym_phys={C_sym:.6f}; C_fullsym_phys={C_full:.6f}; "
        f"C_rank1_lo={C_r1:.6f}; C_†={Cd:.6f}; "
        f"closes={C_full <= Cd}. "
        f"(l0040 Galerkin C_sym={C_gal:.6f} is a different object.)"
    )
    return SprintBResult(
        n=n,
        radii=list(radii),
        D=D,
        C_dagger=Cd,
        max_fft_vs_tensor_err=err_fft,
        max_omega_rel_err=err_om,
        C_sym_physical=C_sym,
        C_fullsym_physical=C_full,
        C_rank1_lower_physical=C_r1,
        C_sym_l0040_galerkin=C_gal,
        physical_closes_band=bool(C_full <= Cd),
        notes=notes,
    )


def main() -> dict:
    res = run()
    out = res.as_dict()
    Path("reports").mkdir(exist_ok=True)
    Path("experiments/exploratory/sprintB_physical_tensor").mkdir(
        parents=True, exist_ok=True
    )
    Path("experiments/exploratory/sprintB_physical_tensor/summary.json").write_text(
        json.dumps(out, indent=2), encoding="utf-8"
    )
    print(json.dumps(out, indent=2, ensure_ascii=True))
    return out


if __name__ == "__main__":
    main()
