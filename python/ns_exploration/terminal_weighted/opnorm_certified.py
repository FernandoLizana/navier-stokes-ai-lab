"""
Certified operator-norm upper bounds via Frobenius accumulation (MPFR upward).
"""

from __future__ import annotations

import hashlib
import json
import math
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
from gmpy2 import mpfr, sqrt

from ns_exploration.conjectures.l0021_quartic_shell import _leray_vec, shell_modes_by_r
from ns_exploration.conjectures.l0040_sym_shor import _sym_col_index
from ns_exploration.conjectures.l0045_streaming_fullsym import _mode_to_basis
from ns_exploration.conjectures.l0048_matrixfree_fullsym import all_dealias_radii
from ns_exploration.experiments.sprintB_physical_tensor import build_hermitian_basis
from ns_exploration.terminal_weighted.constants import (
    DEFAULT_BOUND_METHOD,
    FROBENIUS_STREAMING_ENABLED,
    FROZEN,
)
from ns_exploration.terminal_weighted.intervals import (
    ArithmeticBackend,
    add_up,
    exp_up,
    mpfr_const,
    mul_up,
    sub_up,
    sqrt_up,
    _ctx,
)
from ns_exploration.terminal_weighted.rational_basis import mpfr_leray, mpfr_vdot_real
from ns_exploration.terminal_weighted.opnorm_one_inf import (
    abs_up,
    merge_sum_dicts,
    one_inf_opnorm_hi,
    one_inf_opnorm_hi_rowlist,
)
from ns_exploration.terminal_weighted.tensor import _inner_psi_N_weighted

_CHECKPOINT_PATH = Path("experiments/terminal_weighted/frobenius_checkpoint.json")
_SHELL_ONE_INF_CK_DIR = Path("experiments/terminal_weighted/upgrade_shell_ck")
_MPFR_BATCH = 64
_SHELL_CK_EVERY = 500_000


@dataclass
class CertifiedOpNormBound:
    shell: int | None
    n: int
    D: int
    frobenius_sq_hi: str
    L_op_hi: str
    C_term_hi: str
    n_updates: int
    method: str
    evidence_level: str
    payload_sha256: str

    def as_dict(self) -> dict:
        return asdict(self)


def _accumulate_frobenius_sq_hi(
    n: int,
    radii: tuple[int, ...],
    *,
    output_shell: int | None = None,
    prec: int = 200,
) -> tuple[mpfr, int, int]:
    fro_by, n_by, fro_full, n_full, D, *_ = _accumulate_frobenius_all_shells_single_pass(
        n, radii, prec=prec, workers=1
    )
    if output_shell is None:
        return fro_full, n_full, D
    return fro_by[output_shell], n_by[output_shell], D


def _bound_from_one_inf(
    row_sum: dict[int, mpfr],
    col_sum: dict[int, mpfr],
    n_updates: int,
    n: int,
    D: int,
    *,
    output_shell: int | None,
    prec: int,
) -> CertifiedOpNormBound:
    L_hi = one_inf_opnorm_hi(row_sum, col_sum)
    two_sqrt2 = sqrt_up(mpfr(8, precision=prec))
    C_hi = mul_up(two_sqrt2, L_hi)
    payload = {
        "n": n,
        "shell": output_shell,
        "norm_1_inf": str(L_hi),
        "n_updates": n_updates,
    }
    digest = hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()
    return CertifiedOpNormBound(
        shell=output_shell,
        n=n,
        D=D,
        frobenius_sq_hi="",
        L_op_hi=str(L_hi),
        C_term_hi=str(C_hi),
        n_updates=n_updates,
        method="one_inf_mpfr_up",
        evidence_level="N4",
        payload_sha256=digest,
    )


def _pick_tighter_bound(
    fro: CertifiedOpNormBound, one: CertifiedOpNormBound, *, prec: int
) -> CertifiedOpNormBound:
    if mpfr(one.C_term_hi) < mpfr(fro.C_term_hi):
        return one
    return fro


def _shell_of_mode(b, m: int) -> int:
    for km, _ in b.psis[m]:
        return int(km[0] * km[0] + km[1] * km[1] + km[2] * km[2])
    return 0


def _apply_leray_nloc(
    Nloc: dict[tuple[int, int, int], np.ndarray],
    *,
    rational: bool,
    prec: int,
) -> None:
    for s, vec in list(Nloc.items()):
        if rational:
            Nloc[s] = -mpfr_leray(s, vec, prec=prec)
        else:
            Nloc[s] = -_leray_vec(np.array(s, float), vec)


def _inner_psi_N_weighted_certified(
    b,
    m: int,
    Nloc: dict[tuple[int, int, int], np.ndarray],
    t: float,
    nu: float,
    T: float,
    *,
    output_shell: int | None = None,
    output_shells: frozenset[int] | None = None,
    rational: bool = False,
    prec: int = 128,
) -> float:
    if not rational:
        return _inner_psi_N_weighted(
            b,
            m,
            Nloc,
            t,
            nu,
            T,
            output_shell=output_shell,
            output_shells=output_shells,
        )
    _ctx(prec, up=True)
    sh = _shell_of_mode(b, m)
    if output_shell is not None and sh != output_shell:
        return 0.0
    if output_shells is not None and sh not in output_shells:
        return 0.0
    acc = mpfr(0, precision=prec)
    nu_m = mpfr_const(nu, prec=prec)
    T_m = mpfr_const(T, prec=prec)
    t_m = mpfr(t, precision=prec)
    neg_two = mpfr(-2, precision=prec)
    for km, vm in b.psis[m]:
        if km not in Nloc:
            continue
        r2 = km[0] * km[0] + km[1] * km[1] + km[2] * km[2]
        rf = mpfr(r2, precision=prec)
        expo = mul_up(mul_up(mul_up(neg_two, nu_m), rf), sub_up(T_m, t_m))
        w = mul_up(rf, exp_up(expo))
        re = mpfr_vdot_real(vm, Nloc[km], prec=prec, up=True)
        acc = add_up(acc, mul_up(w, re))
    return float(acc)


def _delta_sq_hi(delta: float, prec: int) -> mpfr:
    with _ctx(prec, up=True):
        d = mpfr(delta)
        return mpfr(d * d)


def _sym_col(p: int, q: int, r: int, D: int) -> int:
    if q == r:
        return q
    lo, hi = (q, r) if q < r else (r, q)
    return _sym_col_index(lo, hi, D)


def _merge_entry_maps(
    dst: dict[tuple[int, int], mpfr],
    src: dict[tuple[int, int], mpfr],
    *,
    prec: int,
) -> None:
    z = mpfr(0, precision=prec)
    for key, val in src.items():
        prev = dst.get(key, z)
        with _ctx(prec):
            dst[key] = mpfr(prev + val)


def _finalize_frobenius_from_entries(
    entry_by: dict[int, dict[tuple[int, int], mpfr]],
    shells: list[int],
    *,
    prec: int,
) -> tuple[dict[int, mpfr], mpfr, dict[int, int]]:
    """Sum |M_ij|^2 from signed entry accumulators (repair P0 #1)."""
    fro_by = {r: mpfr(0, precision=prec) for r in shells}
    n_by = {r: 0 for r in shells}
    fro_full = mpfr(0, precision=prec)
    for sr in shells:
        emap = entry_by.get(sr, {})
        for val in emap.values():
            ad = abs_up(val, prec)
            sq = mul_up(ad, ad)
            fro_by[sr] = add_up(fro_by[sr], sq)
            fro_full = add_up(fro_full, sq)
            n_by[sr] += 1
    return fro_by, fro_full, n_by


def _flush_entry_batch(
    batch: list[tuple[int, int, int, float, int | None]],
    n_full: int,
    row_by: dict[int, dict[int, mpfr]],
    col_by: dict[int, dict[int, mpfr]],
    entry_by: dict[int, dict[tuple[int, int], mpfr]],
    *,
    prec: int,
    D: int,
) -> int:
    z = mpfr(0, precision=prec)
    for p, q, r, delta, sr in batch:
        n_full += 1
        if sr is not None:
            ad = abs_up(delta, prec)
            col = _sym_col(p, q, r, D)
            rs = row_by.setdefault(sr, {})
            cs = col_by.setdefault(sr, {})
            rs[p] = add_up(rs.get(p, z), ad)
            cs[col] = add_up(cs.get(col, z), ad)
            em = entry_by.setdefault(sr, {})
            key = (p, col)
            with _ctx(prec):
                em[key] = mpfr(em.get(key, z) + mpfr(delta))
    batch.clear()
    return n_full


def _flush_entry_batch_one_shell(
    batch: list[tuple[int, int, int, float, int | None]],
    target_shell: int,
    row_sum: list[mpfr],
    col_sum: dict[int, mpfr],
    n_updates: int,
    *,
    prec: int,
    D: int,
) -> int:
    z = mpfr(0, precision=prec)
    for p, q, r, delta, sr in batch:
        if sr != target_shell:
            continue
        n_updates += 1
        ad = abs_up(delta, prec)
        if q == r:
            col = q
        else:
            lo, hi = (q, r) if q < r else (r, q)
            col = _sym_col_index(lo, hi, D)
        row_sum[p] = add_up(row_sum[p], ad)
        prev = col_sum.get(col, z)
        col_sum[col] = add_up(prev, ad)
    batch.clear()
    return n_updates


def _process_pair_range_one_shell(
    b,
    out_set: set,
    hit: dict,
    mode_shell: list[int],
    target_shell: int,
    row_sum: list[mpfr],
    col_sum: dict[int, mpfr],
    n_updates: int,
    *,
    a_start: int,
    a_end: int,
    prec: int,
    progress_cb=None,
    start_pair: int = 0,
    n_pairs: int = 0,
    checkpoint_cb=None,
    rational_coefficients: bool = False,
) -> tuple[int, int]:
    D = b.D
    sq2 = math.sqrt(2.0)
    t, nu, T = FROZEN.T, FROZEN.nu, FROZEN.T
    done = start_pair
    a0 = start_pair // D
    bb0 = start_pair % D
    batch: list[tuple[int, int, int, float, int | None]] = []

    for a in range(max(a_start, a0), a_end):
        bb_rng = range(bb0, D) if a == a0 else range(D)
        for bb in bb_rng:
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
            done += 1
            if not Nloc:
                if progress_cb and n_pairs and done % 250_000 == 0:
                    progress_cb(done, n_pairs)
                continue
            _apply_leray_nloc(Nloc, rational=rational_coefficients, prec=prec)
            ms: set[int] = set()
            for s in Nloc:
                ms.update(hit.get(s, ()))
            for m in ms:
                g = _inner_psi_N_weighted_certified(
                    b, m, Nloc, t, nu, T, output_shell=target_shell,
                    rational=rational_coefficients, prec=prec,
                )
                if g == 0.0:
                    continue
                contrib = g / 6.0
                sr = target_shell
                for p, q, r in (
                    (m, a, bb),
                    (m, bb, a),
                    (a, m, bb),
                    (a, bb, m),
                    (bb, m, a),
                    (bb, a, m),
                ):
                    delta = 0.5 * contrib if q == r else 0.5 * sq2 * contrib
                    batch.append((p, q, r, delta, sr))
                    if len(batch) >= _MPFR_BATCH:
                        n_updates = _flush_entry_batch_one_shell(
                            batch,
                            target_shell,
                            row_sum,
                            col_sum,
                            n_updates,
                            prec=prec,
                            D=D,
                        )
            if progress_cb and n_pairs and done % 250_000 == 0:
                progress_cb(done, n_pairs)
            if checkpoint_cb and n_pairs and done % _SHELL_CK_EVERY == 0:
                checkpoint_cb(done, n_updates, row_sum, col_sum)
    if batch:
        n_updates = _flush_entry_batch_one_shell(
            batch,
            target_shell,
            row_sum,
            col_sum,
            n_updates,
            prec=prec,
            D=D,
        )
    return n_updates, done


def _atomic_json_write(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    tmp.replace(path)


def _save_one_shell_ck(
    path: Path,
    *,
    shell: int,
    n: int,
    radii: tuple[int, ...],
    prec: int,
    done_pairs: int,
    n_pairs: int,
    n_updates: int,
    row_sum: list[mpfr],
    col_sum: dict[int, mpfr],
    rational_coefficients: bool = False,
) -> None:
    _atomic_json_write(
        path,
        {
            "shell": shell,
            "n": n,
            "radii": list(radii),
            "prec": prec,
            "rational_coefficients": rational_coefficients,
            "done_pairs": done_pairs,
            "n_pairs": n_pairs,
            "n_updates": n_updates,
            "row_sum": [str(v) for v in row_sum],
            "col_sum": {str(k): str(v) for k, v in col_sum.items()},
        },
    )


def _load_one_shell_ck(path: Path) -> dict | None:
    if not path.is_file():
        return None
    return json.loads(path.read_text(encoding="utf-8"))


def certified_L_op_one_inf_single_shell_hi(
    n: int,
    radii: tuple[int, ...],
    output_shell: int,
    *,
    prec: int = 200,
    checkpoint_path: Path | None = None,
    heartbeat_path: Path | None = None,
    rational_coefficients: bool = False,
) -> CertifiedOpNormBound:
    """Full D^2 pass tracking 1-inf accumulators for one output shell only (bounded RAM)."""
    _ctx(prec)
    checkpoint_path = checkpoint_path or (
        _SHELL_ONE_INF_CK_DIR / f"shell_{output_shell:03d}.json"
    )
    b = build_hermitian_basis(n, radii)
    by = shell_modes_by_r(n)
    out_set: set[tuple[int, int, int]] = set()
    for r in radii:
        if r in by:
            out_set.update(by[r])
    hit = _mode_to_basis(b)
    D = b.D
    n_pairs = D * D
    start_pair = 0
    n_updates = 0
    row_sum = [mpfr(0, precision=prec) for _ in range(D)]
    col_sum: dict[int, mpfr] = {}

    ck = _load_one_shell_ck(checkpoint_path)
    if (
        ck
        and ck.get("shell") == output_shell
        and ck.get("n") == n
        and tuple(ck.get("radii", [])) == radii
        and ck.get("prec") == prec
        and ck.get("n_pairs") == n_pairs
        and ck.get("rational_coefficients", False) == rational_coefficients
    ):
        start_pair = int(ck["done_pairs"])
        n_updates = int(ck["n_updates"])
        row_sum = [mpfr(v, precision=prec) for v in ck["row_sum"]]
        col_sum = {int(k): mpfr(v, precision=prec) for k, v in ck["col_sum"].items()}
        print(
            f"    shell {output_shell}: resume at pair {start_pair}/{n_pairs}",
            flush=True,
        )

    def _hb(msg: str) -> None:
        if heartbeat_path is None:
            return
        heartbeat_path.parent.mkdir(parents=True, exist_ok=True)
        import time as _time

        heartbeat_path.write_text(
            f"{_time.strftime('%Y-%m-%dT%H:%M:%SZ', _time.gmtime())} shell {output_shell} {msg}\n",
            encoding="utf-8",
        )

    def progress_cb(done: int, total: int) -> None:
        print(f"    shell {output_shell}: pairs {done}/{total}", flush=True)
        if done % 500_000 == 0:
            _hb(f"pairs {done}/{total}")

    def checkpoint_cb(done: int, n_up: int, rows: list[mpfr], cols: dict[int, mpfr]) -> None:
        _save_one_shell_ck(
            checkpoint_path,
            shell=output_shell,
            n=n,
            radii=radii,
            prec=prec,
            done_pairs=done,
            n_pairs=n_pairs,
            n_updates=n_up,
            row_sum=rows,
            col_sum=cols,
            rational_coefficients=rational_coefficients,
        )
        _hb(f"ckpt pairs {done}/{n_pairs}")

    n_updates, _ = _process_pair_range_one_shell(
        b,
        out_set,
        hit,
        [_shell_of_mode(b, m) for m in range(D)],
        output_shell,
        row_sum,
        col_sum,
        n_updates,
        a_start=0,
        a_end=D,
        prec=prec,
        progress_cb=progress_cb,
        start_pair=start_pair,
        n_pairs=n_pairs,
        checkpoint_cb=checkpoint_cb,
        rational_coefficients=rational_coefficients,
    )
    if checkpoint_path.is_file():
        checkpoint_path.unlink()
    return _bound_from_one_inf_rowlist(
        row_sum, col_sum, n_updates, n, D, output_shell=output_shell, prec=prec
    )


def _bound_from_one_inf_rowlist(
    row_sum: list[mpfr],
    col_sum: dict[int, mpfr],
    n_updates: int,
    n: int,
    D: int,
    *,
    output_shell: int | None,
    prec: int,
) -> CertifiedOpNormBound:
    L_hi = one_inf_opnorm_hi_rowlist(row_sum, col_sum)
    two_sqrt2 = sqrt_up(mpfr(8, precision=prec))
    C_hi = mul_up(two_sqrt2, L_hi)
    payload = {
        "n": n,
        "shell": output_shell,
        "norm_1_inf": str(L_hi),
        "n_updates": n_updates,
    }
    digest = hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()
    return CertifiedOpNormBound(
        shell=output_shell,
        n=n,
        D=D,
        frobenius_sq_hi="",
        L_op_hi=str(L_hi),
        C_term_hi=str(C_hi),
        n_updates=n_updates,
        method="one_inf_mpfr_up",
        evidence_level="N4",
        payload_sha256=digest,
    )


def _flush_delta_batch(
    batch: list[float],
    fro_full: mpfr,
    fro_by: dict[int, mpfr],
    n_by: dict[int, int],
    n_full: int,
    shell_tags: list[int | None],
    *,
    prec: int,
) -> tuple[mpfr, dict[int, mpfr], int]:
    for delta, sr in zip(batch, shell_tags):
        val = _delta_sq_hi(delta, prec)
        fro_full = add_up(fro_full, val)
        n_full += 1
        if sr is not None and sr in fro_by:
            fro_by[sr] = add_up(fro_by[sr], val)
            n_by[sr] += 1
    batch.clear()
    shell_tags.clear()
    return fro_full, fro_by, n_full


def _process_pair_range(
    b,
    out_set: set,
    hit: dict,
    mode_shell: list[int],
    shells: list[int],
    n_full: int,
    row_by: dict[int, dict[int, mpfr]],
    col_by: dict[int, dict[int, mpfr]],
    entry_by: dict[int, dict[tuple[int, int], mpfr]],
    *,
    a_start: int,
    a_end: int,
    prec: int,
    progress_cb=None,
    pair_offset: int = 0,
    n_pairs: int = 0,
    rational_coefficients: bool = False,
) -> tuple[dict[int, mpfr], dict[int, int], mpfr, int, int]:
    D = b.D
    sq2 = math.sqrt(2.0)
    t, nu, T = FROZEN.T, FROZEN.nu, FROZEN.T
    done = pair_offset
    batch: list[tuple[int, int, int, float, int | None]] = []

    for a in range(a_start, a_end):
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
            done += 1
            if not Nloc:
                if progress_cb and n_pairs and done % 250_000 == 0:
                    progress_cb(done, n_pairs)
                continue
            _apply_leray_nloc(Nloc, rational=rational_coefficients, prec=prec)
            ms: set[int] = set()
            for s in Nloc:
                ms.update(hit.get(s, ()))
            for m in ms:
                g = _inner_psi_N_weighted_certified(
                    b, m, Nloc, t, nu, T, output_shell=None,
                    rational=rational_coefficients, prec=prec,
                )
                if g == 0.0:
                    continue
                contrib = g / 6.0
                sr = mode_shell[m]
                sr_tag = sr if sr in shells else None
                for p, q, r in (
                    (m, a, bb),
                    (m, bb, a),
                    (a, m, bb),
                    (a, bb, m),
                    (bb, m, a),
                    (bb, a, m),
                ):
                    delta = 0.5 * contrib if q == r else 0.5 * sq2 * contrib
                    batch.append((p, q, r, delta, sr_tag))
                    if len(batch) >= _MPFR_BATCH:
                        n_full = _flush_entry_batch(
                            batch,
                            n_full,
                            row_by,
                            col_by,
                            entry_by,
                            prec=prec,
                            D=D,
                        )
            if progress_cb and n_pairs and done % 250_000 == 0:
                progress_cb(done, n_pairs)
    if batch:
        n_full = _flush_entry_batch(
            batch,
            n_full,
            row_by,
            col_by,
            entry_by,
            prec=prec,
            D=D,
        )
    fro_by, fro_full, n_by = _finalize_frobenius_from_entries(entry_by, shells, prec=prec)
    return fro_by, n_by, fro_full, n_full, done


def _serialize_nested(d: dict[int, dict[int, mpfr]]) -> dict:
    return {str(sh): {str(k): str(v) for k, v in inner.items()} for sh, inner in d.items()}


def _serialize_entry_map(d: dict[int, dict[tuple[int, int], mpfr]]) -> dict:
    return {
        str(sh): {f"{p},{c}": str(v) for (p, c), v in inner.items()} for sh, inner in d.items()
    }


def _fro_worker_chunk(args: tuple) -> dict:
    if len(args) >= 7:
        a_start, a_end, n, radii, prec, shells, rational = args[:7]
    else:
        a_start, a_end, n, radii, prec, shells = args
        rational = False
    _ctx(prec)
    b = build_hermitian_basis(n, radii)
    by = shell_modes_by_r(n)
    out_set: set[tuple[int, int, int]] = set()
    for r in radii:
        if r in by:
            out_set.update(by[r])
    hit = _mode_to_basis(b)
    mode_shell = [_shell_of_mode(b, m) for m in range(b.D)]
    row_by: dict[int, dict[int, mpfr]] = {r: {} for r in shells}
    col_by: dict[int, dict[int, mpfr]] = {r: {} for r in shells}
    entry_by: dict[int, dict[tuple[int, int], mpfr]] = {r: {} for r in shells}
    _, _, _, n_full, _ = _process_pair_range(
        b,
        out_set,
        hit,
        mode_shell,
        shells,
        0,
        row_by,
        col_by,
        entry_by,
        a_start=a_start,
        a_end=a_end,
        prec=prec,
        rational_coefficients=rational,
    )
    return {
        "a_start": a_start,
        "a_end": a_end,
        "n_full": n_full,
        "row_by": _serialize_nested(row_by),
        "col_by": _serialize_nested(col_by),
        "entry_by": _serialize_entry_map(entry_by),
        "D": b.D,
    }


def _save_checkpoint(
    path: Path,
    *,
    done: int,
    n_pairs: int,
    fro_by: dict[int, mpfr],
    n_by: dict[int, int],
    fro_full: mpfr,
    n_full: int,
    prec: int,
    n: int,
    radii: tuple[int, ...],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "done_pairs": done,
        "n_pairs": n_pairs,
        "n": n,
        "radii": list(radii),
        "prec": prec,
        "fro_by": {str(k): str(v) for k, v in fro_by.items()},
        "n_by": {str(k): v for k, v in n_by.items()},
        "fro_full": str(fro_full),
        "n_full": n_full,
    }
    path.write_text(json.dumps(payload), encoding="utf-8")


def _load_checkpoint(path: Path) -> dict | None:
    if not path.is_file():
        return None
    return json.loads(path.read_text(encoding="utf-8"))


def _accumulate_frobenius_all_shells_single_pass(
    n: int,
    radii: tuple[int, ...],
    *,
    prec: int = 200,
    progress_every: int = 250_000,
    workers: int = 1,
    checkpoint_path: Path | None = None,
    rational_coefficients: bool = False,
) -> tuple[dict[int, mpfr], dict[int, int], mpfr, int, int]:
    checkpoint_path = checkpoint_path or _CHECKPOINT_PATH
    shells = sorted(set(radii))

    if workers > 1:
        print(f"    parallel workers={workers}", flush=True)
        b_tmp = build_hermitian_basis(n, radii)
        D = b_tmp.D
        del b_tmp
        chunk_size = max(1, (D + workers - 1) // workers)
        chunks = []
        for i in range(workers):
            a0 = i * chunk_size
            if a0 >= D:
                break
            a1 = min(D, a0 + chunk_size)
            chunks.append((a0, a1, n, radii, prec, shells, rational_coefficients))
        fro_by = {r: mpfr(0, precision=prec) for r in shells}
        n_by = {r: 0 for r in shells}
        fro_full = mpfr(0, precision=prec)
        n_full = 0
        row_by: dict[int, dict[int, mpfr]] = {r: {} for r in shells}
        col_by: dict[int, dict[int, mpfr]] = {r: {} for r in shells}
        entry_by: dict[int, dict[tuple[int, int], mpfr]] = {r: {} for r in shells}
        with ProcessPoolExecutor(max_workers=workers) as ex:
            futs = [ex.submit(_fro_worker_chunk, c) for c in chunks]
            for fut in as_completed(futs):
                res = fut.result()
                print(f"    chunk a=[{res['a_start']},{res['a_end']}) updates={res['n_full']}", flush=True)
                n_full += res["n_full"]
                for rs, rd in res["row_by"].items():
                    sh = int(rs)
                    merge_sum_dicts(row_by[sh], {int(k): mpfr(v, precision=prec) for k, v in rd.items()}, prec=prec)
                for rs, cd in res["col_by"].items():
                    sh = int(rs)
                    merge_sum_dicts(col_by[sh], {int(k): mpfr(v, precision=prec) for k, v in cd.items()}, prec=prec)
                for rs, ed in res["entry_by"].items():
                    sh = int(rs)
                    src = {
                        (int(k.split(",")[0]), int(k.split(",")[1])): mpfr(v, precision=prec)
                        for k, v in ed.items()
                    }
                    _merge_entry_maps(entry_by[sh], src, prec=prec)
                D = res["D"]
        fro_by, fro_full, n_by = _finalize_frobenius_from_entries(entry_by, shells, prec=prec)
        return fro_by, n_by, fro_full, n_full, D, row_by, col_by

    _ctx(prec)
    if checkpoint_path.is_file():
        checkpoint_path.unlink()

    b = build_hermitian_basis(n, radii)
    by = shell_modes_by_r(n)
    out_set: set[tuple[int, int, int]] = set()
    for r in radii:
        if r in by:
            out_set.update(by[r])
    hit = _mode_to_basis(b)
    D = b.D
    mode_shell = [_shell_of_mode(b, m) for m in range(D)]
    n_pairs = D * D
    row_by: dict[int, dict[int, mpfr]] = {r: {} for r in shells}
    col_by: dict[int, dict[int, mpfr]] = {r: {} for r in shells}
    entry_by: dict[int, dict[tuple[int, int], mpfr]] = {r: {} for r in shells}

    def progress_cb(done: int, total: int) -> None:
        if progress_every and done % progress_every == 0:
            print(f"    pairs {done}/{total}", flush=True)

    fro_by, n_by, fro_full, n_full, _ = _process_pair_range(
        b,
        out_set,
        hit,
        mode_shell,
        shells,
        0,
        row_by,
        col_by,
        entry_by,
        a_start=0,
        a_end=D,
        prec=prec,
        progress_cb=progress_cb,
        pair_offset=0,
        n_pairs=n_pairs,
        rational_coefficients=rational_coefficients,
    )
    return fro_by, n_by, fro_full, n_full, D, row_by, col_by


def _bound_from_fro_sq(
    fro_sq: mpfr,
    n_updates: int,
    n: int,
    D: int,
    *,
    output_shell: int | None,
    prec: int,
) -> CertifiedOpNormBound:
    L_hi = sqrt_up(fro_sq)
    two_sqrt2 = sqrt_up(mpfr(8, precision=prec))
    C_hi = mul_up(two_sqrt2, L_hi)
    payload = {"n": n, "shell": output_shell, "fro_sq": str(fro_sq), "n_updates": n_updates}
    digest = hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()
    return CertifiedOpNormBound(
        shell=output_shell,
        n=n,
        D=D,
        frobenius_sq_hi=str(fro_sq),
        L_op_hi=str(L_hi),
        C_term_hi=str(C_hi),
        n_updates=n_updates,
        method="frobenius_mpfr_up",
        evidence_level="N5",
        payload_sha256=digest,
    )


def certified_L_op_frobenius_hi(
    n: int,
    radii: tuple[int, ...],
    *,
    output_shell: int | None = None,
    prec: int = 200,
) -> CertifiedOpNormBound:
    fro_sq, n_updates, D = _accumulate_frobenius_sq_hi(
        n, radii, output_shell=output_shell, prec=prec
    )
    return _bound_from_fro_sq(
        fro_sq, n_updates, n, D, output_shell=output_shell, prec=prec
    )


def certified_shell_blocks(
    n: int = 24,
    radii: tuple[int, ...] | None = None,
    prec: int = 200,
    *,
    workers: int = 1,
    bound_method: str | None = None,
    rational_coefficients: bool = False,
) -> dict:
    bound_method = bound_method or DEFAULT_BOUND_METHOD
    if bound_method in ("frobenius", "best", "min_fro_one_inf"):
        if not FROBENIUS_STREAMING_ENABLED:
            raise ValueError(
                f"bound_method={bound_method!r} disabled: Frobenius streaming sums |delta|^2 "
                "not |M_ij|^2 per entry. Use one_inf_only."
            )
    radii = radii or all_dealias_radii(n)
    shells = sorted(set(radii))
    print(f"  combined pass 1-inf only ({len(shells)} shells)...", flush=True)
    fro_by, n_by, fro_full, n_full, D, row_by, col_by = (
        _accumulate_frobenius_all_shells_single_pass(
            n, radii, prec=prec, workers=workers, rational_coefficients=rational_coefficients
        )
    )
    print(f"  D={D}, full updates={n_full}", flush=True)
    blocks = []
    for r in shells:
        bo = _bound_from_one_inf(row_by[r], col_by[r], n_by[r], n, D, output_shell=r, prec=prec)
        bf = None
        if FROBENIUS_STREAMING_ENABLED and bound_method != "one_inf_only":
            bf = _bound_from_fro_sq(fro_by[r], n_by[r], n, D, output_shell=r, prec=prec)
            bnd = _pick_tighter_bound(bf, bo, prec=prec) if bound_method == "best" else bf
        else:
            bnd = bo
        entry = {**bnd.as_dict(), "C_term_one_inf_hi": bo.C_term_hi}
        if bf is not None:
            entry["C_term_frobenius_hi"] = bf.C_term_hi
        blocks.append(entry)
    full_o = _bound_from_one_inf(
        row_by[shells[0]], col_by[shells[0]], n_by[shells[0]], n, D, output_shell=None, prec=prec
    )
    full_f = _bound_from_fro_sq(fro_full, n_full, n, D, output_shell=None, prec=prec)
    if bound_method == "one_inf_only":
        full_best = full_o
    elif bound_method == "frobenius":
        full_best = full_f
    else:
        full_best = _pick_tighter_bound(full_f, full_o, prec=prec)
    return {
        "n": n,
        "n_shells": len(shells),
        "radii_count": len(radii),
        "bound_method": bound_method,
        "rational_coefficients": rational_coefficients,
        "blocks": blocks,
        "full_at_T_one_inf": full_o.as_dict(),
        "full_at_T_frobenius": full_f.as_dict(),
        "full_at_T_best": full_best.as_dict(),
        "arithmetic_backend": ArithmeticBackend(precision_bits=prec).as_dict(),
        "convention": (
            "M_r: output on shell r only. Repair baseline: ||M_r||_op <= ||M_r||_{1,inf} "
            "(MPFR upward). Frobenius uses entry aggregation then |M_ij|^2 sum."
            + (" Rational mpfr_leray + mpfr_vdot_real coefficients." if rational_coefficients else "")
        ),
    }


def _process_pair_range_cluster_frobenius(
    b,
    out_set: set,
    hit: dict,
    cluster_shells: frozenset[int],
    fro_sq: mpfr,
    n_updates: int,
    *,
    a_start: int,
    a_end: int,
    prec: int,
    rational_coefficients: bool = False,
) -> tuple[mpfr, int]:
    D = b.D
    sq2 = math.sqrt(2.0)
    t, nu, T = FROZEN.T, FROZEN.nu, FROZEN.T
    mode_shell = [_shell_of_mode(b, m) for m in range(D)]
    batch: list[float] = []

    def _flush_batch() -> None:
        nonlocal fro_sq, n_updates
        for delta in batch:
            fro_sq = add_up(fro_sq, _delta_sq_hi(delta, prec))
            n_updates += 1
        batch.clear()

    for a in range(a_start, a_end):
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
            if not Nloc:
                continue
            _apply_leray_nloc(Nloc, rational=rational_coefficients, prec=prec)
            ms: set[int] = set()
            for s in Nloc:
                ms.update(hit.get(s, ()))
            for m in ms:
                if mode_shell[m] not in cluster_shells:
                    continue
                g = _inner_psi_N_weighted_certified(
                    b, m, Nloc, t, nu, T, output_shells=cluster_shells,
                    rational=rational_coefficients, prec=prec,
                )
                if g == 0.0:
                    continue
                contrib = g / 6.0
                for p, q, r in (
                    (m, a, bb),
                    (m, bb, a),
                    (a, m, bb),
                    (a, bb, m),
                    (bb, m, a),
                    (bb, a, m),
                ):
                    batch.append(0.5 * contrib if q == r else 0.5 * sq2 * contrib)
                if len(batch) >= _MPFR_BATCH:
                    _flush_batch()
    if batch:
        _flush_batch()
    return fro_sq, n_updates


def _process_pair_range_cluster_one_inf(
    b,
    out_set: set,
    hit: dict,
    cluster_shells: frozenset[int],
    row_sum: list[mpfr],
    col_sum: dict[int, mpfr],
    n_updates: int,
    *,
    a_start: int,
    a_end: int,
    prec: int,
    rational_coefficients: bool = False,
) -> int:
    """1-inf accumulation for M_cluster = sum_{r in cluster} M_r (combined output)."""
    D = b.D
    sq2 = math.sqrt(2.0)
    t, nu, T = FROZEN.T, FROZEN.nu, FROZEN.T
    mode_shell = [_shell_of_mode(b, m) for m in range(D)]
    batch: list[tuple[int, int, int, float]] = []

    def _flush() -> None:
        nonlocal n_updates
        z = mpfr(0, precision=prec)
        for p, q, r, delta in batch:
            n_updates += 1
            ad = abs_up(delta, prec)
            if q == r:
                col = q
            else:
                lo, hi = (q, r) if q < r else (r, q)
                col = _sym_col_index(lo, hi, D)
            row_sum[p] = add_up(row_sum[p], ad)
            prev = col_sum.get(col, z)
            col_sum[col] = add_up(prev, ad)
        batch.clear()

    for a in range(a_start, a_end):
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
            if not Nloc:
                continue
            _apply_leray_nloc(Nloc, rational=rational_coefficients, prec=prec)
            ms: set[int] = set()
            for s in Nloc:
                ms.update(hit.get(s, ()))
            for m in ms:
                if mode_shell[m] not in cluster_shells:
                    continue
                g = _inner_psi_N_weighted_certified(
                    b, m, Nloc, t, nu, T, output_shells=cluster_shells,
                    rational=rational_coefficients, prec=prec,
                )
                if g == 0.0:
                    continue
                contrib = g / 6.0
                for p, q, r in (
                    (m, a, bb),
                    (m, bb, a),
                    (a, m, bb),
                    (a, bb, m),
                    (bb, m, a),
                    (bb, a, m),
                ):
                    batch.append(
                        (p, q, r, 0.5 * contrib if q == r else 0.5 * sq2 * contrib)
                    )
                if len(batch) >= _MPFR_BATCH:
                    _flush()
    if batch:
        _flush()
    return n_updates


def certified_cluster_one_inf_hi(
    n: int,
    radii: tuple[int, ...],
    cluster_shells: tuple[int, ...],
    *,
    prec: int = 200,
) -> CertifiedOpNormBound:
    """Certified 1-inf upper bound for combined output on cluster shells."""
    _ctx(prec)
    cluster = frozenset(cluster_shells)
    b = build_hermitian_basis(n, radii)
    by = shell_modes_by_r(n)
    out_set: set[tuple[int, int, int]] = set()
    for r in radii:
        if r in by:
            out_set.update(by[r])
    hit = _mode_to_basis(b)
    D = b.D
    row_sum = [mpfr(0, precision=prec) for _ in range(D)]
    col_sum: dict[int, mpfr] = {}
    n_updates = 0
    n_updates = _process_pair_range_cluster_one_inf(
        b,
        out_set,
        hit,
        cluster,
        row_sum,
        col_sum,
        n_updates,
        a_start=0,
        a_end=D,
        prec=prec,
    )
    label = min(cluster) if len(cluster) == 1 else None
    return _bound_from_one_inf_rowlist(
        row_sum, col_sum, n_updates, n, D, output_shell=label, prec=prec
    )


def certified_cluster_frobenius_hi(
    n: int,
    radii: tuple[int, ...],
    cluster_shells: tuple[int, ...],
    *,
    prec: int = 200,
) -> CertifiedOpNormBound:
    """Certified Frobenius upper bound for combined output on cluster shells."""
    _ctx(prec)
    cluster = frozenset(cluster_shells)
    b = build_hermitian_basis(n, radii)
    by = shell_modes_by_r(n)
    out_set: set[tuple[int, int, int]] = set()
    for r in radii:
        if r in by:
            out_set.update(by[r])
    hit = _mode_to_basis(b)
    D = b.D
    fro_sq = mpfr(0, precision=prec)
    n_updates = 0
    fro_sq, n_updates = _process_pair_range_cluster_frobenius(
        b,
        out_set,
        hit,
        cluster,
        fro_sq,
        n_updates,
        a_start=0,
        a_end=D,
        prec=prec,
    )
    label = min(cluster) if len(cluster) == 1 else None
    return _bound_from_fro_sq(
        fro_sq, n_updates, n, D, output_shell=label, prec=prec
    )


def _flush_entry_batch_cluster_best(
    batch: list[tuple[int, int, int, float, int | None]],
    fro_cluster: dict[int, mpfr],
    n_cluster: dict[int, int],
    row_by: dict[int, list[mpfr]],
    col_by: dict[int, dict[int, mpfr]],
    shell_to_cluster: dict[int, int],
    *,
    prec: int,
    D: int,
) -> None:
    z = mpfr(0, precision=prec)
    for p, q, r, delta, sr in batch:
        if sr is None:
            continue
        cid = shell_to_cluster.get(sr)
        if cid is None:
            continue
        val = _delta_sq_hi(delta, prec)
        fro_cluster[cid] = add_up(fro_cluster[cid], val)
        n_cluster[cid] += 1
        ad = abs_up(delta, prec)
        if q == r:
            col = q
        else:
            lo, hi = (q, r) if q < r else (r, q)
            col = _sym_col_index(lo, hi, D)
        rows = row_by[cid]
        cols = col_by[cid]
        rows[p] = add_up(rows[p], ad)
        prev = cols.get(col, z)
        cols[col] = add_up(prev, ad)
    batch.clear()


def _process_pair_range_cluster_best(
    b,
    out_set: set,
    hit: dict,
    mode_shell: list[int],
    shell_to_cluster: dict[int, int],
    fro_cluster: dict[int, mpfr],
    n_cluster: dict[int, int],
    row_by: dict[int, list[mpfr]],
    col_by: dict[int, dict[int, mpfr]],
    *,
    a_start: int,
    a_end: int,
    prec: int,
    progress_cb=None,
    pair_offset: int = 0,
    n_pairs: int = 0,
    rational_coefficients: bool = False,
) -> int:
    D = b.D
    sq2 = math.sqrt(2.0)
    t, nu, T = FROZEN.T, FROZEN.nu, FROZEN.T
    done = pair_offset
    batch: list[tuple[int, int, int, float, int | None]] = []

    for a in range(a_start, a_end):
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
            done += 1
            if not Nloc:
                if progress_cb and n_pairs and done % 250_000 == 0:
                    progress_cb(done, n_pairs)
                continue
            _apply_leray_nloc(Nloc, rational=rational_coefficients, prec=prec)
            ms: set[int] = set()
            for s in Nloc:
                ms.update(hit.get(s, ()))
            for m in ms:
                g = _inner_psi_N_weighted_certified(
                    b, m, Nloc, t, nu, T, output_shell=None,
                    rational=rational_coefficients, prec=prec,
                )
                if g == 0.0:
                    continue
                contrib = g / 6.0
                sr = mode_shell[m]
                sr_tag = sr if sr in shell_to_cluster else None
                for p, q, r in (
                    (m, a, bb),
                    (m, bb, a),
                    (a, m, bb),
                    (a, bb, m),
                    (bb, m, a),
                    (bb, a, m),
                ):
                    delta = 0.5 * contrib if q == r else 0.5 * sq2 * contrib
                    batch.append((p, q, r, delta, sr_tag))
                    if len(batch) >= _MPFR_BATCH:
                        _flush_entry_batch_cluster_best(
                            batch,
                            fro_cluster,
                            n_cluster,
                            row_by,
                            col_by,
                            shell_to_cluster,
                            prec=prec,
                            D=D,
                        )
            if progress_cb and n_pairs and done % 250_000 == 0:
                progress_cb(done, n_pairs)
    if batch:
        _flush_entry_batch_cluster_best(
            batch,
            fro_cluster,
            n_cluster,
            row_by,
            col_by,
            shell_to_cluster,
            prec=prec,
            D=D,
        )
    return done


def _cluster_best_worker_chunk(args: tuple) -> dict:
    a_start, a_end, n, radii, prec, shell_to_cluster = args
    _ctx(prec)
    b = build_hermitian_basis(n, radii)
    by = shell_modes_by_r(n)
    out_set: set[tuple[int, int, int]] = set()
    for r in radii:
        if r in by:
            out_set.update(by[r])
    hit = _mode_to_basis(b)
    mode_shell = [_shell_of_mode(b, m) for m in range(b.D)]
    D = b.D
    cluster_ids = sorted(set(shell_to_cluster.values()))
    fro_cluster = {c: mpfr(0, precision=prec) for c in cluster_ids}
    n_cluster = {c: 0 for c in cluster_ids}
    row_by = {c: [mpfr(0, precision=prec) for _ in range(D)] for c in cluster_ids}
    col_by: dict[int, dict[int, mpfr]] = {c: {} for c in cluster_ids}
    _process_pair_range_cluster_best(
        b,
        out_set,
        hit,
        mode_shell,
        shell_to_cluster,
        fro_cluster,
        n_cluster,
        row_by,
        col_by,
        a_start=a_start,
        a_end=a_end,
        prec=prec,
    )
    return {
        "a_start": a_start,
        "a_end": a_end,
        "fro_cluster": {str(k): str(v) for k, v in fro_cluster.items()},
        "n_cluster": n_cluster,
        "row_by": {str(c): [str(v) for v in rows] for c, rows in row_by.items()},
        "col_by": {
            str(c): {str(k): str(v) for k, v in cols.items()} for c, cols in col_by.items()
        },
        "D": D,
    }


def _accumulate_cluster_bounds_single_pass(
    n: int,
    radii: tuple[int, ...],
    shell_to_cluster: dict[int, int],
    *,
    prec: int = 128,
    workers: int = 4,
) -> tuple[dict[int, mpfr], dict[int, list[mpfr]], dict[int, dict[int, mpfr]], dict[int, int], int]:
    """One D^2 pass: Frobenius + 1-inf accumulators per cluster."""
    _ctx(prec)
    cluster_ids = sorted(set(shell_to_cluster.values()))
    fro_cluster = {c: mpfr(0, precision=prec) for c in cluster_ids}
    n_cluster = {c: 0 for c in cluster_ids}
    row_by: dict[int, list[mpfr]] = {}
    col_by: dict[int, dict[int, mpfr]] = {c: {} for c in cluster_ids}

    if workers > 1:
        print(f"    cluster best pass workers={workers}", flush=True)
        b_tmp = build_hermitian_basis(n, radii)
        D = b_tmp.D
        row_by = {c: [mpfr(0, precision=prec) for _ in range(D)] for c in cluster_ids}
        del b_tmp
        chunk_size = max(1, (D + workers - 1) // workers)
        chunks = []
        stc = dict(shell_to_cluster)
        for i in range(workers):
            a0 = i * chunk_size
            if a0 >= D:
                break
            a1 = min(D, a0 + chunk_size)
            chunks.append((a0, a1, n, radii, prec, stc))
        with ProcessPoolExecutor(max_workers=workers) as ex:
            futs = [ex.submit(_cluster_best_worker_chunk, c) for c in chunks]
            for fut in as_completed(futs):
                res = fut.result()
                print(
                    f"    cluster-best chunk a=[{res['a_start']},{res['a_end']})",
                    flush=True,
                )
                for cid_s, val in res["fro_cluster"].items():
                    cid = int(cid_s)
                    fro_cluster[cid] = add_up(fro_cluster[cid], mpfr(val, precision=prec))
                    n_cluster[cid] += res["n_cluster"][cid]
                for cid_s, rows in res["row_by"].items():
                    cid = int(cid_s)
                    for i, v in enumerate(rows):
                        row_by[cid][i] = add_up(row_by[cid][i], mpfr(v, precision=prec))
                for cid_s, cols in res["col_by"].items():
                    cid = int(cid_s)
                    merge_sum_dicts(
                        col_by[cid],
                        {int(k): mpfr(v, precision=prec) for k, v in cols.items()},
                        prec=prec,
                    )
                D = res["D"]
        return fro_cluster, row_by, col_by, n_cluster, D

    b = build_hermitian_basis(n, radii)
    by = shell_modes_by_r(n)
    out_set: set[tuple[int, int, int]] = set()
    for r in radii:
        if r in by:
            out_set.update(by[r])
    hit = _mode_to_basis(b)
    D = b.D
    row_by = {c: [mpfr(0, precision=prec) for _ in range(D)] for c in cluster_ids}
    mode_shell = [_shell_of_mode(b, m) for m in range(D)]
    n_pairs = D * D

    def progress_cb(done: int, total: int) -> None:
        print(f"    cluster-best pairs {done}/{total}", flush=True)

    _process_pair_range_cluster_best(
        b,
        out_set,
        hit,
        mode_shell,
        shell_to_cluster,
        fro_cluster,
        n_cluster,
        row_by,
        col_by,
        a_start=0,
        a_end=D,
        prec=prec,
        progress_cb=progress_cb,
        n_pairs=n_pairs,
    )
    return fro_cluster, row_by, col_by, n_cluster, D


def _flush_entry_batch_clusters(
    batch: list[tuple[int, int, int, float, int | None]],
    fro_cluster: dict[int, mpfr],
    n_cluster: dict[int, int],
    shell_to_cluster: dict[int, int],
    *,
    prec: int,
) -> None:
    for _p, _q, _r, delta, sr in batch:
        if sr is None:
            continue
        cid = shell_to_cluster.get(sr)
        if cid is None:
            continue
        val = _delta_sq_hi(delta, prec)
        fro_cluster[cid] = add_up(fro_cluster[cid], val)
        n_cluster[cid] += 1
    batch.clear()


def _process_pair_range_clusters(
    b,
    out_set: set,
    hit: dict,
    mode_shell: list[int],
    shell_to_cluster: dict[int, int],
    fro_cluster: dict[int, mpfr],
    n_cluster: dict[int, int],
    *,
    a_start: int,
    a_end: int,
    prec: int,
    progress_cb=None,
    pair_offset: int = 0,
    n_pairs: int = 0,
    rational_coefficients: bool = False,
) -> tuple[dict[int, mpfr], dict[int, int], int]:
    D = b.D
    sq2 = math.sqrt(2.0)
    t, nu, T = FROZEN.T, FROZEN.nu, FROZEN.T
    done = pair_offset
    batch: list[tuple[int, int, int, float, int | None]] = []

    for a in range(a_start, a_end):
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
            done += 1
            if not Nloc:
                if progress_cb and n_pairs and done % 250_000 == 0:
                    progress_cb(done, n_pairs)
                continue
            _apply_leray_nloc(Nloc, rational=rational_coefficients, prec=prec)
            ms: set[int] = set()
            for s in Nloc:
                ms.update(hit.get(s, ()))
            for m in ms:
                g = _inner_psi_N_weighted_certified(
                    b, m, Nloc, t, nu, T, output_shell=None,
                    rational=rational_coefficients, prec=prec,
                )
                if g == 0.0:
                    continue
                contrib = g / 6.0
                sr = mode_shell[m]
                sr_tag = sr if sr in shell_to_cluster else None
                for p, q, r in (
                    (m, a, bb),
                    (m, bb, a),
                    (a, m, bb),
                    (a, bb, m),
                    (bb, m, a),
                    (bb, a, m),
                ):
                    delta = 0.5 * contrib if q == r else 0.5 * sq2 * contrib
                    batch.append((p, q, r, delta, sr_tag))
                    if len(batch) >= _MPFR_BATCH:
                        _flush_entry_batch_clusters(
                            batch, fro_cluster, n_cluster, shell_to_cluster, prec=prec
                        )
            if progress_cb and n_pairs and done % 250_000 == 0:
                progress_cb(done, n_pairs)
    if batch:
        _flush_entry_batch_clusters(
            batch, fro_cluster, n_cluster, shell_to_cluster, prec=prec
        )
    return fro_cluster, n_cluster, done


def _cluster_worker_chunk(args: tuple) -> dict:
    a_start, a_end, n, radii, prec, shell_to_cluster = args
    _ctx(prec)
    b = build_hermitian_basis(n, radii)
    by = shell_modes_by_r(n)
    out_set: set[tuple[int, int, int]] = set()
    for r in radii:
        if r in by:
            out_set.update(by[r])
    hit = _mode_to_basis(b)
    mode_shell = [_shell_of_mode(b, m) for m in range(b.D)]
    cluster_ids = sorted(set(shell_to_cluster.values()))
    fro_cluster = {c: mpfr(0, precision=prec) for c in cluster_ids}
    n_cluster = {c: 0 for c in cluster_ids}
    fro_cluster, n_cluster, _ = _process_pair_range_clusters(
        b,
        out_set,
        hit,
        mode_shell,
        shell_to_cluster,
        fro_cluster,
        n_cluster,
        a_start=a_start,
        a_end=a_end,
        prec=prec,
    )
    return {
        "a_start": a_start,
        "a_end": a_end,
        "fro_cluster": {str(k): str(v) for k, v in fro_cluster.items()},
        "n_cluster": n_cluster,
        "D": b.D,
    }


def _accumulate_frobenius_clusters_single_pass(
    n: int,
    radii: tuple[int, ...],
    shell_to_cluster: dict[int, int],
    *,
    prec: int = 128,
    workers: int = 4,
    progress_every: int = 250_000,
) -> tuple[dict[int, mpfr], dict[int, int], int]:
    """One D^2 pass: Frobenius per cluster (combined output on cluster shells)."""
    _ctx(prec)
    cluster_ids = sorted(set(shell_to_cluster.values()))
    fro_cluster = {c: mpfr(0, precision=prec) for c in cluster_ids}
    n_cluster = {c: 0 for c in cluster_ids}

    if workers > 1:
        print(f"    cluster pass workers={workers}", flush=True)
        b_tmp = build_hermitian_basis(n, radii)
        D = b_tmp.D
        del b_tmp
        chunk_size = max(1, (D + workers - 1) // workers)
        chunks = []
        stc = dict(shell_to_cluster)
        for i in range(workers):
            a0 = i * chunk_size
            if a0 >= D:
                break
            a1 = min(D, a0 + chunk_size)
            chunks.append((a0, a1, n, radii, prec, stc))
        with ProcessPoolExecutor(max_workers=workers) as ex:
            futs = [ex.submit(_cluster_worker_chunk, c) for c in chunks]
            for fut in as_completed(futs):
                res = fut.result()
                print(
                    f"    cluster chunk a=[{res['a_start']},{res['a_end']})",
                    flush=True,
                )
                for cid_s, val in res["fro_cluster"].items():
                    cid = int(cid_s)
                    fro_cluster[cid] = add_up(fro_cluster[cid], mpfr(val, precision=prec))
                    n_cluster[cid] += res["n_cluster"][cid]
                D = res["D"]
        return fro_cluster, n_cluster, D

    b = build_hermitian_basis(n, radii)
    by = shell_modes_by_r(n)
    out_set: set[tuple[int, int, int]] = set()
    for r in radii:
        if r in by:
            out_set.update(by[r])
    hit = _mode_to_basis(b)
    D = b.D
    mode_shell = [_shell_of_mode(b, m) for m in range(D)]
    n_pairs = D * D

    def progress_cb(done: int, total: int) -> None:
        print(f"    cluster pairs {done}/{total}", flush=True)

    _process_pair_range_clusters(
        b,
        out_set,
        hit,
        mode_shell,
        shell_to_cluster,
        fro_cluster,
        n_cluster,
        a_start=0,
        a_end=D,
        prec=prec,
        progress_cb=progress_cb,
        n_pairs=n_pairs,
    )
    return fro_cluster, n_cluster, D


def compare_shell_vs_cluster_band(
    n: int = 24,
    radii: tuple[int, ...] | None = None,
    *,
    prec: int = 128,
) -> dict:
    """
    Band diagnostic: per-shell Frobenius sum vs single-cluster Frobenius.

    Cluster bound on ||sum_{r in C} M_r||_F is always <= sum_r ||M_r||_F
    (triangle inequality). Tighter cluster norms may reduce Route A if
    reweighted per-shell factors are handled (future work).
    """
    radii = radii or tuple(range(1, 7))
    shells = sorted(set(radii))
    per_shell = []
    per_shell_C_sum = mpfr(0, precision=prec)
    for r in shells:
        bnd = certified_L_op_frobenius_hi(n, radii, output_shell=r, prec=prec)
        per_shell.append({"shell": r, "C_term_hi": bnd.C_term_hi, "L_op_hi": bnd.L_op_hi})
        per_shell_C_sum = add_up(per_shell_C_sum, mpfr(bnd.C_term_hi, precision=prec))
    cluster = certified_cluster_frobenius_hi(n, radii, tuple(shells), prec=prec)
    cluster_C = mpfr(cluster.C_term_hi, precision=prec)
    ratio = float(cluster_C / per_shell_C_sum) if per_shell_C_sum > mpfr(0) else 0.0
    return {
        "n": n,
        "radii": list(shells),
        "per_shell_C_sum_hi": str(per_shell_C_sum),
        "cluster_C_term_hi": cluster.C_term_hi,
        "cluster_L_op_hi": cluster.L_op_hi,
        "cluster_tighter_ratio": ratio,
        "per_shell": per_shell,
        "note": (
            "cluster_tighter_ratio = cluster_C / sum(per_shell C); "
            "<1 means combined Frobenius beats triangle sum at t=T"
        ),
    }
