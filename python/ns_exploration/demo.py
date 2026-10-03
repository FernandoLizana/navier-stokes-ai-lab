"""
Public demo: small Taylor–Green run on T³ (exploratory floating-point numerics).

Usage (from repository root, after installing the package):

    python -m ns_exploration.demo

Does not prove regularity, blow-up, or any continuum PDE statement.
"""

from __future__ import annotations

import argparse
import csv
import json
import platform
import sys
import time
from pathlib import Path

import numpy as np

from ns_exploration import __version__
from ns_exploration.initial_conditions.generators import taylor_green
from ns_exploration.spectral.fourier_conventions import ifft_vector
from ns_exploration.spectral.solver import NavierStokesSolver, SolverConfig

# Conservative public defaults (fast on a laptop CPU).
_DEFAULT_N = 8
_DEFAULT_NU = 0.05
_DEFAULT_DT = 0.001
_DEFAULT_T_END = 0.02
_DEFAULT_INTEGRATOR = "rk4"
_MAX_N = 32
_MAX_STEPS = 500
_MAX_WALL_S = 120.0


def _validate(args: argparse.Namespace) -> None:
    if args.n < 4 or args.n > _MAX_N or args.n % 2 != 0:
        raise SystemExit(f"n must be even in [4, {_MAX_N}], got {args.n}")
    if args.nu <= 0:
        raise SystemExit("nu must be positive")
    if args.dt <= 0 or args.t_end <= 0:
        raise SystemExit("dt and t_end must be positive")
    if args.t_end < args.dt:
        raise SystemExit("t_end must be >= dt")
    nsteps = int(np.ceil(args.t_end / args.dt))
    if nsteps > _MAX_STEPS:
        raise SystemExit(
            f"requested {nsteps} steps exceeds demo budget {_MAX_STEPS}; "
            "lower t_end or raise dt"
        )
    if args.integrator not in ("rk4", "etd_rk2", "semi_implicit_euler"):
        raise SystemExit(
            f"integrator must be rk4|etd_rk2|semi_implicit_euler, got {args.integrator}"
        )


def _write_csv(path: Path, history: list[dict]) -> None:
    if not history:
        raise RuntimeError("empty history")
    keys = list(history[0].keys())
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for row in history:
            w.writerow(row)


def _plot(path: Path, history: list[dict], u_slice: np.ndarray, *, n: int) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    t = [h["t"] for h in history]
    energy = [h["energy"] for h in history]
    enstrophy = [h["enstrophy"] for h in history]

    fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.4))
    axes[0].plot(t, energy, color="#1f4e79")
    axes[0].set_xlabel("t (dimensionless)")
    axes[0].set_ylabel("kinetic energy E (dimensionless)")
    axes[0].set_title("Energy vs time")

    axes[1].plot(t, enstrophy, color="#8b3a3a")
    axes[1].set_xlabel("t (dimensionless)")
    axes[1].set_ylabel("enstrophy Ω (dimensionless)")
    axes[1].set_title("Enstrophy vs time")

    # Mid-plane speed |u| at z = n//2 — 2D slice of a 3D field.
    speed = np.sqrt(np.sum(u_slice**2, axis=0))
    im = axes[2].imshow(speed.T, origin="lower", cmap="viridis", aspect="equal")
    axes[2].set_title(f"|u| mid-plane (z={n // 2}), 2D slice of 3D field")
    axes[2].set_xlabel("x index")
    axes[2].set_ylabel("y index")
    fig.colorbar(im, ax=axes[2], fraction=0.046, pad=0.04)

    fig.suptitle(
        "Navier–Stokes AI Lab demo — exploratory floating-point numerics (not a proof)",
        fontsize=10,
    )
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)


def run_demo(args: argparse.Namespace) -> dict:
    _validate(args)
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    cfg = SolverConfig(
        n=args.n,
        nu=args.nu,
        dt=args.dt,
        t_end=args.t_end,
        integrator=args.integrator,
        dealias=True,
        force_hat=None,
        route="B",
        evidence_level="N2",
    )
    uh, meta = taylor_green(args.n, a=args.amplitude)
    t0 = time.perf_counter()
    state = NavierStokesSolver(cfg).run(uh)
    wall = time.perf_counter() - t0
    if wall > _MAX_WALL_S:
        # Soft note only; run already finished.
        pass

    e0 = float(state.history[0]["energy"])
    e1 = float(state.history[-1]["energy"])
    div = float(state.history[-1]["div_l2"])
    finished = state.stopped_reason is None

    summary = {
        "project": "Navier–Stokes AI Lab",
        "package_version": __version__,
        "evidence_level": "N2",
        "disclaimer": (
            "Exploratory floating-point simulation on a truncated Galerkin/"
            "pseudospectral model on T^3. Not a continuum PDE proof."
        ),
        "config": {
            "n": args.n,
            "nu": args.nu,
            "dt": args.dt,
            "t_end": args.t_end,
            "integrator": args.integrator,
            "ic": "taylor_green",
            "amplitude": args.amplitude,
            "dealias": True,
            "force": None,
            "route": "B",
        },
        "ic_metadata": {
            "energy": meta.energy,
            "name": meta.name,
        },
        "result": {
            "finished_normally": finished,
            "stopped_reason": state.stopped_reason,
            "steps": state.step,
            "t_final": state.t,
            "energy_initial": e0,
            "energy_final": e1,
            "div_l2_final": div,
            "wall_seconds": wall,
        },
        "environment": {
            "python": sys.version.split()[0],
            "platform": platform.platform(),
            "numpy": np.__version__,
        },
        "outputs": {
            "diagnostics_csv": "diagnostics.csv",
            "config_json": "config.json",
            "summary_json": "summary.json",
            "figure_png": "diagnostics.png",
        },
    }

    _write_csv(out / "diagnostics.csv", state.history)
    (out / "config.json").write_text(
        json.dumps(summary["config"], indent=2), encoding="utf-8"
    )
    (out / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    u = ifft_vector(state.u_hat)
    mid = args.n // 2
    u_slice = u[:, :, :, mid]
    _plot(out / "diagnostics.png", state.history, u_slice, n=args.n)

    return summary


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Small Taylor–Green Navier–Stokes demo (exploratory, not a proof)"
    )
    p.add_argument("--n", type=int, default=_DEFAULT_N, help="grid size per axis (even)")
    p.add_argument("--nu", type=float, default=_DEFAULT_NU, help="viscosity")
    p.add_argument("--dt", type=float, default=_DEFAULT_DT, help="time step")
    p.add_argument("--t-end", type=float, default=_DEFAULT_T_END, help="final time")
    p.add_argument(
        "--integrator",
        type=str,
        default=_DEFAULT_INTEGRATOR,
        help="rk4 | etd_rk2 | semi_implicit_euler",
    )
    p.add_argument("--amplitude", type=float, default=1.0, help="Taylor–Green amplitude")
    p.add_argument(
        "--out-dir",
        type=str,
        default="demos/out/latest",
        help="output directory (created if missing)",
    )
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        summary = run_demo(args)
    except SystemExit:
        raise
    except Exception as exc:  # noqa: BLE001 — CLI boundary
        print(f"demo failed: {exc}", file=sys.stderr)
        return 1

    r = summary["result"]
    print("Navier-Stokes AI Lab - demo complete (exploratory numerics)")
    print(f"  out: {args.out_dir}")
    print(f"  finished_normally: {r['finished_normally']}")
    if r["stopped_reason"]:
        print(f"  stopped_reason: {r['stopped_reason']}")
    print(f"  steps: {r['steps']}  t_final: {r['t_final']}")
    print(f"  energy: {r['energy_initial']:.10g} -> {r['energy_final']:.10g}")
    print(f"  div_l2_final: {r['div_l2_final']:.6e}")
    print(f"  wall_seconds: {r['wall_seconds']:.3f}")
    print("  NOTE: not a continuum proof; evidence_level=N2")
    # Early safety stop is still a successful demo run (outputs written).
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
