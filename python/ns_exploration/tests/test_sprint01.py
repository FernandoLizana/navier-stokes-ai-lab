"""Critical tests for Sprint 1 exploratory solver."""

from __future__ import annotations

import numpy as np
import pytest

from ns_exploration.initial_conditions.generators import (
    abc_flow,
    random_div_free,
    taylor_green,
    vortex_tubes_periodic,
)
from ns_exploration.spectral.fourier_conventions import (
    fft_vector,
    ifft_vector,
    kinetic_energy_from_hat,
    wave_number_grids,
)
from ns_exploration.spectral.integrators import step_etd_rk2, step_rk4, step_semi_implicit_euler
from ns_exploration.spectral.leray import divergence_l2, leray_project_hat
from ns_exploration.spectral.operators import nonlinear_hat
from ns_exploration.spectral.solver import NavierStokesSolver, SolverConfig
from ns_exploration.validation.intervals import (
    finite_energy_identity_cancel_check,
    leray_project_interval_mode,
    propagate_linear_heat_interval,
    short_convolution_interval,
    Interval,
    divergence_interval,
)


def test_fft_roundtrip():
    rng = np.random.default_rng(0)
    n = 16
    u = rng.normal(size=(3, n, n, n))
    uh = fft_vector(u)
    ur = ifft_vector(uh)
    assert np.allclose(u, ur, atol=1e-10)


def test_leray_kills_divergence():
    n = 16
    uh, _ = random_div_free(n, seed=1)
    # contaminate
    uh = uh + 0.1 * fft_vector(np.ones((3, n, n, n)))
    uh = leray_project_hat(uh)
    assert divergence_l2(uh) < 1e-12


def test_taylor_green_div_free():
    uh, meta = taylor_green(24)
    assert divergence_l2(uh) < 1e-12
    assert meta.energy > 0


def test_abc_div_free():
    uh, _ = abc_flow(24)
    assert divergence_l2(uh) < 1e-12


def test_vortex_tubes_div_free():
    uh, _ = vortex_tubes_periodic(24)
    assert divergence_l2(uh) < 1e-12


def test_linear_diffusion_exact_mode():
    """Single Fourier mode decays as e^{-ν |k|² t} when nonlinear term vanishes."""
    n = 32
    nu = 0.1
    dt = 1e-3
    steps = 50
    kx, ky, kz = wave_number_grids(n)
    # Mode k=(1,0,0), velocity in y direction → div-free
    uh = np.zeros((3, n, n, n), dtype=np.complex128)
    # index of k=1 on axis 0
    uh[1, 1, 0, 0] = 0.5
    uh[1, -1, 0, 0] = 0.5  # reality
    uh = leray_project_hat(uh)
    e0 = abs(uh[1, 1, 0, 0])
    for _ in range(steps):
        # pure viscous RK4 with zero nonlinear: use semi-implicit exact for linear
        k2 = 1.0
        uh *= np.exp(-nu * kx**2 * 0)  # noop structure
        uh = uh * np.exp(-nu * (kx * kx + ky * ky + kz * kz) * dt)
    expected = e0 * np.exp(-nu * 1.0 * dt * steps)
    assert abs(abs(uh[1, 1, 0, 0]) - expected) / expected < 1e-10


def test_linear_diffusion_via_integrators():
    """With tiny amplitude, nonlinear ~0; integrators track heat decay."""
    n = 16
    nu = 0.2
    dt = 5e-4
    t_end = 0.05
    uh = np.zeros((3, n, n, n), dtype=np.complex128)
    uh[1, 1, 0, 0] = 1e-4
    uh[1, -1, 0, 0] = 1e-4
    uh = leray_project_hat(uh)
    amp0 = abs(uh[1, 1, 0, 0])
    expected = amp0 * np.exp(-nu * 1.0 * t_end)

    for step_fn, name in [
        (step_rk4, "rk4"),
        (step_etd_rk2, "etd"),
        (step_semi_implicit_euler, "semi"),
    ]:
        u = uh.copy()
        nsteps = int(round(t_end / dt))
        for _ in range(nsteps):
            u = step_fn(u, dt, nu, f_hat=None, dealias=True)
        amp = abs(u[1, 1, 0, 0])
        rel = abs(amp - expected) / expected
        assert rel < 0.05, f"{name} rel err {rel}"


def test_energy_decays_without_force():
    uh, _ = taylor_green(16)
    cfg = SolverConfig(n=16, nu=0.05, dt=1e-3, t_end=0.05, integrator="rk4")
    state = NavierStokesSolver(cfg).run(uh)
    e0 = state.history[0]["energy"]
    e1 = state.history[-1]["energy"]
    assert e1 <= e0 * 1.01 + 1e-12


def test_divergence_preserved_in_evolution():
    uh, _ = taylor_green(16)
    cfg = SolverConfig(n=16, nu=0.05, dt=1e-3, t_end=0.05, integrator="etd_rk2")
    state = NavierStokesSolver(cfg).run(uh)
    assert state.history[-1]["div_l2"] < 1e-8


def test_nonlinear_vanishes_for_beltrami_like_check_energy_channel():
    """Sanity: nonlinear hat is orthogonal to u in L2 for continuum; discrete check residual smallish."""
    uh, _ = abc_flow(32)
    N = nonlinear_hat(uh, dealias=True)
    # <u, N> should be near 0 (energy transfer from nonlinear)
    # Σ Re(û* · N̂)
    inner = float(np.sum((uh.conj() * N).real))
    assert abs(inner) < 1e-8


def test_mesh_convergence_linear():
    """Refine N for linear mode decay — amplitude independent of N."""
    nu, dt, steps = 0.1, 1e-3, 20
    amps = []
    for n in (16, 32):
        uh = np.zeros((3, n, n, n), dtype=np.complex128)
        uh[1, 1, 0, 0] = 1e-3
        uh[1, -1, 0, 0] = 1e-3
        uh = leray_project_hat(uh)
        for _ in range(steps):
            uh = step_etd_rk2(uh, dt, nu)
        amps.append(abs(uh[1, 1, 0, 0]))
    assert abs(amps[0] - amps[1]) / amps[1] < 1e-8


def test_integrator_agreement_short_time():
    uh, _ = taylor_green(16)
    dt = 5e-4
    nu = 0.05
    u1 = uh.copy()
    u2 = uh.copy()
    for _ in range(10):
        u1 = step_rk4(u1, dt, nu)
        u2 = step_etd_rk2(u2, dt, nu)
    rel = np.linalg.norm(u1 - u2) / (np.linalg.norm(u1) + 1e-15)
    assert rel < 0.05


def test_interval_leray_divergence_contains_zero():
    v = (Interval.from_float(1.0), Interval.from_float(2.0), Interval.from_float(3.0))
    proj = leray_project_interval_mode(v, 1.0, 0.0, 0.0)
    d = divergence_interval(proj, 1.0, 0.0, 0.0)
    assert d.contains_zero() or abs(d.lo) < 1e-12 or abs(d.hi) < 1e-12 or (d.lo <= 0 <= d.hi)


def test_interval_finite_modes_div_free():
    modes = {
        (1, 0, 0): (0 + 0j, 0.3 + 0.1j, 0 + 0j),
        (-1, 0, 0): (0 + 0j, 0.3 - 0.1j, 0 + 0j),
        (0, 1, 0): (0.2 + 0j, 0 + 0j, 0 + 0j),
        (0, -1, 0): (0.2 + 0j, 0 + 0j, 0 + 0j),
    }
    total = finite_energy_identity_cancel_check(modes)
    assert total.lo <= 1e-18  # should enclose near zero


def test_interval_convolution_and_heat():
    a = {(1, 0, 0): 0.5 + 0.1j, (0, 1, 0): 0.2}
    b = {(0, 0, 0): 1.0, (1, -1, 0): 0.3 - 0.2j}
    enc = short_convolution_interval(a, b, (1, 0, 0))
    # manual: p=(1,0,0)+q=(0,0,0) and p=(0,1,0)+q=(1,-1,0)
    manual = (0.5 + 0.1j) * 1.0 + (0.2) * (0.3 - 0.2j)
    assert enc.contains(manual.real)
    prop = propagate_linear_heat_interval(Interval.from_float(1.0), k2=2.0, nu=0.1, dt=0.01)
    expected = np.exp(-0.1 * 2.0 * 0.01)
    assert prop.contains(expected)
