"""Tests for adjoint, gate battery pieces, and conjecture C-0001."""

from __future__ import annotations

from ns_exploration.conjectures.c0001_enstrophy_bound import evaluate_conjecture
from ns_exploration.initial_conditions.generators import random_div_free, taylor_green
from ns_exploration.optimization.adjoint import (
    enstrophy_gradient_hat,
    nonlinear_jvp,
    nonlinear_vjp,
    real_inner,
    verify_adjoint_vs_fd,
)
from ns_exploration.optimization.candidate_pipeline import run_gate_battery
from ns_exploration.spectral.leray import leray_project_hat


def test_nonlinear_vjp_matches_jvp_inner_product():
    uh, _ = taylor_green(16)
    v, _ = random_div_free(16, seed=3, energy_target=0.1)
    lam, _ = random_div_free(16, seed=4, energy_target=0.1)
    v = leray_project_hat(v)
    lam = leray_project_hat(lam)
    jvp = nonlinear_jvp(uh, v, eps=1e-6)
    vjp = nonlinear_vjp(uh, lam)
    lhs = real_inner(lam, jvp)
    rhs = real_inner(vjp, v)
    rel = abs(lhs - rhs) / (abs(lhs) + abs(rhs) + 1e-30)
    assert rel < 0.05, f"VJP/JVP mismatch rel={rel} lhs={lhs} rhs={rhs}"


def test_enstrophy_gradient_directional():
    uh, _ = taylor_green(16)
    d, _ = random_div_free(16, seed=5, energy_target=0.05)
    d = leray_project_hat(d)
    g = enstrophy_gradient_hat(uh)
    eps = 1e-6
    from ns_exploration.diagnostics.metrics import compute_diagnostics

    jp = compute_diagnostics(uh + eps * d, 0.1).enstrophy
    jm = compute_diagnostics(uh - eps * d, 0.1).enstrophy
    fd = (jp - jm) / (2 * eps)
    adj = real_inner(g, d)
    rel = abs(fd - adj) / (abs(fd) + abs(adj) + 1e-30)
    assert rel < 1e-4


def test_adjoint_vs_fd_short_time():
    u, _ = random_div_free(12, seed=0, energy_target=0.5)
    d, _ = random_div_free(12, seed=17, energy_target=1.0)
    rep = verify_adjoint_vs_fd(u, d, nu=0.1, dt=1e-3, t_end=0.008, tol=0.2)
    assert rep.agreed, f"rel_err={rep.relative_error}"


def test_gate_battery_small():
    rows = run_gate_battery(seeds=[0, 1], n_low=12, gate_tol=0.5, t_end=0.015)
    assert len(rows) == 2
    assert all("gate_passed" in r for r in rows)


def test_conjecture_c0001_evaluates():
    c = evaluate_conjecture(M=1e6, n=12)
    assert c.id == "C-0001"
    assert c.numerical_support_max is not None
    assert c.clay_implication.startswith("None")
    c2 = evaluate_conjecture(M=1e-8, n=12)
    assert c2.refuted
