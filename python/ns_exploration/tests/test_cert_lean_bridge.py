"""
Cross-check that the Python N5 certificate constants match the values proved
exactly in the Lean bridge module `NSGalerkin.Certificates`
(shell_k3_cap = 18450, full_N12_cap = 51450, and Ω_eq = 150·M for the fixed
(E0,ν)=(1/2,1/10) regime).

This does NOT run Lean; it pins the Python side to the Lean-verified integers,
so a drift on either side breaks a test. Lean compilation is exercised by the
experiment module.
"""

from __future__ import annotations

from ns_exploration.validation.l0003_certificate import (
    build_l0003_certificate,
    build_l0003_full_dealias_certificate,
)

# Values proved in lean/NSGalerkin/NSGalerkin/Certificates.lean (native_decide + decide).
LEAN_SHELL_K3_CAP = 18450
LEAN_FULL_N12_CAP = 51450
LEAN_FULL_N16_CAP = 199650


def test_python_shell_matches_lean():
    cert = build_l0003_certificate(n=12, k_shell=3, E0=0.5, nu=0.1)
    assert abs(cert.Omega_bound_hi - LEAN_SHELL_K3_CAP) < 1e-6
    # cleared-denominator relation Ω_eq = 150·M
    assert abs(cert.Omega_eq_hi - 150 * cert.M) < 1e-6


def test_python_full_matches_lean():
    cert = build_l0003_full_dealias_certificate(n=12, E0=0.5, nu=0.1)
    assert abs(cert.Omega_bound_hi - LEAN_FULL_N12_CAP) < 1e-6
    assert abs(cert.Omega_eq_hi - 150 * cert.M) < 1e-6


def test_python_full_n16_matches_lean():
    cert = build_l0003_full_dealias_certificate(n=16, E0=0.5, nu=0.1)
    assert cert.K_squared_int == 75
    assert cert.M == 1331
    assert abs(cert.Omega_bound_hi - LEAN_FULL_N16_CAP) < 1e-6
    assert abs(cert.Omega_eq_hi - 150 * cert.M) < 1e-6


def test_growth_n12_to_n16():
    """Documents ∼M/ν² growth: N16 cap / N12 cap = M16/M12."""
    c12 = build_l0003_full_dealias_certificate(n=12)
    c16 = build_l0003_full_dealias_certificate(n=16)
    assert c16.M > c12.M
    assert c16.Omega_bound_hi > c12.Omega_bound_hi
    # exact ratio of cleared caps equals ratio of M (same E0,ν)
    assert abs(c16.Omega_eq_hi / c12.Omega_eq_hi - c16.M / c12.M) < 1e-9


def test_cleared_relation_general():
    """Ω_eq = 6 M E0²/ν² reduces to 150 M exactly when E0=1/2, ν=1/10."""
    E0, nu = 0.5, 0.1
    for M in (33, 123, 343, 1331):
        omega_eq = 6 * M * E0**2 / nu**2
        assert abs(omega_eq - 150 * M) < 1e-6


def test_n16_certificate_verifies():
    from ns_exploration.validation.l0003_certificate_verify import verify_l0003_certificate

    cert = build_l0003_full_dealias_certificate(n=16)
    ok, checks = verify_l0003_certificate(cert.as_dict())
    assert ok
    assert all(checks.values())


def test_l0018_envelope_matches_lean_omega0():
    """L-0018 / C-0003 cap = 75/2, already proved as full_N16_omega0 in Lean."""
    from ns_exploration.conjectures.l0018_envelope import lemma_l0018
    from ns_exploration.validation.l0018_certificate import build_l0018_envelope_certificate

    LEAN_N16_OMEGA0 = 75 / 2  # theorem full_N16_omega0 / l0018_envelope_N16
    b = lemma_l0018()
    assert b.K2_worst == 75
    assert abs(b.Omega_cap - LEAN_N16_OMEGA0) < 1e-15
    cert = build_l0018_envelope_certificate(R=LEAN_N16_OMEGA0)
    assert cert.closes
    assert abs(cert.Omega_cap_hi - LEAN_N16_OMEGA0) < 1e-15
