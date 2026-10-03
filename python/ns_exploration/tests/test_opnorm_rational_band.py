"""Production opnorm_certified rational Leray vs float reference."""

from __future__ import annotations

from gmpy2 import mpfr

from ns_exploration.terminal_weighted.opnorm_certified import certified_shell_blocks


def test_rational_band123_encloses_float():
    radii = (1, 2, 3)
    prec = 128
    fl = certified_shell_blocks(24, radii, prec=prec, workers=1, rational_coefficients=False)
    rat = certified_shell_blocks(24, radii, prec=prec, workers=1, rational_coefficients=True)
    assert rat.get("rational_coefficients") is True
    for bf, br in zip(fl["blocks"], rat["blocks"]):
        assert mpfr(br["C_term_hi"]) >= mpfr(bf["C_term_hi"]) * mpfr("0.999999")
    fl_sum = sum(mpfr(b["C_term_hi"]) for b in fl["blocks"])
    rat_sum = sum(mpfr(b["C_term_hi"]) for b in rat["blocks"])
    assert float(rat_sum) >= float(fl_sum) * (1 - 1e-9)
