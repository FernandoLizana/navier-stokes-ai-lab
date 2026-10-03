"""Band pilot: compare Route A integral frobenius vs best (1-inf)."""
from __future__ import annotations

import json

from gmpy2 import mpfr

from ns_exploration.terminal_weighted.certify_shell_integral import shell_integral_certificate

m = json.load(open("experiments/terminal_weighted/shell_manifest.json", encoding="utf-8"))
bf = [{**b, "C_term_hi": b["C_term_frobenius_hi"], "method": "frobenius_mpfr_up"} for b in m["blocks"]]
mf = dict(m)
mf["blocks"] = bf
mf["bound_method"] = "frobenius"
cf = shell_integral_certificate(24, prec=128, manifest=mf)
cb = shell_integral_certificate(24, prec=128, manifest=m)
ratio = float(mpfr(cb["integral_hi"])) / float(mpfr(cf["integral_hi"]))
print("band Route A frobenius I_hi", cf["integral_hi"])
print("band Route A best I_hi", cb["integral_hi"])
print("ratio best/fro", ratio)
