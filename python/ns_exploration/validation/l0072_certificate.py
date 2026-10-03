"""L-0072 route-refutation certificate (terminal full-dealias, C-0008)."""

from __future__ import annotations

import json
from pathlib import Path

from ns_exploration.terminal_weighted.constants import FROZEN

DEFAULT_CERT = Path("certificates/CERT-L0072-C0008-terminal-full-dealias.json")
DEFAULT_MANIFEST = Path("experiments/terminal_weighted/shell_manifest_best_full.json")
REGISTRY = Path("conjectures/active/L-0072.json")


def verify_l0072_route_refuted(
    cert_path: str | Path = DEFAULT_CERT,
    manifest_path: str | Path = DEFAULT_MANIFEST,
) -> tuple[bool, list[str]]:
    """
    Expected outcome: verifier FAIL (route cannot close C-0008).

    Returns (route_refuted, messages). route_refuted=True when certificate
    is well-formed and rigorously fails the closure checks.
    """
    import sys

    root = Path(__file__).resolve().parents[3]
    if str(root) not in sys.path:
        sys.path.insert(0, str(root))
    from verify_c0008_terminal_certificate import verify_certificate

    cert_path = Path(cert_path)
    manifest_path = Path(manifest_path)
    msgs: list[str] = []

    if not cert_path.is_file():
        return False, [f"Missing certificate: {cert_path}"]
    if not manifest_path.is_file():
        return False, [f"Missing manifest: {manifest_path}"]

    cert = json.loads(cert_path.read_text(encoding="utf-8"))
    ok, vmsgs = verify_certificate(cert_path, manifest_path=manifest_path)
    msgs.extend(vmsgs)

    route_refuted = (
        not ok
        and cert.get("closes_strict") is False
        and float(cert.get("parameters", {}).get("M_target", 0)) == FROZEN.c0008_target_M
        and any("I_term_hi" in m or "omega_T_hi" in m for m in vmsgs)
    )
    if route_refuted:
        msgs.append("ROUTE_REFUTED: terminal shell integral cannot close C-0008")
    else:
        msgs.append("UNEXPECTED: expected FAIL closure checks for L-0072 route")
    return route_refuted, msgs


def build_l0072_registry_snapshot(
    cert_path: str | Path = DEFAULT_CERT,
    manifest_path: str | Path = DEFAULT_MANIFEST,
) -> dict:
    cert = json.loads(Path(cert_path).read_text(encoding="utf-8"))
    manifest = json.loads(Path(manifest_path).read_text(encoding="utf-8"))
    intv = cert.get("integral_interval", ["0", "0"])
    omega = cert.get("omega_T_interval", ["0", "0"])
    istar = cert.get("i_star_interval", ["0", "0"])
    return {
        "lemma_id": "L-0072",
        "route": "terminal",
        "status": "route_refuted",
        "evidence_level": "N5",
        "n": 24,
        "parent": "L-0071",
        "claim": "C-0008",
        "phase_e_best_route_D_I_hi": float(intv[1]),
        "phase_e_omega_T_hi": float(omega[1]),
        "I_star_lo": float(istar[0]),
        "verifier": "FAIL",
        "certificate": str(cert_path).replace("\\", "/"),
        "manifest_best": str(manifest_path).replace("\\", "/"),
        "manifest_sha256": manifest.get("manifest_sha256", cert.get("source_hashes", {}).get("shell_manifest")),
        "refutation_type": "route_refuted",
        "refutes_conjecture": False,
        "closes_c0008": False,
        "clay_implication": "None. Route refutation on finite Galerkin relaxation only; not Clay.",
        "notes": cert.get("limitations", ["per-shell triangle"])[0] if cert.get("limitations") else "",
    }


def save_l0072_registry(path: str | Path = REGISTRY) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(build_l0072_registry_snapshot(), indent=2), encoding="utf-8")
    return path
