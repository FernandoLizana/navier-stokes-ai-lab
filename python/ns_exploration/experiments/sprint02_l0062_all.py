"""Sprint: L-0062 ladder + L-0063 greedy export + refresh maps."""

from __future__ import annotations

from ns_exploration.conjectures.l0062_cosmo_ladder_n5 import lemma_l0062, save_lemma_l0062
from ns_exploration.conjectures.l0063_greedy_cluster_export import (
    lemma_l0063,
    run_greedy_export,
    save_lemma_l0063,
)
from ns_exploration.conjectures.l0057_greedy_cluster_sos import (
    lemma_l0057,
    save_lemma_l0057,
    write_cluster_runbook,
    write_slurm_template,
)
from ns_exploration.validation.l0062_certificate import (
    build_l0062_certificates,
    save_l0062_certificates,
)


def main() -> None:
    # L-0062: COSMO bands {1..9}, {1..10}, {1..12}
    l62 = lemma_l0062()
    save_lemma_l0062(l62)
    certs = build_l0062_certificates()
    paths = save_l0062_certificates(certs)
    print(f"L-0062: {len(l62.rows or [])} bands validated")
    for r in l62.rows or []:
        print(f"  {{1..{r.band_radii[-1]}}} D={r.D} C_ub_hi={r.C_ub_hi:.4f} ratio={r.ratio_sos_over_shor:.3f}")
    for p in paths:
        print(f"  cert={p.name}")

    # L-0063: greedy export if missing
    l63 = lemma_l0063()
    if not l63.export_complete:
        print("L-0063: running greedy sparse export (may take 30+ min)...")
        exp = run_greedy_export()
        print(f"  export returncode={exp['returncode']}")
        l63 = lemma_l0063()
    save_lemma_l0063(l63)
    print(f"L-0063: status={l63.status} export_complete={l63.export_complete}")

    # Refresh L-0057 cluster pack
    l57 = lemma_l0057()
    save_lemma_l0057(l57)
    write_cluster_runbook(l57)
    write_slurm_template(l57)
    print(f"L-0057: cluster_ready D={l57.greedy_D} extrap_SOS={l57.C_ub_sos_extrap:.3f}")


if __name__ == "__main__":
    main()
