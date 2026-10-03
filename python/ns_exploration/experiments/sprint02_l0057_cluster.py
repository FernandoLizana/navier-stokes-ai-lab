"""Sprint 02: L-0057 greedy cluster SOS job pack."""

from __future__ import annotations

from ns_exploration.conjectures.l0057_greedy_cluster_sos import (
    lemma_l0057,
    save_lemma_l0057,
    write_cluster_runbook,
    write_slurm_template,
)


def main() -> None:
    l57 = lemma_l0057()
    save_lemma_l0057(l57)
    runbook = write_cluster_runbook(l57)
    slurm = write_slurm_template(l57)
    print(f"L-0057: {l57.lemma_id} status={l57.status}")
    print(f"  D={l57.greedy_D} extrap_SOS={l57.C_ub_sos_extrap:.3f}")
    print(f"  runbook={runbook}")
    print(f"  slurm={slurm}")


if __name__ == "__main__":
    main()
