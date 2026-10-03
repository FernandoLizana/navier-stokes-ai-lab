# Publication notes

**Status:** repository is **PUBLIC** (`https://github.com/FernandoLizana/navier-stokes-ai-lab`).

## Pre-public checklist (completed locally before visibility change)

1. Scanned tracked tree for tokens (`ghp_` / `gho_` / `github_pat_` / private keys / AWS-style keys): only documentation mentions of scan patterns.
2. No personal absolute paths (`Users\\…`) and no personal username leftovers in tracked content.
3. No `.env` / credential / `.pem` files tracked.
4. Workflows use `actions/checkout` + pip/pytest only; no `secrets.*` usage.
5. Large regenerable dump `reports/l0049_full_M.npz` kept out of git; remaining large files are scientific JSON (~5–9 MB).

## Suggested GitHub blurb

**Description:**  
AI-assisted numerical experiments, verification tooling, and lessons from an ambitious Navier–Stokes exploration. Not a Millennium Prize solution.

**Topics (suggested):** `navier-stokes`, `pseudospectral`, `scientific-computing`, `verification`, `lean4`, `open-science`, `ai-assisted-research`

**Slug:** `navier-stokes-ai-lab` (rename remote when you choose; keep `ns_exploration` imports).

## Large / duplicate artifacts

| Item | Role | Suggestion |
|------|------|------------|
| `reports/l0049_full_M.npz` | Large binary (~hundreds of MB) | Git LFS or release asset; document checksum if published |
| `python/experiments/terminal_weighted/shell_one_inf_ck/` | Intra-shell checkpoints | Local / `compute_pack` only; gitignored pattern |
| `compute_pack/data/` | Portable compute state | Transfer folder; not required for clone+demo |
| `artifacts/C0008_HANDOFF_PACK/` | Dated snapshot | Keep as archive; point newcomers to live `python/` |
| Root `*.rar` dumps | Offline transfer | Keep out of git (`.gitignore`) |

`.gitignore` does **not** remove already-tracked large files. If something huge is already in history, plan a maintainer-approved history cleanup or LFS migration—do not force-push without agreement.

## Secrets scan

Before publish: search for `.env`, API keys, private absolute paths in notebooks, and credential JSON. Do not copy secret values into issues or this doc.

## Third-party notices

MIT license retained (`LICENSE`). Dependencies (NumPy, SciPy, Matplotlib, optional gmpy2, Lean mathlib if used) keep their own licenses—cite when redistributing binaries.

## CI note

`demo-ci.yml` is configured in-tree. Until it runs on GitHub Actions for your fork, treat it as **pending remote execution**. Local demo + pytest smoke were exercised in the repositioning review.

## What not to claim in the About box

- Solved / nearly solved Clay NS.
- Affiliation with Clay Math Institute or OpenAI.
- Percentage progress toward the Millennium problem from finite I_hi gaps.
