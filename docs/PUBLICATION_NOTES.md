# Publication notes

**Default policy: keep the GitHub repository PRIVATE** until the maintainer explicitly switches visibility to public.

## Private → public checklist

1. Confirm no secrets, tokens, or personal absolute paths remain (`rg -i 'fliza|api_key|ghp_|sk-'` or equivalent).
2. Confirm large binaries policy (LFS / release assets) for items like `reports/l0049_full_M.npz`.
3. Update About + topics on GitHub.
4. Only then: GitHub → Settings → Change repository visibility → Public.

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
