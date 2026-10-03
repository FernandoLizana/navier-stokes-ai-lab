# Contributing

Thanks for interest in **Navier–Stokes AI Lab** (historical: NS-MRL / bastardus2).

## Setup

```bash
python3 -m venv .venv
source .venv/bin/activate   # Windows: .\.venv\Scripts\Activate.ps1
pip install -e ".[dev]"
python -m ns_exploration.demo
pytest python/ns_exploration/tests/test_sprint01.py python/ns_exploration/tests/test_demo.py -q
```

Heavy MPFR / certificate campaigns and Julia/Lean stacks are optional. See `python/requirements-c0008-repair.txt` and `docs/PUBLICATION_NOTES.md`.

## Reporting bugs

Include: OS, Python version, exact command, expected vs actual behavior, and whether you used the demo or a heavier path. Do not paste secrets or private paths.

## Scientific claims

When proposing a mathematical or numerical claim:

1. State the **domain** (e.g. Galerkin N≤24 on \(\mathbb{T}^3\), not continuum PDE).
2. Give an **evidence label** (`docs/SCIENTIFIC_INTEGRITY.md`).
3. Point to **artifactsable artifacts** (test, script, certificate + verifier).
4. Do not equate verifier PASS with continuum proof.

## Pull requests

- Prefer small, reviewable changes.
- Do not “fix” failing scientific checks by weakening thresholds without documentation.
- Do not delete historical failed certificates or audit reports to clean the narrative.
- Keep Route B (unforced) and forced/blow-up tracks clearly separated in claims.

## License

Contributions are under the MIT license (`LICENSE`).
