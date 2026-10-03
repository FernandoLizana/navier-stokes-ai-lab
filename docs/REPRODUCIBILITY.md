# Reproducibility

Every validated result must ship: source, compiler versions, dependency pins, architecture, parameters, hashes, data, intervals, logs, certificate, minimal verifier, container.

Verifier ≪ generator. Do not require GPU, NN, dashboard, or network to verify.

## Public demo (entry point)

From a clean checkout with the root package installed (`pip install -e ".[dev]"`):

```bash
python -m ns_exploration.demo --out-dir demos/out/latest
```

Outputs: `diagnostics.csv`, `config.json`, `summary.json`, `diagnostics.png`.  
Evidence level **N2** (exploratory float). Not a continuum proof.

See also `README.md` quick start and `docs/PUBLICATION_NOTES.md` for large artifacts.
