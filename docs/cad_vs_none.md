# CAD vs none pipeline

This workflow runs paired simulations (no CAD import vs tessellated CAD), produces
paper-ready outputs, generates CAD/none comparison plots with uncertainties, and
writes a markdown summary for a paper appendix.

For xook/SLURM step-by-step instructions, see `docs/cluster_run_cad_vs_none.md`.

## Entry point

Run the cluster entrypoint from the repo root:

```bash
python3 analysis/cluster_run_cad_vs_none_prod.py \
  --config configs/cad_vs_none_prod.yaml
```

Key flags:

- `--sanity-only`: run only the sanity gate (10k/10k).
- `--skip-sanity`: skip the sanity gate and go straight to production.
- `--resume`: skip completed batches and existing analysis outputs.
- `--dry-run`: print commands only (no files created).

## Sanity gate (automatic)

Before production, the pipeline runs:

- 10k thrown (none) + 10k thrown (cad)
- batch size = 5000, max workers = 4
- paper outputs for both modes
- comparison plots in `paper_outputs/<tag>/sanity/compare/`
- summary in `docs/cad_vs_none_summary_<tag>_sanity.md`

The sanity gate fails (and production is blocked) if any of:

- merged ROOT files are missing/empty
- cutflow hits/throughgoing are zero
- Lcos median deviates by >2% from thickness
- required comparison plots are missing

## Production run

After sanity passes, production proceeds:

- no CAD (cad-mode none): 1,000,000 thrown
- CAD tessellated: 500,000 thrown
- merged ROOTs in `outputs/<tag>/none/merged.root` and `outputs/<tag>/cad/merged.root`
- paper outputs in `paper_outputs/<tag>/none/` and `paper_outputs/<tag>/cad/`
- comparison plots in `paper_outputs/<tag>/compare/`
- summary in `docs/cad_vs_none_summary_<tag>.md`

## Output layout

```
outputs/<tag>/
  none/
    batches/
    logs/
    macros/
    merged.root
  cad/
    batches/
    logs/
    macros/
    merged.root
  run_metadata.json

paper_outputs/<tag>/
  none/
  cad/
  compare/
  sanity/
    none/
    cad/
    compare/

docs/
  cad_vs_none_summary_<tag>_sanity.md
  cad_vs_none_summary_<tag>.md
```

## Dependencies

- `uproot` (already required by `analysis/make_paper_outputs.py`)
- `PyYAML` (for the cluster entrypoint config loader)
