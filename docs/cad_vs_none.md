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
- `--only-compare`: reuse existing merged roots and regenerate paper outputs + compare plots only.

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

Core/tail plot files (per mode) for paper readability:

- `fig_edep_ccd_core.pdf`, `fig_edep_ccd_tail.pdf`
- `fig_dEdx_core.pdf`, `fig_dEdx_tail.pdf`
- `fig_charge_vs_coszen.pdf` (core y-range) and `fig_charge_vs_coszen_tail.pdf`

Compare folder highlights (`paper_outputs/<tag>/compare/`):

- `compare_hit_efficiency_vs_coszen.pdf`
- `compare_through_fraction_per_thrown_vs_coszen.pdf`
- `compare_through_fraction_of_hits_vs_coszen.pdf`
- `compare_edep_core.pdf`, `compare_edep_tail.pdf`
- `compare_dedx_core.pdf`, `compare_dedx_tail.pdf`
- `compare_lcos_distribution.pdf`
- `compare_charge_hits_core.pdf`, `compare_charge_hits_tail.pdf`, plus throughgoing variants when available

Regenerate comparisons only (reuse existing paper outputs, both layouts supported):

```bash
python analysis/compare_modes.py \
  --tag cad_vs_none_prod_v1 \
  --paper-none paper_outputs/cad_vs_none_prod_v1_none \
  --paper-cad paper_outputs/cad_vs_none_prod_v1_cad \
  --out paper_outputs/cad_vs_none_prod_v1/compare
```

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

Bundle everything for sharing:

```bash
python analysis/package_paper_outputs.py --tag cad_vs_none_prod_v1
```

## Dependencies

- `uproot` (already required by `analysis/make_paper_outputs.py`)
- `PyYAML` (for the cluster entrypoint config loader)
