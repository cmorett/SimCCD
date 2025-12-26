# Paper-ready figures and tables

## Dependencies

- Python 3 with `uproot`, `numpy`, `matplotlib` (and `pandas` optional).
- Uses `pixelization/pixelize_helpers.py` (standalone) for offline track images.

## Run the analysis on any B02ntuples ROOT file

```
python analysis/make_paper_outputs.py \
  --input build/vs2022/main/Release/B02ntuples.root \
  --output paper_outputs --tag demo \
  --examples 24 --dist-events 50000 \
  --pixel-size-microns 15 --thickness-microns 725 \
  --canvas-mode adaptive --margin-pix 24 \
  --quality-only
```

Key CLI knobs: `--canvas-mode adaptive|fixed`, `--canvas-size` (when fixed), `--quality-only` (plots use non-truncated events), `--pixel-size-microns`, `--thickness-microns`, `--ev-per-electron`, `--margin-pix`, `--seed`, `--examples`, `--dist-events`.

Outputs (under `paper_outputs/<tag>/`):

- Figures:
  - Source/geometry: `fig_energy_spectrum.pdf`, `fig_coszenith_down.pdf`, `fig_coszen_all_vs_hits.pdf`, `fig_hit_efficiency_vs_coszen.pdf`, `fig_xyImpact.pdf`, `fig_z0.pdf`, `fig_xy0.pdf`
  - Energy validation: `validate_energy_hits_logy.pdf`, `validate_energy_hits_zoom.pdf`
  - CCD hits: `fig_edep_ccd.pdf`, `fig_trackLen_ccd.pdf`, `fig_edep_vs_trackLen.pdf`, `fig_edep_vs_len_hits.pdf`, `fig_costheta_vs_trackLen.pdf`
  - dE/dx: `fig_dEdx_hits.pdf`, `fig_dEdx.pdf` (zoom), `fig_dEdx_tail.pdf` (log)
  - Pixel examples: `fig_pixelized_examples.pdf` + per-event PNGs in `images/` (quality events with scale bars and annotations)
  - Pixel metrics (quality overlaid on "all" unless `--quality-only`): `fig_cluster_size.pdf`, `fig_cluster_charge.pdf`, `fig_sigma_trans.pdf`, `fig_sigma_long.pdf`, `fig_elongation.pdf`, `fig_sigma_x_axis.pdf`, `fig_sigma_y_axis.pdf`, `fig_width_vs_coszen.pdf`, `fig_charge_vs_coszen.pdf`
  - Length diagnostics: `fig_length_pix_geom.pdf`, `fig_length_pix_img.pdf`, `fig_length_pix.pdf` (PCA vs expected 2D), `fig_length_pix_legacy.pdf` (legacy 3D vs PCA)
- Tables:
  - `tables/validation_summary.csv` (energy, clamp counts, CCD stats, dE/dx, impact uniformity, truncation fraction, sigma/length stats)
  - `tables/pixel_metrics_all.csv` and `tables/pixel_metrics_quality.csv` (pixel-level metrics with truncation flag)
  - `tables/units_and_conventions.txt`
- Config: `run_config.json`

Quality cuts applied to pixel metrics: `EdepCCD>0`, `trackLenCCD>0`, and `is_truncated=False` (charge touching image edge). PCA-based sigmas/lengths are the physics-facing metrics; axis-aligned sigmas are kept as diagnostics only.

Quick regression sanity (optional):

```
python analysis/regression_check.py --summary paper_outputs/<tag>/tables/validation_summary.csv
```

## One-command end-to-end (optional)

If you ran the 1M script with `--make-paper-outputs`, it will call the analysis automatically on the merged ROOT file (if `hadd` is available). Otherwise, point `--input` to the merged or single-chunk ROOT file.
