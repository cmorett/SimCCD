# Paper-ready figures and tables

## Dependencies

- Python 3 with `uproot`, `numpy`, `matplotlib` (and `pandas` optional).
- Uses `pixelization/pixelize_helpers.py` (standalone) for offline track images.

## Run the analysis on any B02ntuples ROOT file

```
python analysis/make_paper_outputs.py \
  --input build/vs2022/main/Release/B02ntuples.root \
  --output paper_outputs --tag demo \
  --examples 200 --dist-events 50000
```

Outputs (under `paper_outputs/<tag>/`):

- Figures:
  - `fig_energy_spectrum.pdf`
  - `fig_coszenith_down.pdf`
  - `fig_xyImpact.pdf`
  - `fig_z0.pdf`
  - `fig_edep_ccd.pdf`
  - `fig_trackLen_ccd.pdf`
  - `fig_edep_vs_trackLen.pdf`
  - `fig_costheta_vs_trackLen.pdf`
  - `fig_pixelized_examples.pdf` (plus per-event PNGs in `images/`)
  - Pixel metrics distributions: `fig_cluster_size.pdf`, `fig_cluster_charge.pdf`, `fig_sigma_x.pdf`, `fig_sigma_y.pdf`, `fig_length_pix.pdf`
- Tables:
  - `tables/validation_summary.csv` (energy, clamp counts, CCD stats, impact uniformity)
  - `tables/pixel_metrics.csv` (if pixelization samples exist)
  - `tables/units_and_conventions.txt`
- Config: `run_config.json`

Key assumptions (also written to `units_and_conventions.txt` and `run_config.json`):

- `EevtPri` in GeV; `thetaPri`/`phiPri` are Geant4 angles; downward cos(zenith) = `-cos(thetaPri)`.
- Impact-plane coordinates `muonXImp/muonYImp/muonZImp` in cm.
- CCD summaries:
  - `EdepCCD` in GeV (total deposited energy in sensitive silicon).
  - `trackLenCCD` chord length between first/last CCD step (cm).
  - `dirX/Y/Z` primary direction unit vector.
- Pixelization model defaults:
  - Pixel size 15 µm, CCD thickness 725 µm.
  - 3.7 eV per electron; Gaussian transverse diffusion with a simple depth-dependent sigma.

## One-command end-to-end (optional)

If you ran the 1M script with `--make-paper-outputs`, it will call the analysis automatically on the merged ROOT file (if `hadd` is available). Otherwise, point `--input` to the merged or single-chunk ROOT file.
