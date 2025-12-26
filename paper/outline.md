# SimCCD Tier-B Paper Outline (scaffold)

## 1. Abstract
- One-paragraph summary of the Tier-B flux sampling, CCD hit efficiency, and validation against geometric expectations.

## 2. Introduction
- Motivation for realistic flux sampling vs. forced footprint.
- Brief overview of CCD geometry and readout goals.

## 3. Simulation and Geometry
- Geant4 setup, primitive geometry, and 725 um CCD thickness.
- Scoring volume definition and provenance tracking.

## 4. Generator and Event Weights
- Modified Gaisser (Guan 2015) flux sampling.
- Source plane configuration and event livetime weights.

## 5. Pixelization and Reconstruction
- Offline pixelization, PCA metrics, and quality selection.
- Definitions of L3D, L2D_expected, and L2D_image.

## 6. Validation Results
- Length consistency (projected vs. image).
- dE/dx and energy spectra for hits.
- Angular dependence of efficiency and widths.

## 7. Conclusions / Future Work
- Summary of validation and readiness for CAD geometry.
- Planned improvements and targeted sampling mode.

## Figures (paper outputs)
- `fig_energy_spectrum.pdf`: primary energy spectrum (all thrown).
- `validate_energy_hits_logy.pdf`: primary energy, all vs hits (log y).
- `validate_energy_hits_zoom.pdf`: primary energy zoom with fine binning.
- `fig_coszenith_down.pdf`: downward cos(zenith) distribution.
- `fig_coszen_all_vs_hits.pdf`: cos(zenith) all vs hits (normalized).
- `fig_hit_efficiency_vs_coszen.pdf`: hit efficiency vs cos(zenith).
- `fig_xyImpact.pdf`: impact plane distribution (all thrown).
- `fig_edep_ccd.pdf`: CCD deposited energy (hits).
- `fig_trackLen_ccd.pdf`: CCD track length (hits).
- `fig_edep_vs_len_hits.pdf`: Edep vs L3D with mean dE/dx guide.
- `fig_dEdx_hits.pdf`: dE/dx distribution for hits.
- `fig_length_pix.pdf`: L2D_image vs L2D_expected + ratio.
- `fig_length_pix_legacy.pdf`: legacy 3D vs 2D length comparison.
- `fig_width_vs_coszen.pdf`: cluster width vs cos(zenith).
- `fig_charge_vs_coszen.pdf`: cluster charge vs cos(zenith).
