# Cosmic muon sampling modes

SimCCD now supports two generator modes:

- **forced_footprint** (legacy): back-projects muons through the CCD footprint; nearly every generated muon hits the CCD. Useful for fast pixel-statistics studies.
- **tierB_plane_flux** (new, Tier B): samples muons from a realistic Modified Gaisser flux on a source plane larger than the CCD, producing natural hit/miss fractions and an event weight (effective livetime) per event.

## Flux model (Guan 2015 Modified Gaisser)

We use the sea-level parameterisation from Guan et al. (2015) which extends the standard PDG Gaisser form to low energies and large zenith angles:

- Earth-curvature correction:  
  `cos(theta*) = sqrt( cos^2(theta) + P1^2 + P2 cos^P3(theta) + P4 cos^P5(theta) )`  
  with `P1=0.102573`, `P2=-0.068287`, `P3=0.958633`, `P4=0.0407253`, `P5=0.817285`.
- Modified intensity (units: cm^-2 s^-1 sr^-1 GeV^-1):  
  `dI/dE/dΩ = 0.14 E^{-2.7} ( 1/(1 + 1.1 E cosθ*/115) + 0.054/(1 + 1.1 E cosθ*/850) )`
  multiplied by the low-energy correction `(1 + 3.64 / (E * cosθ*^{1.29}))^{-1}`.

Sampling accounts for the projected rate through a horizontal plane: the joint density is proportional to `I(E,θ) * cosθ * sinθ` with `φ` uniform in `[0, 2π)`.

## Source plane geometry

- The source plane sits at `z = sourcePlaneZ` above the impact plane (default `z=0` at the CCD plane).
- By default (`sourcePlaneAutoSize=false`), `(x0,y0)` is uniform over a manual rectangle `sourcePlaneLx x sourcePlaneLy` (defaults: 10 cm x 10 cm). This already produces a nontrivial miss fraction for a sub-cm^2 CCD.
- Auto-size (`sourcePlaneAutoSize=true`) enlarges the plane to cover the CCD footprint for a chosen `thetaMax`:  
  `Lx = max(sourcePlaneLx, px*cm + 2*z0*tan(thetaMax) + margin)` and same for `Ly`. This is convenient when tilting the flux to larger zenith angles; expect a smaller hit fraction as the plane grows.
- Event weights: for Tier B with flux sampling (`useFixedEnergy=false`), each event stores `muonWeight_s` and `eventLivetime_s = 1 / (flux_integral * plane_area)` so `N` generated events represent `N * eventLivetime_s` seconds of effective exposure. Weights are set to zero when `useFixedEnergy=true`.
- Generator configuration is written into the ntuple (`B02Evts`): `muonModeCode`, `fluxModelCode`, `cfg_sourcePlaneZ_cm`, `cfg_sourcePlaneLx_cm`, `cfg_sourcePlaneLy_cm`, `cfg_thetaMax_deg`, `cfg_EminGeV_eff`, `cfg_EmaxGeV_eff` for reproducibility.

## Overburden and charge ratio (realism knobs)

- Overburden slab (optional): `/sim/overburden/enable true|false`, `/sim/overburden/thickness <value> <unit>`, `/sim/overburden/zTop <value> <unit>`, `/sim/overburden/material <G4_NAME>`. Default is disabled. The slab sits below `zTop` (top surface), so a muon generated on the source plane at `zTop` immediately traverses the thickness.
- Muon charge control: `/generator/muonChargeMode equal|fixedRatio|plusOnly|minusOnly` and `/generator/muonPlusToMinusRatio <r>` (default 1.25). Event ntuple stores `muonPDG` and `muonChargeSign`.
- CCD region controls: `/sim/cuts/ccdGammaCut`, `/sim/cuts/ccdElectronCut`, `/sim/cuts/ccdPositronCut`, and `/sim/ccd/maxStep` (0 disables). CCD thickness is 0.05 cm by default.

## Provenance / reproducibility

- Run-level settings are attached to the output via constant per-event columns (prefix `prov_`): seeds (`prov_seed1/2`, `prov_useTimeSeed`), overburden toggles and dimensions, CCD cuts/step/max thickness, muon charge ratio, and hashed identifiers for git/macro/physics list (`prov_gitHashCode`, `prov_macroHashCode`, `prov_macroPathHash`, `prov_physicsListHash`).
- `analysis/make_paper_outputs.py` surfaces these in `tables/validation_summary.csv` and `run_config.json`. If available, the dedicated `B02RunInfo` tree will also carry the same information (string columns may be disabled on some builds).
- Log line `[RunInfo] git=... macro=... seeds=(... ...)` is printed once per run to make the provenance explicit in batch logs.

## Quick start

Batch run using the provided macro:
```
main/build/Release/b02_executable.exe --no-vis main/macros/run_tierB_cosmic_muons.mac
```

Key UI commands (all under `/generator/`):
- `muonMode tierB_plane_flux` (or `forced_footprint`)
- `fluxModel guan2015|pdg_gaisser`
- `useFixedEnergy true|false` plus `muonEnergyGeV` (fixed) or `EminGeV`/`EmaxGeV` (flux)
- `thetaMax <deg>`
- `sourcePlaneZ <value> cm`
- `sourcePlaneAutoSize true|false`
- `sourcePlaneLx <value> cm`, `sourcePlaneLy <value> cm`, `sourcePlaneMargin <value> cm`

Outputs (ROOT ntuple `B02Evts`) now include Tier B diagnostics: `muonCosTheta`, `muonWeight_s`, `eventLivetime_s`, `muonModeCode`, and `muonEnergySampledGeV` plus the existing positions (`muonX0/Y0/Z0`, `muonXImp/YImp/ZImp`) and angles.

Validation: use `analysis/validate_muon_tierB.py --input <B02ntuples.root>` to dump cos(theta), energy, and `(x0,y0)` histograms alongside the hit fraction and effective livetime.

## References
- Guan, M. et al., *"A parametrization of the cosmic-ray muon flux at sea-level"*, 2015.
- Particle Data Group, *Cosmic Rays* review, standard Gaisser parameterisation.
