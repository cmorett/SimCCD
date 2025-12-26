# Paper v1 reproducibility (Tier-B, primitive geometry)

## Baseline thickness
- CCD thickness: 725 um (0.0725 cm).
- Defined in `main/src/B02DetectorConstruction.cc`:
  - `BuildPrimitiveGeometry`: `halfZ = 0.3625 * mm` (0.725 mm total).
  - `BuildCCDOverlay`: `halfZ = 0.3625 * mm` (0.725 mm total).
- Debug print: `/sim/ccd/printInfo true` (enabled in `main/macros/run_tierB_paper_default_1M.mac`).

## Geant4 run (50k throws)
Macro used:
- Base macro: `main/macros/run_tierB_paper_default_1M.mac`
- Driver macro: `run/paper_v1_driver_50k.mac`

Command line:
```
mkdir outputs\paper_v1_50k
main\build\Release\b02_executable.exe --geometry primitive --no-vis run\paper_v1_driver_50k.mac
```

## Analysis (paper outputs)
```
python analysis/make_paper_outputs.py --input outputs/paper_v1_50k/B02ntuples_paper_v1_50k.root --output paper_outputs --tag paper_v1_50k --examples 24 --dist-events 50000 --pixel-size-microns 15 --thickness-microns 725 --canvas-mode adaptive --margin-pix 24 --quality-only
```

## Outputs
- Figures/tables live under `paper_outputs/paper_v1_50k/`.

## Environment notes
- Geant4: 11.3.2 (from run log).
- Date: 2025-12-26.
- Commit: paper-v1 tag (resolve with `git rev-parse paper-v1`).
