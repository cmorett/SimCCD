# Muon generator validation (ROOT)

## 1) Produce the ntuple (Geant4)

Build:

- `cmake --preset vs2022`
- `cmake --build --preset vs2022-release`

Run the validation macro (fast primitive geometry, no vis):

- From `build/vs2022/main/Release`:
  - `.\b02_executable.exe --geometry primitive --no-vis ..\..\..\..\run\run_muongen_validation.mac`

This produces `build/vs2022/main/Release/B02ntuples.root` with `TTree` `B02Evts` and (at least) these branches:

- `thetaPri`, `EevtPri`
- `muonX0`, `muonY0`, `muonZ0` (primary vertex)
- `muonXImp`, `muonYImp`, `muonZImp` (impact-plane sampling point from generator)

Notes:

- The generator defaults to **fixed energy** (`/generator/useFixedEnergy true`), but the validation macro explicitly sets `/generator/useFixedEnergy false` to exercise the spectrum sampling.
- The validation macro sets fixed seeds via `/runAction/*` so the run is repeatable.

## 2) Generate PDFs (ROOT, batch)

From the repo root:

- `root -l -b -q ROOTanalysis/run_validate_muongen.C`

This writes (in the repo root):

- `validate_energy.pdf`
- `validate_costheta.pdf`
- `validate_xy0.pdf`
- `validate_z0.pdf`

Optional impact-plane validation:

- `set SIMCCD_VALIDATE_IMPACT=1` (cmd.exe) or `$env:SIMCCD_VALIDATE_IMPACT='1'` (PowerShell)
- Run the same ROOT command again: `root -l -b -q ROOTanalysis/run_validate_muongen.C`

This additionally writes:

- `validate_xyImpact.pdf`

## What each plot means

- `validate_costheta.pdf`: histogram of **downward zenith cosine** `cos_zenith^↓ = -cos(thetaPri)` clamped into `[0,1]`.
- `validate_energy.pdf`: primary muon kinetic energy `EevtPri` in **GeV**.
- `validate_xy0.pdf`: `muonX0/muonY0` = **primary vertex** position (not the impact-plane point).
- `validate_z0.pdf`: `muonZ0` = primary vertex height; if constant, the plot uses a tight axis range.
- `validate_xyImpact.pdf`: `muonXImp/muonYImp` = **impact-plane sampling point** (should be uniform over the configured plane).

