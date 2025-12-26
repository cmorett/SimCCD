# SimCCD

Geant4-based SimCCD with CAD import. The build and runtime no longer depend on the working directory; assets and CAD paths are resolved via CLI, env vars, or the CMake-configured defaults.

## Build
- Requires Geant4 11.x, CMake, Assimp (via vcpkg or system), compiler toolchain.
- Configure from the project root (Windows example):
  ```
  cmake -S main -B main/build -DCMAKE_BUILD_TYPE=Release
  cmake --build main/build --config Release
  ```
  CMake generates `main/build/Release/b02_executable.exe` (or `b02_executable` on Linux). Assets are copied into the build tree automatically.

## Runtime options
- `--geometry cad|primitive` (default: cad)
- `--cad-mode merged|tessellated|parts` (default: merged)
  - `merged`: stable bounding-box multi-union (no overlap warnings)
  - `tessellated`: merged raw meshes (may be slower)
  - `parts`: individual meshes placed separately (may produce overlaps)
- `--cad-file <path>` or `--assets-dir <dir>`
- Env vars: `SIMCCD_CAD_FILE`, `SIMCCD_ASSETS_DIR` (or use defaults compiled via `SimCCDConfig.hh`).
- `--no-vis` to skip visualization in batch mode.

Example:
```
main/build/Release/b02_executable.exe --no-vis --cad-mode merged main/macros/run_with_cad.mac
main/build/Release/b02_executable.exe --no-vis --geometry primitive main/macros/run_without_cad.mac
```

## Macros
- `main/macros/geom_check.mac` – `/run/initialize` + `/geometry/test/run` for quick geometry validation.
- `main/macros/run_with_cad.mac` – production macro for CAD.
- `main/macros/run_without_cad.mac` – production macro for primitive geometry.
- `main/macros/run_tierB_cosmic_muons.mac` – Tier B cosmic muons with Modified Gaisser sampling on an extended source plane (nontrivial hit/miss fraction).
- `main/macros/run_tierB_paper_default_1M.mac` – paper defaults for Tier B (10×10 cm plane at z=8 cm, Emin=1 GeV, muon charge ratio 1.25, overburden/cuts off).

## Cosmic muon generator modes
- Legacy (default): `forced_footprint` keeps the back-projected, mostly-hits CCD generator.
- Tier B: `tierB_plane_flux` samples direction+energy from Guan2015 Modified Gaisser with a configurable source plane and optional auto-sizing. Use the macros above or set `/generator/muonMode tierB_plane_flux` plus plane/energy/charge commands. Stores per-event weights/livetime for rate conversion. See `docs/muon_sampling.md` for details.

## Scripts
- `scripts/validate_geometry.py` – fails on overlap/navigation warnings.
- `scripts/benchmark_geometry.py` – runs CAD vs primitive and writes `benchmark_geometry.pdf` (matplotlib required).
- `scripts/run_simccd_1M.sh` – production runner that now drives Tier B by default using `run_tierB_paper_default_1M.mac`; accepts `--events/--chunks`, `--geometry cad|primitive`, and `--make-paper-outputs`.

## Cluster
- See `docs/cluster.md` for ICN-UNAM notes and the SLURM helper `slurm/run_geometry.sbatch`.
