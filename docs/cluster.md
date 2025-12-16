# Running on the ICN-UNAM Cluster

The defaults now work from any working directory as long as the CAD assets can be found. Use the environment variables below to pin the asset location when running on the cluster.

## Quick steps

1. Connect from PowerShell: `ssh carlos.morett@tochtli.icn.unam.mx`
2. Load tools (example module set – adjust to site names):
   ```bash
   module load gcc/12 cmake/3.26 geant4/11.1.1 root/6.30
   ```
3. Set up a build directory on the shared storage (avoid `$HOME` quotas):
   ```bash
   cd /storage/icn/carlos.morett/simccd
   export SIMCCD_ASSETS_DIR=$PWD/main/assets
   cmake -S main -B build -DCMAKE_BUILD_TYPE=Release
   cmake --build build --config Release
   ```
4. Run a quick validation (no visualization):
   ```bash
   ./build/b02_executable --no-vis --cad-mode merged main/macros/geom_check.mac
   ```
5. Submit batch jobs with the provided SLURM script:
   ```bash
   sbatch slurm/run_geometry.sbatch
   ```

## CAD path controls

- `--cad-file /path/to/assembly.dae`
- `--assets-dir /path/to/assets`
- `SIMCCD_CAD_FILE` or `SIMCCD_ASSETS_DIR`

If none are provided, the code falls back to the configured source/build asset paths from CMake, so it still works when launched from inside the build tree.

## Useful scripts

- `scripts/validate_geometry.py` – fails if overlap/navigation warnings appear.
- `scripts/benchmark_geometry.py` – compares CAD vs primitive runtime and writes `benchmark_geometry.pdf` (requires matplotlib).

