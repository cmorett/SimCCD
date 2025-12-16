#!/bin/bash
#SBATCH --job-name=simccd_geom
#SBATCH --output=%x_%j.log
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --ntasks=1

set -euo pipefail

# Adjust module names to the cluster environment.
module load gcc/12 cmake/3.26 geant4/11.1.1 root/6.30

ROOT_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
BUILD_DIR="$ROOT_DIR/build"
RUN_DIR="$ROOT_DIR/run_${SLURM_JOB_ID:-local}"
mkdir -p "$BUILD_DIR" "$RUN_DIR"

export SIMCCD_ASSETS_DIR="${SIMCCD_ASSETS_DIR:-$ROOT_DIR/main/assets}"

cmake -S "$ROOT_DIR/main" -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE=Release
cmake --build "$BUILD_DIR" --config Release

BINARY="$BUILD_DIR/b02_executable"
if [ ! -x "$BINARY" ] && [ -x "$BUILD_DIR/Release/b02_executable" ]; then
  BINARY="$BUILD_DIR/Release/b02_executable"
fi

cd "$RUN_DIR"
echo "Running in $RUN_DIR with assets at $SIMCCD_ASSETS_DIR"

"$BINARY" --no-vis --geometry cad --cad-mode merged "$ROOT_DIR/main/macros/run_with_cad.mac"

