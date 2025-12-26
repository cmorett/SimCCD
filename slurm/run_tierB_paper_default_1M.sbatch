#!/usr/bin/env bash
#SBATCH --job-name=simccd_1M
#SBATCH --output=%x_%j.log
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --ntasks=1

set -euo pipefail

module purge
module load lamod/GEANT4/11.1.1
module load lamod/cmake/3.31

REPO_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
BUILD_DIR="${BUILD_DIR:-$REPO_DIR/build/xook}"
OUT_BASE="${OUT_BASE:-$REPO_DIR/outputs}"
RUN_TAG="${RUN_TAG:-tierB_1M_${SLURM_JOB_ID:-local}_$(date +%Y%m%d_%H%M%S)}"
OUT_DIR="${OUT_BASE}/${RUN_TAG}"

EVENTS="${EVENTS:-1000000}"
GEOMETRY="${GEOMETRY:-primitive}"
SEED1="${SEED1:-12345}"
SEED2="${SEED2:-67890}"
MACRO_SRC="${MACRO_SRC:-$REPO_DIR/main/macros/run_tierB_paper_default_1M.mac}"

mkdir -p "$BUILD_DIR" "$OUT_DIR"

export SIMCCD_ASSETS_DIR="${SIMCCD_ASSETS_DIR:-$REPO_DIR/main/assets}"

if [[ ! -f "$MACRO_SRC" ]]; then
  echo "[simccd_1M] Macro not found: $MACRO_SRC"
  exit 1
fi

if [[ ! -x "$BUILD_DIR/b02_executable" ]]; then
  cmake -S "$REPO_DIR/main" -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE=Release
  cmake --build "$BUILD_DIR" --config Release -j "${SIMCCD_BUILD_JOBS:-8}"
fi

EXE="$BUILD_DIR/b02_executable"
if [[ ! -x "$EXE" ]] && [[ -x "$BUILD_DIR/Release/b02_executable" ]]; then
  EXE="$BUILD_DIR/Release/b02_executable"
elif [[ ! -x "$EXE" ]] && [[ -x "$BUILD_DIR/main/b02_executable" ]]; then
  EXE="$BUILD_DIR/main/b02_executable"
fi
if [[ ! -x "$EXE" ]]; then
  echo "[simccd_1M] Could not find built executable under $BUILD_DIR"
  exit 1
fi

MACRO_COPY="$OUT_DIR/run_tierB_paper_default_1M.mac"
cp "$MACRO_SRC" "$MACRO_COPY"

RUN_MACRO="$OUT_DIR/run_driver_1M.mac"
cat > "$RUN_MACRO" <<EOF
/control/execute $MACRO_COPY
/analysis/setFileName B02
/analysis/ntuple/setFileName 0 B02ntuples
/analysis/ntuple/setFileName 1 B02ntuples
/analysis/ntuple/setFileName 2 B02ntuples
/runAction/useTimeSeed false
/runAction/useFixedSeeds true
/runAction/seed1 $SEED1
/runAction/seed2 $SEED2
/runAction/macroPath $RUN_MACRO
/runAction/provenanceTag $RUN_TAG
/run/printProgress 10000
/run/initialize
/run/beamOn $EVENTS
EOF

cd "$OUT_DIR"
echo "[simccd_1M] Running $EXE --geometry $GEOMETRY --no-vis $RUN_MACRO"
"$EXE" --geometry "$GEOMETRY" --no-vis "$RUN_MACRO" > "$OUT_DIR/run_${RUN_TAG}.log" 2>&1

if [[ -f "$OUT_DIR/B02ntuples.root" ]]; then
  mv "$OUT_DIR/B02ntuples.root" "$OUT_DIR/B02ntuples_${RUN_TAG}.root"
fi
if [[ -f "$OUT_DIR/B02.root" ]]; then
  mv "$OUT_DIR/B02.root" "$OUT_DIR/B02_${RUN_TAG}.root"
fi

echo "[simccd_1M] DONE. Outputs in $OUT_DIR"
