#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_DIR="${BUILD_DIR:-$ROOT_DIR/build/xook}"
RUN_DIR="${RUN_DIR:-$ROOT_DIR/outputs/cluster_smoke_$(date +%Y%m%d_%H%M%S)}"
GEOMETRY="${GEOMETRY:-primitive}"
EVENTS="${EVENTS:-100}"
SEED1="${SEED1:-12345}"
SEED2="${SEED2:-67890}"
MACRO_BASE="${MACRO_BASE:-$ROOT_DIR/main/macros/run_tierB_paper_default_1M.mac}"

module purge
module load lamod/GEANT4/11.1.1
module load lamod/cmake/3.31

if [[ ! -f "$MACRO_BASE" ]]; then
  echo "[cluster_smoke_test] Macro not found: $MACRO_BASE"
  exit 1
fi

mkdir -p "$BUILD_DIR" "$RUN_DIR"

export SIMCCD_ASSETS_DIR="${SIMCCD_ASSETS_DIR:-$ROOT_DIR/main/assets}"

cmake -S "$ROOT_DIR/main" -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE=Release
cmake --build "$BUILD_DIR" --config Release -j "${SIMCCD_BUILD_JOBS:-8}"

EXE="$BUILD_DIR/b02_executable"
if [[ ! -x "$EXE" ]] && [[ -x "$BUILD_DIR/Release/b02_executable" ]]; then
  EXE="$BUILD_DIR/Release/b02_executable"
elif [[ ! -x "$EXE" ]] && [[ -x "$BUILD_DIR/main/b02_executable" ]]; then
  EXE="$BUILD_DIR/main/b02_executable"
fi
if [[ ! -x "$EXE" ]]; then
  echo "[cluster_smoke_test] Could not find built executable under $BUILD_DIR"
  exit 1
fi

MACRO_BASE_ABS="$(cd "$(dirname "$MACRO_BASE")" && pwd)/$(basename "$MACRO_BASE")"
RUN_MACRO="$RUN_DIR/run_smoke.mac"

cat > "$RUN_MACRO" <<EOF
/control/execute $MACRO_BASE_ABS
/analysis/setFileName B02
/analysis/ntuple/setFileName 0 B02ntuples
/analysis/ntuple/setFileName 1 B02ntuples
/analysis/ntuple/setFileName 2 B02ntuples
/runAction/useTimeSeed false
/runAction/useFixedSeeds true
/runAction/seed1 $SEED1
/runAction/seed2 $SEED2
/runAction/macroPath $RUN_MACRO
/runAction/provenanceTag cluster_smoke
/run/printProgress 50
/run/initialize
/run/beamOn $EVENTS
EOF

cd "$RUN_DIR"
echo "[cluster_smoke_test] Running $EXE --geometry $GEOMETRY --no-vis $RUN_MACRO"
"$EXE" --geometry "$GEOMETRY" --no-vis "$RUN_MACRO" > "$RUN_DIR/run_smoke.log" 2>&1

if [[ ! -f "$RUN_DIR/B02ntuples.root" ]]; then
  echo "[cluster_smoke_test] Missing output: $RUN_DIR/B02ntuples.root"
  exit 1
fi
if [[ ! -f "$RUN_DIR/B02.root" ]]; then
  echo "[cluster_smoke_test] Missing output: $RUN_DIR/B02.root"
  exit 1
fi

echo "[cluster_smoke_test] OK. Outputs in $RUN_DIR"
