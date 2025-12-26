#!/usr/bin/env bash
set -euo pipefail

EXE="./b02_executable.exe"
GEOMETRY="primitive"
OUT_BASE="outputs"
TAG=""
NEVENTS=1000000
CHUNKS=1
SEED1=12345
SEED2=67890
MACRO_BASE="main/macros/run_tierB_paper_default_1M.mac"
MAKE_PAPER=0

usage() {
  cat <<EOF
Usage: $0 [--exe PATH] [--geometry primitive|cad] [--out DIR] [--tag NAME] [--events N] [--chunks N]
          [--seed1 N] [--seed2 N] [--macro PATH] [--make-paper-outputs]

Defaults: exe=$EXE, geometry=$GEOMETRY, out=$OUT_BASE, events=$NEVENTS, chunks=$CHUNKS
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --exe) EXE="$2"; shift 2;;
    --geometry) GEOMETRY="$2"; shift 2;;
    --out) OUT_BASE="$2"; shift 2;;
    --tag) TAG="$2"; shift 2;;
    --events) NEVENTS="$2"; shift 2;;
    --chunks) CHUNKS="$2"; shift 2;;
    --seed1) SEED1="$2"; shift 2;;
    --seed2) SEED2="$2"; shift 2;;
    --macro) MACRO_BASE="$2"; shift 2;;
    --make-paper-outputs) MAKE_PAPER=1; shift;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown arg: $1"; usage; exit 1;;
  esac
done

if [[ -z "$TAG" ]]; then
  TAG="run_$(date +%Y%m%d_%H%M%S)"
fi

OUT_DIR="${OUT_BASE}/${TAG}"
mkdir -p "$OUT_DIR"

# events per chunk (ceil)
events_per_chunk=$(( (NEVENTS + CHUNKS - 1) / CHUNKS ))

declare -a chunk_files

for ((i=0; i<CHUNKS; ++i)); do
  chunk_tag=$(printf "%s_chunk%02d" "$TAG" "$i")
  outfile="$OUT_DIR/B02ntuples_${chunk_tag}.root"
  logfile="$OUT_DIR/run_${chunk_tag}.log"
  macro="$OUT_DIR/macro_${chunk_tag}.mac"
  seed1=$((SEED1 + i))
  seed2=$((SEED2 + i))
  macro_abs=$(python - <<'PY'
import os,sys
path=sys.argv[1]
print(os.path.abspath(path))
PY
"$macro")
  cat > "$macro" <<EOF
/control/execute $MACRO_BASE
/analysis/setFileName $outfile
/analysis/ntuple/setFileName 0 $outfile
/analysis/ntuple/setFileName 1 $outfile
/runAction/useTimeSeed false
/runAction/useFixedSeeds true
/runAction/seed1 $seed1
/runAction/seed2 $seed2
/runAction/macroPath $macro_abs
/runAction/provenanceTag $chunk_tag
/run/printProgress 10000
/generator/muonMode tierB_plane_flux
/generator/useFixedEnergy false
/run/initialize
/run/beamOn $events_per_chunk
EOF
  echo "[run_simccd_1M] chunk $i: launching $EXE --geometry $GEOMETRY --no-vis $macro"
  set +e
  "$EXE" --geometry "$GEOMETRY" --no-vis "$macro" > "$logfile" 2>&1
  status=$?
  set -e
  if [[ $status -ne 0 ]]; then
    echo "[run_simccd_1M] chunk $i failed (status $status), see $logfile"
    exit $status
  fi
  if [[ ! -f "$outfile" ]]; then
    echo "[run_simccd_1M] chunk $i produced no ROOT file ($outfile)"
    exit 1
  fi
  chunk_files+=("$outfile")
done

merged=""
if [[ ${#chunk_files[@]} -gt 1 ]]; then
  if command -v hadd >/dev/null 2>&1; then
    merged="$OUT_DIR/B02ntuples_${TAG}_merged.root"
    echo "[run_simccd_1M] Merging chunks into $merged"
    hadd -f "$merged" "${chunk_files[@]}"
  else
    echo "[run_simccd_1M] hadd not found; skipping merge."
  fi
else
  merged="${chunk_files[0]}"
fi

cat > "$OUT_DIR/run_metadata.json" <<EOF
{
  "exe": "$EXE",
  "geometry": "$GEOMETRY",
  "events": $NEVENTS,
  "chunks": $CHUNKS,
  "events_per_chunk": $events_per_chunk,
  "seed1": $SEED1,
  "seed2": $SEED2,
  "macro_base": "$MACRO_BASE",
  "output_dir": "$OUT_DIR",
  "merged_file": "$merged"
}
EOF

if [[ $MAKE_PAPER -eq 1 && -n "$merged" ]]; then
  python analysis/make_paper_outputs.py --input "$merged" --output paper_outputs --tag "$TAG"
fi

echo "[run_simccd_1M] DONE. Outputs in $OUT_DIR"
