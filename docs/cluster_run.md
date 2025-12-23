# Cluster / batch runs (10k quick test, 1M production)

## Quick 10k smoke test

```
# build once
cmake --preset vs2022
cmake --build --preset vs2022-release

# run 10k events (deterministic seeds, spectrum mode)
./build/vs2022/main/Release/b02_executable.exe --geometry primitive --no-vis run/run_simccd_10k.mac
```

Outputs: `B02ntuples_10k.root` (plus log on stdout).

## One-file, cluster-friendly 1M run (Linux)

```
bash scripts/run_simccd_1M.sh --exe ./build/vs2022/main/Release/b02_executable.exe \
  --geometry primitive \
  --out outputs --tag runA \
  --events 1000000 --chunks 10
```

What it does:

- Generates per-chunk macros with deterministic seeds (`seed1/seed2` increment per chunk).
- Runs in batch (`--no-vis`), prints progress every 10k.
- Outputs per-chunk ROOT files + logs under `outputs/runA/`.
- Writes `run_metadata.json` summarizing config.
- If `hadd` is available, merges into `B02ntuples_runA_merged.root`.

Flags:

- `--geometry cad` to run CAD geometry.
- `--events N` to change total events (default 1,000,000).
- `--chunks N` to split into N jobs; events per chunk is `ceil(events/chunks)`.
- `--seed1/--seed2` to override deterministic seeds.
- `--make-paper-outputs` will run `analysis/make_paper_outputs.py` on the merged ROOT file (if produced).

## Windows batch (PowerShell)

```
.\scripts\run_simccd_1M.ps1 -Exe .\build\vs2022\main\Release\b02_executable.exe `
  -Geometry primitive -Out outputs -Tag runA -Events 1000000 -Chunks 10 -MakePaperOutputs
```

## Notes

- Default macros for big runs: `run/run_simccd_1M.mac` (1M) and `run/run_simccd_10k.mac` (10k).
- Both macros force the muon spectrum (`/generator/useFixedEnergy false`) and deterministic seeds.
- Logs land next to the outputs; failures return non-zero.
