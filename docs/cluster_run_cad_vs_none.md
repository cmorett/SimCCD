# Xook CAD vs none runs

This guide runs the paired pipeline on the UNAM xook cluster (SLURM).

## Environment setup (login shell)

```bash
module purge
module load lamod/GEANT4/11.1.1 lamod/cmake/3.31 lamod/Anaconda/202309
source /software/lamod/bin/Anaconda/202309/etc/profile.d/conda.sh
conda activate /storage/icn/carlos.morett/opt/conda_envs/simccd-py311
cd /storage/icn/carlos.morett/repos/SimCCD
```

## Build (only if needed)

```bash
cmake -S . -B /storage/icn/carlos.morett/build/simccd
cmake --build /storage/icn/carlos.morett/build/simccd -j 8
```

## Smoke test (fast)

```bash
export SIMCCD_EXE=/storage/icn/carlos.morett/build/simccd/b02_executable
export SIMCCD_PYTHON="$(which python3)"
python3 analysis/cluster_run_cad_vs_none_prod.py \
  --config configs/cad_vs_none_prod.yaml \
  --skip-sanity \
  --thrown-none 2000 \
  --thrown-cad 2000 \
  --batch-size 1000 \
  --max-workers 2 \
  --tag cad_vs_none_smoke
```

## Production run (paired none vs tessellated CAD)

```bash
export SIMCCD_EXE=/storage/icn/carlos.morett/build/simccd/b02_executable
export SIMCCD_PYTHON="$(which python3)"
python3 analysis/cluster_run_cad_vs_none_prod.py \
  --config configs/cad_vs_none_prod.yaml \
  --resume
```

## SLURM submission (one-line)

```bash
sbatch --job-name=cad_vs_none_prod \
  --output=/storage/icn/carlos.morett/logs/%x_%j.out \
  --error=/storage/icn/carlos.morett/logs/%x_%j.err \
  --wrap="bash -lc 'module purge && \
    module load lamod/GEANT4/11.1.1 lamod/cmake/3.31 lamod/Anaconda/202309 && \
    source /software/lamod/bin/Anaconda/202309/etc/profile.d/conda.sh && \
    conda activate /storage/icn/carlos.morett/opt/conda_envs/simccd-py311 && \
    cd /storage/icn/carlos.morett/repos/SimCCD && \
    export SIMCCD_EXE=/storage/icn/carlos.morett/build/simccd/b02_executable && \
    export SIMCCD_PYTHON=\$(which python3) && \
    python3 analysis/cluster_run_cad_vs_none_prod.py --config configs/cad_vs_none_prod.yaml --resume'"
```

## Outputs and logs

- SLURM logs: `/storage/icn/carlos.morett/logs/%x_%j.out` and `.err`
- Paired outputs: `outputs/<tag>/none/merged.root` and `outputs/<tag>/cad/merged.root`
- Per-batch logs: `outputs/<tag>/<mode>/logs/*.log`
- Paper outputs: `paper_outputs/<tag>/none/`, `paper_outputs/<tag>/cad/`, `paper_outputs/<tag>/compare/`
- Summaries: `docs/cad_vs_none_summary_<tag>_sanity.md` and `docs/cad_vs_none_summary_<tag>.md`

Notes:
- `SIMCCD_EXE` overrides the `exe:` entry in `configs/cad_vs_none_prod.yaml`.
- `SIMCCD_PYTHON` controls which Python interpreter is used for all pipeline subprocesses.
