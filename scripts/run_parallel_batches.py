#!/usr/bin/env python3
"""
Run multiple single-threaded Geant4 jobs in parallel and merge ROOT outputs.
"""
from __future__ import annotations

import argparse
import math
import subprocess
import time
from pathlib import Path
from typing import List

import uproot


def split_events(total: int, jobs: int) -> List[int]:
    base = total // jobs
    remainder = total % jobs
    counts = [base + (1 if i < remainder else 0) for i in range(jobs)]
    return [c for c in counts if c > 0]


def write_macro(path: Path, events: int, output_root: Path, seed1: int, seed2: int, tag: str) -> None:
    lines = [
        "/control/execute main/macros/run_tierB_paper_default_1M.mac",
        f"/analysis/setFileName {output_root.as_posix()}",
        f"/analysis/ntuple/setFileName 0 {output_root.as_posix()}",
        f"/analysis/ntuple/setFileName 1 {output_root.as_posix()}",
        "/runAction/useTimeSeed false",
        "/runAction/useFixedSeeds true",
        f"/runAction/seed1 {seed1}",
        f"/runAction/seed2 {seed2}",
        f"/runAction/macroPath {path.as_posix()}",
        f"/runAction/provenanceTag {tag}",
        "/run/printProgress 1000",
        "/run/initialize",
        f"/run/beamOn {events}",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def merge_trees(output_root: Path, inputs: List[Path], tree_name: str) -> None:
    arrays = uproot.concatenate([f"{p}:{tree_name}" for p in inputs], library="np")
    with uproot.recreate(output_root) as fout:
        fout[tree_name] = arrays


def merge_outputs(output_root: Path, inputs: List[Path]) -> None:
    trees = ["B02Evts", "B02Hits", "B02RunInfo"]
    with uproot.recreate(output_root) as fout:
        for tree in trees:
            try:
                arrays = uproot.concatenate([f"{p}:{tree}" for p in inputs], library="np")
            except Exception:
                continue
            fout[tree] = arrays


def main() -> int:
    parser = argparse.ArgumentParser(description="Run parallel Geant4 jobs and merge outputs.")
    parser.add_argument("--jobs", type=int, default=4, help="Number of parallel jobs.")
    parser.add_argument("--events", type=int, default=5000, help="Total events to run.")
    parser.add_argument("--geometry", default="cad", choices=["cad", "primitive"])
    parser.add_argument("--cad-mode", default="tessellated")
    parser.add_argument("--cad-file", default="main/assets/cad/new_assembly_cleaned_assimp.dae")
    parser.add_argument("--output-tag", default="tierB_parallel")
    parser.add_argument("--exe", default="build/vs2022/main/Release/b02_executable.exe")
    args = parser.parse_args()

    jobs = max(1, args.jobs)
    counts = split_events(args.events, jobs)
    if not counts:
        raise SystemExit("No events to run.")

    out_dir = Path("build/vs2022/main/Release")
    macro_dir = Path("run/auto")
    log_dir = Path("outputs/parallel_logs")
    macro_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)

    base_seed1 = 12345
    base_seed2 = 67890
    chunk_paths = []
    procs = []

    for idx, count in enumerate(counts):
        tag = f"{args.output_tag}_chunk{idx}"
        out_root = out_dir / f"{tag}.root"
        macro_path = macro_dir / f"{tag}.mac"
        seed1 = base_seed1 + idx * 101
        seed2 = base_seed2 + idx * 103
        write_macro(macro_path, count, out_root, seed1, seed2, tag)
        log_path = log_dir / f"{tag}.log"
        cmd = [
            args.exe,
            "--geometry",
            args.geometry,
            "--cad-mode",
            args.cad_mode,
            "--cad-file",
            args.cad_file,
            "--no-vis",
            str(macro_path),
        ]
        procs.append(subprocess.Popen(cmd, stdout=log_path.open("w"), stderr=subprocess.STDOUT))
        chunk_paths.append(out_root)

    start = time.time()
    for proc in procs:
        proc.wait()
    elapsed = time.time() - start

    merged_root = out_dir / f"{args.output_tag}.root"
    merge_outputs(merged_root, chunk_paths)

    print(f"elapsed_seconds={elapsed:.3f}")
    print(f"merged_root={merged_root.as_posix()}")
    print(f"chunks={len(chunk_paths)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
