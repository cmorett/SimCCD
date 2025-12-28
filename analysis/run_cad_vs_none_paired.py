#!/usr/bin/env python3
"""
Run paired CAD-vs-none batches deterministically, in parallel, and merge outputs.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Sequence

try:
    import uproot
except ImportError as exc:
    raise SystemExit("This script requires the 'uproot' package. Install via pip.") from exc


SKIP_PREFIXES = (
    "/analysis/setFileName",
    "/analysis/ntuple/setFileName",
    "/runAction/useTimeSeed",
    "/runAction/useFixedSeeds",
    "/runAction/seed1",
    "/runAction/seed2",
    "/runAction/macroPath",
    "/runAction/provenanceTag",
    "/run/printProgress",
    "/run/initialize",
    "/run/beamOn",
)


@dataclass
class ModeConfig:
    name: str
    cad_mode: str
    thrown: int


@dataclass
class BatchJob:
    mode: str
    index: int
    events: int
    seed1: int
    seed2: int
    macro_path: Path
    output_root: Path
    log_path: Path
    provenance_tag: str


@dataclass
class ModeResult:
    mode: str
    batches: int
    events: int
    elapsed_seconds: float
    per_batch_seconds: float
    merged_root: Path
    batch_files: List[Path]


@dataclass
class RunResult:
    tag: str
    out_root: Path
    none: ModeResult
    cad: ModeResult
    metadata_path: Path


def split_batches(total: int, batch_size: int) -> List[int]:
    if batch_size <= 0:
        raise ValueError("batch_size must be positive.")
    counts = []
    remaining = total
    while remaining > 0:
        counts.append(min(batch_size, remaining))
        remaining -= batch_size
    return counts


def read_macro_lines(path: Path) -> List[str]:
    return path.read_text(encoding="utf-8").splitlines()


def filter_macro_lines(lines: Sequence[str]) -> List[str]:
    filtered = []
    for line in lines:
        stripped = line.strip()
        if not stripped:
            filtered.append(line)
            continue
        if any(stripped.startswith(prefix) for prefix in SKIP_PREFIXES):
            continue
        filtered.append(line)
    return filtered


def write_macro(
    path: Path,
    base_lines: Sequence[str],
    output_root: Path,
    seed1: int,
    seed2: int,
    provenance_tag: str,
    events: int,
    print_progress: int,
) -> None:
    lines = []
    lines.append(f"# Auto-generated macro: {path.name}")
    lines.extend(base_lines)
    out_path = output_root.as_posix()
    lines.extend(
        [
            f"/analysis/setFileName {out_path}",
            f"/analysis/ntuple/setFileName 0 {out_path}",
            f"/analysis/ntuple/setFileName 1 {out_path}",
            "/runAction/useTimeSeed false",
            "/runAction/useFixedSeeds true",
            f"/runAction/seed1 {seed1}",
            f"/runAction/seed2 {seed2}",
            f"/runAction/macroPath {path.as_posix()}",
            f"/runAction/provenanceTag {provenance_tag}",
            f"/run/printProgress {print_progress}",
            "/run/initialize",
            f"/run/beamOn {events}",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def is_valid_root(path: Path) -> bool:
    if not path.exists() or path.stat().st_size == 0:
        return False
    try:
        with uproot.open(path) as root:
            if "B02Evts" not in root:
                return False
            if len(root["B02Evts"]) == 0:
                return False
    except Exception:
        return False
    return True


def run_batch(
    exe: Path,
    job: BatchJob,
    geometry_mode: str,
    cad_file: Optional[Path],
    assets_dir: Optional[Path],
    dry_run: bool,
) -> float:
    cmd = [str(exe), "--geometry", geometry_mode, "--cad-mode", job.mode]
    if cad_file:
        cmd += ["--cad-file", str(cad_file)]
    if assets_dir:
        cmd += ["--assets-dir", str(assets_dir)]
    cmd += ["--no-vis", str(job.macro_path)]
    if dry_run:
        print("[dry-run] " + " ".join(cmd))
        return 0.0
    job.log_path.parent.mkdir(parents=True, exist_ok=True)
    start = time.time()
    with job.log_path.open("w", encoding="utf-8") as log:
        subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=True)
    return time.time() - start


def merge_outputs(output_root: Path, inputs: Sequence[Path]) -> None:
    trees = ["B02Evts", "B02Hits", "B02RunInfo"]
    with uproot.recreate(output_root) as fout:
        for tree in trees:
            try:
                arrays = uproot.concatenate([f"{p}:{tree}" for p in inputs], library="np")
            except Exception:
                continue
            fout[tree] = arrays


def build_jobs(
    mode: ModeConfig,
    counts: Sequence[int],
    seed_base: int,
    out_root: Path,
    base_lines: Sequence[str],
    print_progress: int,
    tag: str,
    write_files: bool,
) -> List[BatchJob]:
    jobs = []
    macros_dir = out_root / mode.name / "macros"
    logs_dir = out_root / mode.name / "logs"
    batches_dir = out_root / mode.name / "batches"
    if write_files:
        macros_dir.mkdir(parents=True, exist_ok=True)
        logs_dir.mkdir(parents=True, exist_ok=True)
        batches_dir.mkdir(parents=True, exist_ok=True)

    for idx, events in enumerate(counts):
        seed1 = seed_base + idx * 2
        seed2 = seed_base + idx * 2 + 1
        batch_tag = f"{tag}_{mode.name}_batch{idx:05d}"
        macro_path = macros_dir / f"batch_{idx:05d}.mac"
        output_root = batches_dir / f"batch_{idx:05d}.root"
        log_path = logs_dir / f"batch_{idx:05d}.log"
        if write_files:
            write_macro(
                macro_path,
                base_lines,
                output_root,
                seed1,
                seed2,
                batch_tag,
                events,
                print_progress,
            )
        jobs.append(
            BatchJob(
                mode=mode.cad_mode,
                index=idx,
                events=events,
                seed1=seed1,
                seed2=seed2,
                macro_path=macro_path,
                output_root=output_root,
                log_path=log_path,
                provenance_tag=batch_tag,
            )
        )
    return jobs


def run_mode(
    exe: Path,
    mode: ModeConfig,
    counts: Sequence[int],
    seed_base: int,
    out_root: Path,
    base_lines: Sequence[str],
    tag: str,
    geometry_mode: str,
    cad_file: Optional[Path],
    assets_dir: Optional[Path],
    max_workers: int,
    resume: bool,
    dry_run: bool,
) -> ModeResult:
    print_progress = 10000 if max(counts) >= 10000 else 1000
    jobs = build_jobs(
        mode,
        counts,
        seed_base,
        out_root,
        base_lines,
        print_progress,
        tag,
        write_files=not dry_run,
    )
    to_run: List[BatchJob] = []
    batch_files = []
    for job in jobs:
        batch_files.append(job.output_root)
        if resume and is_valid_root(job.output_root):
            continue
        to_run.append(job)

    if dry_run:
        for job in to_run:
            run_batch(exe, job, geometry_mode, cad_file, assets_dir, dry_run=True)
        merged_root = out_root / mode.name / "merged.root"
        return ModeResult(
            mode=mode.name,
            batches=len(counts),
            events=mode.thrown,
            elapsed_seconds=0.0,
            per_batch_seconds=0.0,
            merged_root=merged_root,
            batch_files=batch_files,
        )

    start = time.time()
    if to_run:
        from concurrent.futures import ThreadPoolExecutor, as_completed

        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            future_map = {
                executor.submit(run_batch, exe, job, geometry_mode, cad_file, assets_dir, False): job
                for job in to_run
            }
            for future in as_completed(future_map):
                future.result()
    elapsed = time.time() - start

    missing = [p for p in batch_files if not is_valid_root(p)]
    if missing:
        missing_list = "\n".join(str(p) for p in missing)
        raise SystemExit(f"Missing or invalid batch outputs for {mode.name}:\n{missing_list}")

    merged_root = out_root / mode.name / "merged.root"
    merge_outputs(merged_root, batch_files)

    per_batch = elapsed / max(1, len(to_run))
    return ModeResult(
        mode=mode.name,
        batches=len(counts),
        events=mode.thrown,
        elapsed_seconds=elapsed,
        per_batch_seconds=per_batch,
        merged_root=merged_root,
        batch_files=batch_files,
    )


def run_paired(
    tag: str,
    exe: Path,
    macro: Path,
    geometry_none: str,
    geometry_cad: str,
    thrown_none: int,
    thrown_cad: int,
    batch_size: int,
    max_workers: int,
    seed_base: int,
    out_root: Path,
    cad_file: Optional[Path],
    assets_dir: Optional[Path],
    geometry_mode: str,
    resume: bool,
    dry_run: bool,
) -> RunResult:
    if not exe.exists():
        raise SystemExit(f"Executable not found: {exe}")
    if not macro.exists():
        raise SystemExit(f"Macro not found: {macro}")

    base_lines = filter_macro_lines(read_macro_lines(macro))
    counts_none = split_batches(thrown_none, batch_size)
    counts_cad = split_batches(thrown_cad, batch_size)

    if dry_run:
        print("[dry-run] Paired run plan")
        print(f"  tag={tag}")
        print(f"  exe={exe}")
        print(f"  macro={macro}")
        print(f"  none={thrown_none} events, {len(counts_none)} batches")
        print(f"  cad={thrown_cad} events, {len(counts_cad)} batches")

    max_workers = max(1, max_workers)
    none_cfg = ModeConfig(name="none", cad_mode=geometry_none, thrown=thrown_none)
    cad_cfg = ModeConfig(name="cad", cad_mode=geometry_cad, thrown=thrown_cad)

    if not dry_run:
        out_root.mkdir(parents=True, exist_ok=True)
    none_result = run_mode(
        exe,
        none_cfg,
        counts_none,
        seed_base,
        out_root,
        base_lines,
        tag,
        geometry_mode,
        cad_file,
        assets_dir,
        max_workers,
        resume,
        dry_run,
    )
    cad_result = run_mode(
        exe,
        cad_cfg,
        counts_cad,
        seed_base,
        out_root,
        base_lines,
        tag,
        geometry_mode,
        cad_file,
        assets_dir,
        max_workers,
        resume,
        dry_run,
    )

    metadata = {
        "tag": tag,
        "exe": str(exe),
        "macro": str(macro),
        "geometry_mode": geometry_mode,
        "geometry_none": geometry_none,
        "geometry_cad": geometry_cad,
        "cad_file": str(cad_file) if cad_file else "",
        "assets_dir": str(assets_dir) if assets_dir else "",
        "thrown_none": thrown_none,
        "thrown_cad": thrown_cad,
        "batch_size": batch_size,
        "max_workers": max_workers,
        "seed_base": seed_base,
        "none": {
            "batches": none_result.batches,
            "events": none_result.events,
            "elapsed_seconds": none_result.elapsed_seconds,
            "per_batch_seconds": none_result.per_batch_seconds,
            "merged_root": str(none_result.merged_root),
        },
        "cad": {
            "batches": cad_result.batches,
            "events": cad_result.events,
            "elapsed_seconds": cad_result.elapsed_seconds,
            "per_batch_seconds": cad_result.per_batch_seconds,
            "merged_root": str(cad_result.merged_root),
        },
    }
    metadata_path = out_root / "run_metadata.json"
    if not dry_run:
        metadata_path.write_text(json.dumps(metadata, indent=2), encoding="utf-8")

    return RunResult(
        tag=tag,
        out_root=out_root,
        none=none_result,
        cad=cad_result,
        metadata_path=metadata_path,
    )


def main() -> int:
    parser = argparse.ArgumentParser(description="Run paired CAD-vs-none batches.")
    parser.add_argument("--tag", required=True, help="Run tag.")
    parser.add_argument("--exe", required=True, help="Path to b02_executable.")
    parser.add_argument("--mac", required=True, help="Macro file with generator settings.")
    parser.add_argument("--geometry-none", default="none", help="CAD mode for no-geometry run.")
    parser.add_argument("--geometry-cad", default="tessellated", help="CAD mode for CAD run.")
    parser.add_argument("--thrown-none", type=int, required=True, help="Total events for none.")
    parser.add_argument("--thrown-cad", type=int, required=True, help="Total events for CAD.")
    parser.add_argument("--batch-size", type=int, default=5000, help="Events per batch.")
    parser.add_argument("--max-workers", type=int, default=4, help="Parallel workers.")
    parser.add_argument("--seed-base", type=int, default=12345, help="Seed base for batches.")
    parser.add_argument("--out-root", default="", help="Output root directory.")
    parser.add_argument("--cad-file", default="", help="CAD file path.")
    parser.add_argument("--assets-dir", default="", help="Assets directory.")
    parser.add_argument("--geometry-mode", default="cad", help="Geometry mode (cad or primitive).")
    parser.add_argument("--resume", action="store_true", help="Skip batches with valid outputs.")
    parser.add_argument("--dry-run", action="store_true", help="Print actions without running.")
    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parent.parent
    exe = (repo_root / args.exe).resolve() if not Path(args.exe).is_absolute() else Path(args.exe)
    macro = (repo_root / args.mac).resolve() if not Path(args.mac).is_absolute() else Path(args.mac)
    cad_file = None
    if args.cad_file:
        cad_file = (repo_root / args.cad_file).resolve() if not Path(args.cad_file).is_absolute() else Path(args.cad_file)
    assets_dir = None
    if args.assets_dir:
        assets_dir = (
            (repo_root / args.assets_dir).resolve()
            if not Path(args.assets_dir).is_absolute()
            else Path(args.assets_dir)
        )

    out_root = Path(args.out_root) if args.out_root else repo_root / "outputs" / args.tag
    result = run_paired(
        tag=args.tag,
        exe=exe,
        macro=macro,
        geometry_none=args.geometry_none,
        geometry_cad=args.geometry_cad,
        thrown_none=args.thrown_none,
        thrown_cad=args.thrown_cad,
        batch_size=args.batch_size,
        max_workers=args.max_workers,
        seed_base=args.seed_base,
        out_root=out_root,
        cad_file=cad_file,
        assets_dir=assets_dir,
        geometry_mode=args.geometry_mode,
        resume=args.resume,
        dry_run=args.dry_run,
    )

    if args.dry_run:
        return 0
    print(f"[done] none merged: {result.none.merged_root}")
    print(f"[done] cad merged: {result.cad.merged_root}")
    print(f"[done] metadata: {result.metadata_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
