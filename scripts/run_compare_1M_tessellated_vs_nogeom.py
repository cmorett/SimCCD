#!/usr/bin/env python3
"""
Run 1M events for no-geometry and tessellated CAD, generate paper outputs,
and create side-by-side PDF comparisons.
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
import time
from pathlib import Path
from typing import Dict

try:
    from pypdf import PdfReader, PdfWriter, PageObject, Transformation
except ImportError as exc:
    raise SystemExit(
        "This script requires pypdf. Install it with: pip install pypdf"
    ) from exc


def default_jobs() -> int:
    cpu = os.cpu_count() or 4
    return max(1, min(12, cpu))


def run_cmd(cmd, cwd: Path) -> float:
    print("[run] " + " ".join(str(c) for c in cmd))
    start = time.time()
    subprocess.run(cmd, cwd=cwd, check=True)
    elapsed = time.time() - start
    print(f"[run] done in {elapsed:.1f} s")
    return elapsed


def run_parallel(
    repo_root: Path,
    exe: str,
    events: int,
    jobs: int,
    cad_mode: str,
    cad_file: str,
    output_tag: str,
) -> float:
    cmd = [
        sys.executable,
        "scripts/run_parallel_batches.py",
        "--events",
        str(events),
        "--jobs",
        str(jobs),
        "--geometry",
        "cad",
        "--cad-mode",
        cad_mode,
        "--cad-file",
        cad_file,
        "--output-tag",
        output_tag,
        "--exe",
        exe,
    ]
    return run_cmd(cmd, cwd=repo_root)


def make_paper_outputs(repo_root: Path, input_root: Path, tag: str, output_base: str) -> float:
    cmd = [
        sys.executable,
        "analysis/make_paper_outputs.py",
        "--input",
        str(input_root),
        "--output",
        output_base,
        "--tag",
        tag,
    ]
    return run_cmd(cmd, cwd=repo_root)


def list_pdfs(folder: Path) -> Dict[str, Path]:
    return {p.name: p for p in folder.glob("*.pdf") if p.is_file()}


def merge_side_by_side(left_path: Path, right_path: Path, out_path: Path) -> None:
    left_reader = PdfReader(str(left_path))
    right_reader = PdfReader(str(right_path))
    writer = PdfWriter()
    max_pages = max(len(left_reader.pages), len(right_reader.pages))
    for idx in range(max_pages):
        left_page = left_reader.pages[idx] if idx < len(left_reader.pages) else None
        right_page = right_reader.pages[idx] if idx < len(right_reader.pages) else None
        if left_page is None and right_page is None:
            continue
        if left_page is None:
            left_page = PageObject.create_blank_page(
                width=right_page.mediabox.width,
                height=right_page.mediabox.height,
            )
        if right_page is None:
            right_page = PageObject.create_blank_page(
                width=left_page.mediabox.width,
                height=left_page.mediabox.height,
            )
        left_w = float(left_page.mediabox.width)
        left_h = float(left_page.mediabox.height)
        right_w = float(right_page.mediabox.width)
        right_h = float(right_page.mediabox.height)
        merged = PageObject.create_blank_page(width=left_w + right_w, height=max(left_h, right_h))
        merged.merge_page(left_page)
        try:
            merged.merge_translated_page(right_page, left_w, 0, expand=False)
        except AttributeError:
            import copy

            right_copy = copy.copy(right_page)
            right_copy.add_transformation(Transformation().translate(tx=left_w, ty=0))
            merged.merge_page(right_copy)
        writer.add_page(merged)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("wb") as handle:
        writer.write(handle)


def write_compare_readme(compare_dir: Path, left_label: str, right_label: str) -> None:
    readme = compare_dir / "README.txt"
    readme.write_text(
        "Side-by-side plots.\n"
        f"Left: {left_label}\n"
        f"Right: {right_label}\n",
        encoding="ascii",
    )


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run 1M events with no-geometry and tessellated CAD, then compare plots."
    )
    parser.add_argument("--events", type=int, default=1_000_000, help="Total events per run.")
    parser.add_argument("--jobs", type=int, default=default_jobs(), help="Parallel jobs to launch.")
    parser.add_argument(
        "--exe",
        default="build/vs2022/main/Release/b02_executable.exe",
        help="Path to b02_executable.exe",
    )
    parser.add_argument(
        "--cad-file",
        default="main/assets/cad/new_assembly_cleaned_assimp.dae",
        help="CAD file to use for tessellated run.",
    )
    parser.add_argument(
        "--tessellated-tag",
        default="tierB_paper_default_1M_tessellated",
        help="Tag for tessellated outputs.",
    )
    parser.add_argument(
        "--nogeom-tag",
        default="tierB_paper_default_1M_nogeom",
        help="Tag for no-geometry outputs.",
    )
    parser.add_argument(
        "--compare-tag",
        default="compare_tessellated_vs_nogeom_1M",
        help="Tag for side-by-side outputs.",
    )
    parser.add_argument(
        "--output-base",
        default="paper_outputs",
        help="Base directory for paper outputs.",
    )
    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parent.parent
    exe = args.exe
    cad_file = args.cad_file

    print("[info] Running no-geometry first, then tessellated.")
    run_parallel(repo_root, exe, args.events, args.jobs, "none", cad_file, args.nogeom_tag)
    run_parallel(repo_root, exe, args.events, args.jobs, "tessellated", cad_file, args.tessellated_tag)

    nogeom_root = repo_root / "build/vs2022/main/Release" / f"{args.nogeom_tag}.root"
    tess_root = repo_root / "build/vs2022/main/Release" / f"{args.tessellated_tag}.root"
    if not nogeom_root.exists():
        raise SystemExit(f"Missing expected ROOT file: {nogeom_root}")
    if not tess_root.exists():
        raise SystemExit(f"Missing expected ROOT file: {tess_root}")

    make_paper_outputs(repo_root, nogeom_root, args.nogeom_tag, args.output_base)
    make_paper_outputs(repo_root, tess_root, args.tessellated_tag, args.output_base)

    output_base = repo_root / args.output_base
    left_dir = output_base / args.tessellated_tag
    right_dir = output_base / args.nogeom_tag
    compare_dir = output_base / args.compare_tag

    left_pdfs = list_pdfs(left_dir)
    right_pdfs = list_pdfs(right_dir)
    common = sorted(set(left_pdfs.keys()) & set(right_pdfs.keys()))
    if not common:
        raise SystemExit("No matching PDF files found to merge.")

    compare_dir.mkdir(parents=True, exist_ok=True)
    for name in common:
        merge_side_by_side(left_pdfs[name], right_pdfs[name], compare_dir / name)
    write_compare_readme(compare_dir, args.tessellated_tag, args.nogeom_tag)

    print(f"[done] No-geom outputs: {right_dir}")
    print(f"[done] Tessellated outputs: {left_dir}")
    print(f"[done] Side-by-side outputs: {compare_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
