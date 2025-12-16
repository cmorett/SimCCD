#!/usr/bin/env python3
"""
Benchmark CAD vs primitive geometry and generate a PDF bar chart.
"""

import argparse
import re
import subprocess
import sys
import time
import tempfile
from pathlib import Path
from typing import Dict, List


def run_case(label: str, binary: Path, macro: Path, extra_args: List[str]) -> Dict[str, float]:
    cmd = [str(binary), "--no-vis", *extra_args, str(macro)]
    start_wall = time.perf_counter()
    start_cpu = time.process_time()
    result = subprocess.run(cmd, capture_output=True, text=True)
    wall = time.perf_counter() - start_wall
    cpu = time.process_time() - start_cpu
    output = (result.stdout or "") + (result.stderr or "")

    if result.returncode != 0:
        sys.stderr.write(output)
        raise SystemExit(f"{label} run failed with code {result.returncode}")

    return {"label": label, "wall": wall, "cpu": cpu, "output": output}


def make_plot(data: List[Dict[str, float]], pdf_path: Path) -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise SystemExit("matplotlib is required to generate the benchmark plot.") from exc

    labels = [d["label"] for d in data]
    wall = [d["wall"] for d in data]

    fig, ax = plt.subplots(figsize=(6, 4))
    bars = ax.bar(labels, wall, color=["#4c78a8", "#72b7b2"])
    ax.set_ylabel("Wall time (s)")
    ax.set_title("Geometry runtime comparison")
    ax.grid(axis="y", linestyle="--", alpha=0.4)

    for bar, value in zip(bars, wall):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height(), f"{value:.2f}s",
                ha="center", va="bottom", fontsize=9)

    fig.tight_layout()
    pdf_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(pdf_path, format="pdf")
    plt.close(fig)


def main() -> None:
    repo_root = Path(__file__).resolve().parents[1]
    default_binary = repo_root / "main" / "build" / "Release" / "b02_executable.exe"
    cad_macro = repo_root / "main" / "macros" / "run_with_cad.mac"
    prim_macro = repo_root / "main" / "macros" / "run_without_cad.mac"
    default_pdf = repo_root / "benchmark_geometry.pdf"

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, default=default_binary, help="Path to b02_executable")
    parser.add_argument("--cad-macro", type=Path, default=cad_macro, help="Macro for CAD run")
    parser.add_argument("--primitive-macro", type=Path, default=prim_macro, help="Macro for primitive run")
    parser.add_argument("--pdf", type=Path, default=default_pdf, help="Output PDF for the plot")
    parser.add_argument("--events", type=int, default=10000, help="Number of events to run in each macro")
    args = parser.parse_args()

    binary = args.binary if args.binary.is_absolute() else (repo_root / args.binary)
    cad_macro = args.cad_macro if args.cad_macro.is_absolute() else (repo_root / args.cad_macro)
    prim_macro = args.primitive_macro if args.primitive_macro.is_absolute() else (repo_root / args.primitive_macro)
    pdf_path = args.pdf if args.pdf.is_absolute() else (repo_root / args.pdf)

    if not binary.exists():
        for cand in [
            repo_root / "main" / "build" / "b02_executable",
            repo_root / "main" / "build" / "Release" / "b02_executable",
            repo_root / "build" / "b02_executable",
        ]:
            if cand.exists():
                binary = cand
                break

    for path in (binary, cad_macro, prim_macro):
        if not path.exists():
            raise SystemExit(f"Missing dependency: {path}")

    def prepare_macro(source: Path, events: int, dest: Path) -> Path:
        text = source.read_text()
        if "/run/beamOn" in text:
            text = re.sub(r"/run/beamOn\s+\d+", f"/run/beamOn {events}", text)
        else:
            text += f"\n/run/beamOn {events}\n"
        dest.write_text(text)
        return dest

    runs = []
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)
        cad_macro_run = prepare_macro(cad_macro, args.events, tmpdir_path / "cad.mac")
        prim_macro_run = prepare_macro(prim_macro, args.events, tmpdir_path / "prim.mac")

        runs.append(run_case("cad-merged", binary, cad_macro_run, ["--geometry", "cad", "--cad-mode", "merged"]))
        runs.append(run_case("primitive", binary, prim_macro_run, ["--geometry", "primitive"]))

    print("Geometry runtime comparison:")
    print(f"{'mode':<12} {'wall (s)':>10} {'cpu (s)':>10}")
    for r in runs:
        print(f"{r['label']:<12} {r['wall']:>10.2f} {r['cpu']:>10.2f}")

    make_plot(runs, pdf_path)
    print(f"Saved plot to: {pdf_path}")


if __name__ == "__main__":
    main()
