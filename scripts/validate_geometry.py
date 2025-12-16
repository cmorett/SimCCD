#!/usr/bin/env python3
"""
Run a short geometry validation and fail if overlap/navigation warnings appear.
"""

import argparse
import re
import subprocess
import sys
from pathlib import Path


FORBIDDEN_PATTERNS = [
    r"GeomNav1002",
    r"GeomNav0003",
    r"GeomSolids1001",
    r"GeomVol1002",
    r"Defects in solid",
    r"wrong orientation",
    r"Overlap with volume",
    r"Stuck Track",
]


def run_validation(binary: Path, macro: Path, cad_mode: str) -> str:
    cmd = [
        str(binary),
        "--no-vis",
        "--geometry",
        "cad",
        "--cad-mode",
        cad_mode,
        str(macro),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    output = (result.stdout or "") + (result.stderr or "")
    if result.returncode != 0:
        sys.stderr.write(output)
        raise SystemExit(f"Validation run failed with code {result.returncode}")
    return output


def main() -> None:
    repo_root = Path(__file__).resolve().parents[1]
    default_binary = repo_root / "main" / "build" / "Release" / "b02_executable.exe"
    default_macro = repo_root / "main" / "macros" / "geom_check.mac"

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, default=default_binary, help="Path to b02_executable")
    parser.add_argument("--macro", type=Path, default=default_macro, help="Macro to drive the geometry test")
    parser.add_argument("--cad-mode", default="merged", help="CAD mode to validate (merged, tessellated, parts)")
    args = parser.parse_args()

    binary = args.binary if args.binary.is_absolute() else (repo_root / args.binary)
    macro = args.macro if args.macro.is_absolute() else (repo_root / args.macro)

    if not binary.exists():
        candidates = [
            repo_root / "main" / "build" / "b02_executable",
            repo_root / "main" / "build" / "Release" / "b02_executable",
            repo_root / "build" / "b02_executable",
        ]
        for cand in candidates:
            if cand.exists():
                binary = cand
                break

    if not binary.exists():
        raise SystemExit(f"Binary not found: {binary}")
    if not macro.exists():
        raise SystemExit(f"Macro not found: {macro}")

    output = run_validation(binary, macro, args.cad_mode)

    problems = []
    for pattern in FORBIDDEN_PATTERNS:
        if re.search(pattern, output, flags=re.IGNORECASE):
            problems.append(pattern)

    if problems:
        sys.stderr.write("Geometry validation failed; detected patterns:\n")
        for pat in problems:
            sys.stderr.write(f"  - {pat}\n")
        raise SystemExit(1)

    print("Geometry validation passed with no overlap/navigation warnings.")


if __name__ == "__main__":
    main()
