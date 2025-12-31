#!/usr/bin/env python3
"""
Package per-mode paper outputs and comparisons into a single zip bundle.

Usage:
  python analysis/package_paper_outputs.py --tag cad_vs_none_prod_v1
  python analysis/package_paper_outputs.py --tag cad_vs_none_prod_v1 \
      --paper-none paper_outputs/cad_vs_none_prod_v1_none \
      --paper-cad paper_outputs/cad_vs_none_prod_v1_cad \
      --compare paper_outputs/cad_vs_none_prod_v1/compare \
      --summary docs/cad_vs_none_summary_cad_vs_none_prod_v1.md
"""
from __future__ import annotations

import argparse
from pathlib import Path
import zipfile


def discover_dir(tag: str, mode: str, override: str | None) -> Path:
    if override:
        return Path(override)
    candidates = [Path("paper_outputs") / tag / mode, Path(f"paper_outputs/{tag}_{mode}")]
    for cand in candidates:
        if cand.exists():
            return cand
    raise SystemExit(f"Could not find paper output directory for {mode}. Tried: {candidates}")


def main() -> int:
    parser = argparse.ArgumentParser(description="Package paper outputs + compare into a zip bundle.")
    parser.add_argument("--tag", required=True, help="Tag name (e.g., cad_vs_none_prod_v1).")
    parser.add_argument("--paper-none", dest="paper_none", default=None, help="Paper outputs dir for none mode.")
    parser.add_argument("--paper-cad", dest="paper_cad", default=None, help="Paper outputs dir for cad mode.")
    parser.add_argument("--compare", default=None, help="Compare output dir (default paper_outputs/<tag>/compare).")
    parser.add_argument("--summary", default=None, help="Summary markdown path (default docs/cad_vs_none_summary_<tag>.md).")
    parser.add_argument("--output", default=None, help="Output zip path (default paper_outputs/<tag>_bundle.zip).")
    args = parser.parse_args()

    paper_none = discover_dir(args.tag, "none", args.paper_none)
    paper_cad = discover_dir(args.tag, "cad", args.paper_cad)
    compare_dir = Path(args.compare) if args.compare else Path("paper_outputs") / args.tag / "compare"
    if not compare_dir.exists():
        raise SystemExit(f"Compare directory missing: {compare_dir}")

    summary_path = Path(args.summary) if args.summary else Path("docs") / f"cad_vs_none_summary_{args.tag}.md"
    comparison_summary = compare_dir / "comparison_summary.csv"
    manifest = compare_dir / "compare_manifest.txt"

    output_zip = Path(args.output) if args.output else Path("paper_outputs") / f"{args.tag}_bundle.zip"
    output_zip.parent.mkdir(parents=True, exist_ok=True)

    members = [
        (paper_none, f"{paper_none.name}/"),
        (paper_cad, f"{paper_cad.name}/"),
        (compare_dir, f"{compare_dir.name}/"),
    ]
    with zipfile.ZipFile(output_zip, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for src_dir, arc_prefix in members:
            for path in src_dir.rglob("*"):
                arcname = Path(arc_prefix) / path.relative_to(src_dir)
                zf.write(path, arcname)
        if summary_path.exists():
            zf.write(summary_path, summary_path.name)
        if comparison_summary.exists():
            zf.write(comparison_summary, comparison_summary.name)
        if manifest.exists():
            zf.write(manifest, manifest.name)

    print(f"Wrote bundle: {output_zip}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
