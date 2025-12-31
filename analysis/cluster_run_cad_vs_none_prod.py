#!/usr/bin/env python3
"""
Cluster-friendly entrypoint for CAD-vs-none paired runs with sanity gate.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, List

try:
    import yaml
except ImportError as exc:
    raise SystemExit("This script requires PyYAML. Install via pip.") from exc


REPO_ROOT = Path(__file__).resolve().parent.parent


def resolve_python_executable() -> str:
    override = os.environ.get("SIMCCD_PYTHON")
    if not override:
        return sys.executable
    resolved = shutil.which(override)
    if resolved:
        return resolved
    path = Path(override).expanduser()
    if path.is_file():
        if os.name != "nt" and not os.access(path, os.X_OK):
            raise SystemExit(f"SIMCCD_PYTHON is not executable: {path}")
        return str(path)
    raise SystemExit(f"SIMCCD_PYTHON not found: {override}")


def load_yaml(path: Path) -> Dict:
    with path.open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle) or {}


def resolve_path(value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else (REPO_ROOT / path).resolve()


def run_cmd(cmd: List[str], dry_run: bool) -> None:
    if dry_run:
        print("[dry-run] " + " ".join(cmd))
        return
    subprocess.run(cmd, check=True, cwd=REPO_ROOT)


def run_paired(
    python_exe: str,
    tag: str,
    exe: Path,
    macro: Path,
    cad_file: Path,
    geometry_none: str,
    geometry_cad: str,
    thrown_none: int,
    thrown_cad: int,
    batch_size: int,
    max_workers: int,
    seed_base: int,
    out_root: Path,
    resume: bool,
    dry_run: bool,
) -> Path:
    cmd = [
        python_exe,
        "analysis/run_cad_vs_none_paired.py",
        "--tag",
        tag,
        "--exe",
        str(exe),
        "--mac",
        str(macro),
        "--geometry-none",
        geometry_none,
        "--geometry-cad",
        geometry_cad,
        "--thrown-none",
        str(thrown_none),
        "--thrown-cad",
        str(thrown_cad),
        "--batch-size",
        str(batch_size),
        "--max-workers",
        str(max_workers),
        "--seed-base",
        str(seed_base),
        "--out-root",
        str(out_root),
        "--cad-file",
        str(cad_file),
    ]
    if resume:
        cmd.append("--resume")
    if dry_run:
        cmd.append("--dry-run")
    run_cmd(cmd, dry_run)
    return out_root / "run_metadata.json"


def run_make_paper_outputs(
    python_exe: str,
    merged_root: Path,
    output_base: Path,
    tag: str,
    analysis_cfg: Dict,
    dry_run: bool,
    resume: bool,
) -> None:
    cutflow = output_base / tag / "tables" / "cutflow.csv"
    if resume and cutflow.exists():
        return
    cmd = [
        python_exe,
        "analysis/make_paper_outputs.py",
        "--input",
        str(merged_root),
        "--output",
        str(output_base),
        "--tag",
        tag,
        "--pixel-size-microns",
        str(analysis_cfg.get("pixel_size_microns", 15)),
        "--thickness-microns",
        str(analysis_cfg.get("thickness_microns", 725)),
        "--min-charge-e",
        str(analysis_cfg.get("min_charge_e", 1.0e4)),
        "--examples",
        str(analysis_cfg.get("examples", 200)),
        "--dist-events",
        str(analysis_cfg.get("dist_events", 50000)),
    ]
    if analysis_cfg.get("quality_only", True):
        cmd.append("--quality-only")
    run_cmd(cmd, dry_run)


def run_compare_modes(
    python_exe: str,
    cad_dir: Path,
    none_dir: Path,
    out_dir: Path,
    analysis_cfg: Dict,
    dry_run: bool,
    resume: bool,
    summary_path: Path | None = None,
    config_path: Path | None = None,
    tag: str | None = None,
) -> None:
    must_have = [
        out_dir / "compare_hit_efficiency_vs_coszen.pdf",
        out_dir / "compare_through_fraction_per_thrown_vs_coszen.pdf",
        out_dir / "compare_edep_core.pdf",
    ]
    if resume and all(p.exists() for p in must_have) and (summary_path is None or summary_path.exists()):
        return
    cmd = [
        python_exe,
        "analysis/compare_modes.py",
        "--cad",
        str(cad_dir),
        "--none",
        str(none_dir),
        "--out",
        str(out_dir),
        "--min-charge-e",
        str(analysis_cfg.get("min_charge_e", 1.0e4)),
        "--thickness-microns",
        str(analysis_cfg.get("thickness_microns", 725)),
    ]
    if summary_path is not None:
        cmd.extend(["--summary", str(summary_path)])
    if config_path is not None:
        cmd.extend(["--config", str(config_path)])
    if tag is not None:
        cmd.extend(["--tag", tag])
    run_cmd(cmd, dry_run)


def read_csv_first_row(path: Path) -> Dict[str, str]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            return row
    return {}


def read_cutflow(path: Path) -> Dict[str, Dict[str, str]]:
    rows = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            rows[row["stage"]] = row
    return rows


def check_sanity(
    none_root: Path,
    cad_root: Path,
    none_cutflow: Path,
    cad_cutflow: Path,
    none_validation: Path,
    cad_validation: Path,
    compare_dir: Path,
    thickness_cm: float,
) -> None:
    if not none_root.exists() or none_root.stat().st_size == 0:
        raise SystemExit(f"Sanity fail: missing none merged ROOT {none_root}")
    if not cad_root.exists() or cad_root.stat().st_size == 0:
        raise SystemExit(f"Sanity fail: missing cad merged ROOT {cad_root}")

    for label, path in [("none", none_cutflow), ("cad", cad_cutflow)]:
        if not path.exists():
            raise SystemExit(f"Sanity fail: missing cutflow {label}: {path}")
        cutflow = read_cutflow(path)
        if int(float(cutflow.get("hits", {}).get("count", 0))) <= 0:
            raise SystemExit(f"Sanity fail: no hits in {label} cutflow.")
        if int(float(cutflow.get("throughgoing", {}).get("count", 0))) <= 0:
            raise SystemExit(f"Sanity fail: no throughgoing events in {label} cutflow.")

    for label, path in [("none", none_validation), ("cad", cad_validation)]:
        if not path.exists():
            raise SystemExit(f"Sanity fail: missing validation summary {label}: {path}")
        row = read_csv_first_row(path)
        if not row:
            raise SystemExit(f"Sanity fail: empty validation summary {label}.")
        lcos_ratio = row.get("Lcos_ratio_to_thickness", "")
        if lcos_ratio:
            ratio = float(lcos_ratio)
            if not (0.98 <= ratio <= 1.02):
                raise SystemExit(f"Sanity fail: {label} Lcos ratio out of bounds: {ratio}")
        else:
            median = float(row.get("Lcos_median_cm", "0"))
            if abs(median - thickness_cm) > 0.02 * thickness_cm:
                raise SystemExit(f"Sanity fail: {label} Lcos median out of bounds: {median}")

    must_have = [
        compare_dir / "compare_hit_efficiency_vs_coszen.pdf",
        compare_dir / "compare_through_fraction_per_thrown_vs_coszen.pdf",
        compare_dir / "compare_edep_core.pdf",
    ]
    for path in must_have:
        if not path.exists():
            raise SystemExit(f"Sanity fail: missing comparison plot {path}")


def list_figures(folder: Path) -> List[str]:
    return sorted(p.name for p in folder.glob("*.pdf") if p.is_file())


def git_info() -> Dict[str, str]:
    def run_git(args: List[str]) -> str:
        try:
            return subprocess.check_output(["git"] + args, cwd=REPO_ROOT, text=True).strip()
        except Exception:
            return "unknown"

    return {
        "hash": run_git(["rev-parse", "HEAD"]),
        "dirty": "yes" if run_git(["status", "--porcelain"]) else "no",
    }


def write_summary(
    path: Path,
    tag: str,
    config: Dict,
    metadata_path: Path,
    none_cutflow: Path,
    cad_cutflow: Path,
    compare_dir: Path,
    compare_summary: Path,
) -> None:
    meta = {}
    if metadata_path.exists():
        meta = json.loads(metadata_path.read_text(encoding="utf-8"))
    git = git_info()

    def format_cutflow(cutflow_path: Path) -> Dict[str, str]:
        rows = read_cutflow(cutflow_path)
        return {stage: rows[stage]["count"] for stage in rows}

    none_cf = format_cutflow(none_cutflow)
    cad_cf = format_cutflow(cad_cutflow)

    comparison_rows = []
    if compare_summary.exists():
        with compare_summary.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            comparison_rows = list(reader)

    figures = list_figures(compare_dir)

    lines = []
    lines.append(f"# CAD vs none summary ({tag})")
    lines.append("")
    lines.append("## Provenance")
    lines.append(f"- git_hash: {git['hash']}")
    lines.append(f"- git_dirty: {git['dirty']}")
    lines.append(f"- exe: {config.get('exe', '')}")
    lines.append(f"- macro: {config.get('mac', '')}")
    lines.append(f"- cad_file: {config.get('cad_file', '')}")
    lines.append(f"- geometry_none: {config.get('geometry', {}).get('none', '')}")
    lines.append(f"- geometry_cad: {config.get('geometry', {}).get('cad', '')}")
    if meta:
        lines.append(f"- thrown_none: {meta.get('thrown_none', '')}")
        lines.append(f"- thrown_cad: {meta.get('thrown_cad', '')}")
        lines.append(f"- batch_size: {meta.get('batch_size', '')}")
        lines.append(f"- max_workers: {meta.get('max_workers', '')}")
        lines.append(f"- seed_base: {meta.get('seed_base', '')}")
    lines.append("")
    lines.append("## Runtime")
    if meta:
        lines.append(f"- none_elapsed_s: {meta.get('none', {}).get('elapsed_seconds', 0):.1f}")
        lines.append(f"- cad_elapsed_s: {meta.get('cad', {}).get('elapsed_seconds', 0):.1f}")
        lines.append(f"- none_per_batch_s: {meta.get('none', {}).get('per_batch_seconds', 0):.2f}")
        lines.append(f"- cad_per_batch_s: {meta.get('cad', {}).get('per_batch_seconds', 0):.2f}")
    lines.append("")
    lines.append("## Cutflow (counts)")
    stages = sorted(set(none_cf.keys()) | set(cad_cf.keys()))
    lines.append("| stage | none | cad |")
    lines.append("| --- | --- | --- |")
    for stage in stages:
        lines.append(f"| {stage} | {none_cf.get(stage, '0')} | {cad_cf.get(stage, '0')} |")
    lines.append("")
    lines.append("## Key effect sizes (CAD vs none)")
    if comparison_rows:
        lines.append("| metric | cad | none | ratio |")
        lines.append("| --- | --- | --- | --- |")
        for row in comparison_rows:
            cad_val = f"{row['cad_value']} ± {row['cad_err']}"
            none_val = f"{row['none_value']} ± {row['none_err']}"
            ratio_val = f"{row['ratio']} ± {row['ratio_err']}"
            lines.append(f"| {row['metric']} | {cad_val} | {none_val} | {ratio_val} |")
    else:
        lines.append("- comparison_summary.csv missing or empty")
    lines.append("")
    lines.append("## Figures")
    if figures:
        for fig in figures:
            lines.append(f"- {fig}")
    else:
        lines.append("- (none found)")
    lines.append("")
    lines.append("## Conservative interpretation (template)")
    lines.append(
        "The CAD-vs-none differences should be interpreted as acceptance effects from added "
        "material/geometry rather than changes in the generator or physics list. Use the ratio "
        "panels with their uncertainties to assess significance. If deviations are at the few-percent "
        "level, additional statistics are required before claiming any robust effect."
    )

    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="Run CAD-vs-none paired production with sanity gate.")
    parser.add_argument("--config", default="configs/cad_vs_none_prod.yaml", help="Config YAML.")
    parser.add_argument("--tag", default=None, help="Override tag.")
    parser.add_argument("--exe", default=None, help="Override executable path.")
    parser.add_argument("--mac", default=None, help="Override macro path.")
    parser.add_argument("--cad-file", default=None, help="Override CAD file.")
    parser.add_argument("--thrown-none", type=int, default=None, help="Override none thrown.")
    parser.add_argument("--thrown-cad", type=int, default=None, help="Override cad thrown.")
    parser.add_argument("--batch-size", type=int, default=None, help="Override batch size.")
    parser.add_argument("--max-workers", type=int, default=None, help="Override max workers.")
    parser.add_argument("--seed-base", type=int, default=None, help="Override seed base.")
    parser.add_argument("--resume", action="store_true", help="Resume existing outputs.")
    parser.add_argument("--dry-run", action="store_true", help="Print commands only.")
    parser.add_argument("--sanity-only", action="store_true", help="Run sanity gate only.")
    parser.add_argument("--skip-sanity", action="store_true", help="Skip sanity gate.")
    parser.add_argument("--only-compare", action="store_true", help="Regenerate analysis/compare only (no simulation).")
    args = parser.parse_args()

    config_path = resolve_path(args.config)
    config = load_yaml(config_path)
    if os.environ.get("SIMCCD_EXE"):
        config["exe"] = os.environ["SIMCCD_EXE"]
    if args.tag:
        config["tag"] = args.tag
    if args.exe:
        config["exe"] = args.exe
    if args.mac:
        config["mac"] = args.mac
    if args.cad_file:
        config["cad_file"] = args.cad_file
    if args.thrown_none is not None:
        config.setdefault("thrown", {})["none"] = args.thrown_none
    if args.thrown_cad is not None:
        config.setdefault("thrown", {})["cad"] = args.thrown_cad
    if args.batch_size is not None:
        config["batch_size"] = args.batch_size
    if args.max_workers is not None:
        config["max_workers"] = args.max_workers
    if args.seed_base is not None:
        config["seed_base"] = args.seed_base
    if args.resume:
        config["resume"] = True

    python_exe = resolve_python_executable()

    tag = config.get("tag", "cad_vs_none")
    exe = resolve_path(config["exe"])
    macro = resolve_path(config["mac"])
    cad_file = resolve_path(config["cad_file"])
    geometry_none = config.get("geometry", {}).get("none", "none")
    geometry_cad = config.get("geometry", {}).get("cad", "tessellated")
    thrown_none = int(config.get("thrown", {}).get("none", 1_000_000))
    thrown_cad = int(config.get("thrown", {}).get("cad", 500_000))
    batch_size = int(config.get("batch_size", 5000))
    max_workers = int(config.get("max_workers", 4))
    seed_base = int(config.get("seed_base", 12345))
    resume = bool(config.get("resume", False) or args.resume)
    analysis_cfg = config.get("analysis", {})

    out_root = resolve_path(config.get("out_root", f"outputs/{tag}"))
    paper_root = resolve_path(config.get("paper_outputs", f"paper_outputs/{tag}"))

    thickness_cm = float(analysis_cfg.get("thickness_microns", 725)) * 1.0e-4

    if args.only_compare:
        none_root = out_root / "none" / "merged.root"
        cad_root = out_root / "cad" / "merged.root"
        if not none_root.exists() or not cad_root.exists():
            raise SystemExit(f"--only-compare requested but merged roots missing: {none_root}, {cad_root}")
        run_make_paper_outputs(
            python_exe, none_root, paper_root, "none", analysis_cfg, args.dry_run, resume=False
        )
        run_make_paper_outputs(
            python_exe, cad_root, paper_root, "cad", analysis_cfg, args.dry_run, resume=False
        )
        compare_dir = paper_root / "compare"
        summary_path = REPO_ROOT / f"docs/cad_vs_none_summary_{tag}.md"
        run_compare_modes(
            python_exe,
            cad_dir=paper_root / "cad",
            none_dir=paper_root / "none",
            out_dir=compare_dir,
            analysis_cfg=analysis_cfg,
            dry_run=args.dry_run,
            resume=False,
            summary_path=summary_path,
            config_path=config_path,
            tag=tag,
        )
        return 0

    # Sanity gate
    if not args.skip_sanity:
        sanity_tag = f"{tag}_sanity"
        sanity_out_root = out_root / "sanity"
        sanity_paper_root = paper_root / "sanity"

        sanity_meta = run_paired(
            python_exe=python_exe,
            tag=sanity_tag,
            exe=exe,
            macro=macro,
            cad_file=cad_file,
            geometry_none=geometry_none,
            geometry_cad=geometry_cad,
            thrown_none=10_000,
            thrown_cad=10_000,
            batch_size=5000,
            max_workers=min(4, max_workers),
            seed_base=seed_base,
            out_root=sanity_out_root,
            resume=resume,
            dry_run=args.dry_run,
        )

        none_root = sanity_out_root / "none" / "merged.root"
        cad_root = sanity_out_root / "cad" / "merged.root"
        run_make_paper_outputs(
            python_exe, none_root, sanity_paper_root, "none", analysis_cfg, args.dry_run, resume
        )
        run_make_paper_outputs(
            python_exe, cad_root, sanity_paper_root, "cad", analysis_cfg, args.dry_run, resume
        )

        compare_dir = sanity_paper_root / "compare"
        run_compare_modes(
            python_exe,
            cad_dir=sanity_paper_root / "cad",
            none_dir=sanity_paper_root / "none",
            out_dir=compare_dir,
            analysis_cfg=analysis_cfg,
            dry_run=args.dry_run,
            resume=resume,
            summary_path=REPO_ROOT / f"docs/cad_vs_none_summary_{tag}_sanity.md",
            config_path=config_path,
            tag=f"{tag}_sanity",
        )

        if not args.dry_run:
            none_cutflow = sanity_paper_root / "none" / "tables" / "cutflow.csv"
            cad_cutflow = sanity_paper_root / "cad" / "tables" / "cutflow.csv"
            none_val = sanity_paper_root / "none" / "tables" / "validation_summary.csv"
            cad_val = sanity_paper_root / "cad" / "tables" / "validation_summary.csv"
            check_sanity(
                none_root,
                cad_root,
                none_cutflow,
                cad_cutflow,
                none_val,
                cad_val,
                compare_dir,
                thickness_cm,
            )

        if args.sanity_only:
            return 0

    # Production run
    prod_meta = run_paired(
        python_exe=python_exe,
        tag=tag,
        exe=exe,
        macro=macro,
        cad_file=cad_file,
        geometry_none=geometry_none,
        geometry_cad=geometry_cad,
        thrown_none=thrown_none,
        thrown_cad=thrown_cad,
        batch_size=batch_size,
        max_workers=max_workers,
        seed_base=seed_base,
        out_root=out_root,
        resume=resume,
        dry_run=args.dry_run,
    )

    none_root = out_root / "none" / "merged.root"
    cad_root = out_root / "cad" / "merged.root"
    run_make_paper_outputs(python_exe, none_root, paper_root, "none", analysis_cfg, args.dry_run, resume)
    run_make_paper_outputs(python_exe, cad_root, paper_root, "cad", analysis_cfg, args.dry_run, resume)

    compare_dir = paper_root / "compare"
    summary_path = REPO_ROOT / f"docs/cad_vs_none_summary_{tag}.md"
    run_compare_modes(
        python_exe,
        cad_dir=paper_root / "cad",
        none_dir=paper_root / "none",
        out_dir=compare_dir,
        analysis_cfg=analysis_cfg,
        dry_run=args.dry_run,
        resume=resume,
        summary_path=summary_path,
        config_path=config_path,
        tag=tag,
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
