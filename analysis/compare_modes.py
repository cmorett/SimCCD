#!/usr/bin/env python3
"""
Paper-ready comparison of two simulation modes (e.g., CAD vs none) using paper output folders.

Inputs:
- Two mode-specific output folders produced by analysis/make_paper_outputs.py (each contains tables/*.csv and run_config.json).
- Optional direct merged ROOT overrides (if run_config inputs are missing or moved).

Outputs:
- Comparison PDFs with ratio panels and uncertainties under the requested output folder.
- comparison_summary.csv with headline metrics and ratios.
- Optional markdown summary when --summary is provided.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path
import sys
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import colors  # noqa: E402
import numpy as np  # noqa: E402

try:
    import uproot  # noqa: E402
except ImportError as exc:  # pragma: no cover - dependency guard
    raise SystemExit("This script requires the 'uproot' package. Install via pip.") from exc

try:
    import pandas as pd  # noqa: E402
except ImportError:  # pragma: no cover - optional
    pd = None

try:
    import yaml  # noqa: E402
except ImportError:  # pragma: no cover - optional
    yaml = None

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from pixelization.pixelize_helpers import EV_PER_ELECTRON  # noqa: E402

THROUGHGOING_LCOS_TOL = 0.02
THROUGHGOING_MIN_L3D_CM = 0.01
COS_BINS_DEFAULT = 20
LFS_SIZE_BYTES = 2048


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def load_json(path: Path) -> Dict:
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {}


def load_run_config(mode_dir: Path) -> Dict:
    return load_json(mode_dir / "run_config.json")


def is_lfs_pointer(path: Path) -> bool:
    if not path.exists():
        return False
    try:
        if path.stat().st_size < LFS_SIZE_BYTES:
            head = path.read_text(encoding="utf-8", errors="ignore").splitlines()
            if head and "git-lfs" in head[0]:
                return True
            return "git-lfs" in path.read_text(encoding="utf-8", errors="ignore")
    except Exception:
        return False
    return False


def assert_not_lfs_pointer(path: Path) -> None:
    if not path.exists():
        return
    size = path.stat().st_size
    first_line = ""
    try:
        first_line = path.read_text(encoding="utf-8", errors="ignore").splitlines()[0]
    except Exception:
        first_line = ""
    if size < LFS_SIZE_BYTES or "git-lfs" in first_line:
        raise SystemExit(
            f"{path} looks like a Git LFS pointer (size {size} bytes). "
            "Point to real merged ROOTs or rerun merge."
        )


def load_csv_first_row(path: Path) -> Dict[str, float]:
    if not path.exists():
        return {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            parsed: Dict[str, float] = {}
            for k, v in row.items():
                try:
                    parsed[k] = float(v)
                except (TypeError, ValueError):
                    parsed[k] = v
            return parsed
    return {}


def load_validation_row(mode_dir: Path) -> Dict[str, float]:
    return load_csv_first_row(mode_dir / "tables" / "validation_summary.csv")


def load_cutflow(mode_dir: Path) -> Dict[str, int]:
    path = mode_dir / "tables" / "cutflow.csv"
    rows: Dict[str, int] = {}
    if not path.exists():
        return rows
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            try:
                rows[row["stage"]] = int(float(row.get("count", 0)))
            except Exception:
                continue
    return rows


def percentile_from_pixel_metrics(mode_dir: Path, column: str, percentile: float = 99.0) -> Optional[float]:
    path = mode_dir / "tables" / "pixel_metrics_quality.csv"
    if not path.exists():
        return None
    if pd is not None:
        try:
            series = pd.read_csv(path, usecols=[column])[column].dropna()
            if series.empty:
                return None
            return float(np.percentile(series.to_numpy(), percentile))
        except Exception:
            return None
    values: List[float] = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            if column not in row:
                break
            try:
                values.append(float(row[column]))
            except Exception:
                continue
    if not values:
        return None
    return float(np.percentile(np.asarray(values, dtype=float), percentile))


def load_events(root_path: Path, fields: Sequence[str]) -> Dict[str, np.ndarray]:
    assert_not_lfs_pointer(root_path)
    with uproot.open(root_path) as root_file:
        if "B02Evts" not in root_file:
            raise RuntimeError(f"B02Evts tree missing in {root_path}")
        tree = root_file["B02Evts"]
        events: Dict[str, np.ndarray] = {}
        for name in fields:
            if name in tree:
                events[name] = tree[name].array(library="np")
        return events


def compute_cos_down(events: Dict[str, np.ndarray]) -> np.ndarray:
    if "muonCosTheta" in events:
        cos_down = np.asarray(events["muonCosTheta"], dtype=float)
    elif "thetaPri" in events:
        cos_down = -np.cos(np.asarray(events["thetaPri"], dtype=float))
    else:
        raise RuntimeError("Neither muonCosTheta nor thetaPri found; cannot compute cos(zenith).")
    return np.clip(cos_down, 0.0, 1.0)


def compute_hit_mask(edep: np.ndarray, track_len: np.ndarray, n_steps: Optional[np.ndarray]) -> np.ndarray:
    mask = (edep > 0) | (track_len > 0)
    if n_steps is not None:
        mask |= n_steps > 0
    return mask


def compute_throughgoing_mask(
    edep: np.ndarray,
    track_len: np.ndarray,
    cos_down: np.ndarray,
    thickness_cm: float,
    min_charge_e: float,
    mask_hit: np.ndarray,
) -> np.ndarray:
    lcos = track_len * np.abs(cos_down)
    charge = edep * 1.0e9 / EV_PER_ELECTRON
    tol = THROUGHGOING_LCOS_TOL * thickness_cm
    return (
        mask_hit
        & (track_len > THROUGHGOING_MIN_L3D_CM)
        & (np.abs(lcos - thickness_cm) < tol)
        & (charge >= min_charge_e)
    )


def proportion_with_error(numerator: int, denominator: int) -> Tuple[float, float]:
    if denominator <= 0:
        return 0.0, 0.0
    p = numerator / denominator
    err = math.sqrt(p * (1.0 - p) / denominator)
    return p, err


def ratio_with_uncertainty(
    num: np.ndarray, den: np.ndarray, num_err: np.ndarray, den_err: np.ndarray
) -> Tuple[np.ndarray, np.ndarray]:
    ratio = np.divide(num, den, out=np.full_like(num, np.nan, dtype=float), where=den > 0)
    err = np.full_like(ratio, np.nan, dtype=float)
    mask = (den > 0) & (num > 0)
    if np.any(mask):
        err[mask] = ratio[mask] * np.sqrt(
            np.square(num_err[mask] / num[mask]) + np.square(den_err[mask] / den[mask])
        )
    return ratio, err


def safe_binomial(success: np.ndarray, total: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    success = np.asarray(success, dtype=float)
    total = np.asarray(total, dtype=float)
    eff = np.divide(success, total, out=np.zeros_like(total, dtype=float), where=total > 0)
    err = np.zeros_like(eff)
    mask = total > 0
    err[mask] = np.sqrt(np.clip(eff[mask] * (1.0 - eff[mask]) / total[mask], 0.0, None))
    return eff, err


def centers_from_edges(edges: np.ndarray) -> np.ndarray:
    return 0.5 * (edges[:-1] + edges[1:])


def logy_limits(vals: np.ndarray, floor: float = 1.0e-12) -> Tuple[float, float]:
    v = vals[np.isfinite(vals) & (vals > 0)]
    if v.size == 0:
        return floor, floor * 10.0
    ymin = max(np.min(v) * 0.5, floor)
    ymax = np.max(v) * 1.2
    if ymax <= ymin:
        ymax = ymin * 10.0
    return ymin, ymax


def make_ratio_series_plot(
    x: np.ndarray,
    y_none: np.ndarray,
    y_cad: np.ndarray,
    err_none: np.ndarray,
    err_cad: np.ndarray,
    xlabel: str,
    ylabel: str,
    title: str,
    out_path: Path,
    label_none: str,
    label_cad: str,
) -> None:
    ratio, ratio_err = ratio_with_uncertainty(y_cad, y_none, err_cad, err_none)
    fig, (ax, axr) = plt.subplots(
        2, 1, figsize=(7.5, 6.5), sharex=True, gridspec_kw={"height_ratios": [3, 1]}
    )
    ax.errorbar(x, y_none, yerr=err_none, fmt="o", ms=4, label=label_none)
    ax.errorbar(x, y_cad, yerr=err_cad, fmt="s", ms=4, label=label_cad)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(alpha=0.3)
    ax.legend()

    axr.axhline(1.0, color="black", linewidth=1)
    axr.errorbar(x, ratio, yerr=ratio_err, fmt="o", ms=4, color="tab:purple")
    axr.set_xlabel(xlabel)
    axr.set_ylabel("CAD/none")
    axr.grid(alpha=0.3)

    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path)
    plt.close(fig)


def make_hist_ratio_plot(
    data_none: np.ndarray,
    data_cad: np.ndarray,
    bins: np.ndarray,
    norm_none: float,
    norm_cad: float,
    xlabel: str,
    ylabel: str,
    title: str,
    out_path: Path,
    label_none: str,
    label_cad: str,
    logy: bool = False,
) -> None:
    fig, (ax, axr) = plt.subplots(
        2, 1, figsize=(7.5, 6.5), sharex=True, gridspec_kw={"height_ratios": [3, 1]}
    )
    if data_none.size == 0 and data_cad.size == 0:
        ax.text(0.5, 0.5, "No entries", transform=ax.transAxes, ha="center", va="center")
        axr.axis("off")
        fig.tight_layout()
        out_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_path)
        plt.close(fig)
        return

    counts_none, edges = np.histogram(data_none, bins=bins)
    counts_cad, _ = np.histogram(data_cad, bins=edges)
    centers = centers_from_edges(edges)

    rate_none = counts_none / max(norm_none, 1.0)
    rate_cad = counts_cad / max(norm_cad, 1.0)
    err_none = np.sqrt(counts_none) / max(norm_none, 1.0)
    err_cad = np.sqrt(counts_cad) / max(norm_cad, 1.0)

    ratio, ratio_err = ratio_with_uncertainty(rate_cad, rate_none, err_cad, err_none)

    ax.step(centers, rate_none, where="mid", label=label_none)
    ax.step(centers, rate_cad, where="mid", label=label_cad)
    ax.errorbar(centers, rate_none, yerr=err_none, fmt="o", ms=3)
    ax.errorbar(centers, rate_cad, yerr=err_cad, fmt="s", ms=3)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(alpha=0.3)
    if logy:
        y_top_values = np.concatenate([rate_none, rate_cad])
        ymin, ymax = logy_limits(y_top_values)
        pos = y_top_values[np.isfinite(y_top_values) & (y_top_values > 0)]
        data_min = float(np.min(pos)) if pos.size else float("nan")
        data_max = float(np.max(pos)) if pos.size else float("nan")
        print(f"[logy] {out_path.stem}: data_min={data_min:.3g}, data_max={data_max:.3g}, ylim=({ymin:.3g}, {ymax:.3g})")
        ax.set_yscale("log")
        ax.set_ylim(ymin, ymax)
    ax.legend()

    axr.axhline(1.0, color="black", linewidth=1)
    axr.errorbar(centers, ratio, yerr=ratio_err, fmt="o", ms=3, color="tab:purple")
    axr.set_xlabel(xlabel)
    axr.set_ylabel("CAD/none")
    axr.grid(alpha=0.3)
    ax.set_xlim(edges[0], edges[-1])

    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path)
    plt.close(fig)


def combined_percentile(values: Iterable[np.ndarray], pct: float) -> float:
    arrays = []
    for v in values:
        arr = np.asarray(v, dtype=float).ravel()
        if arr.size:
            arrays.append(arr[np.isfinite(arr)])
    if not arrays:
        return float("nan")
    merged = np.concatenate(arrays)
    if merged.size == 0:
        return float("nan")
    return float(np.percentile(merged, pct))


def find_matching_column(columns: Sequence[str], candidates: Sequence[str]) -> Optional[str]:
    lower_map = {c.lower(): c for c in columns}
    for cand in candidates:
        if cand in columns:
            return cand
    for cand in candidates:
        key = cand.lower()
        if key in lower_map:
            return lower_map[key]
    return None


def load_pixel_metrics_from_csv(
    mode_dir: Path,
) -> Optional[Tuple[Path, Dict[str, np.ndarray], int, int]]:
    candidates = [
        mode_dir / "tables" / "pixel_metrics_quality.csv",
        mode_dir / "tables" / "pixel_metrics_all.csv",
    ]
    path = next((p for p in candidates if p.exists()), None)
    if path is None:
        return None

    df = None
    rows: List[Dict[str, str]] = []
    columns: List[str] = []
    n_rows = 0
    try:
        if pd is not None:
            df = pd.read_csv(path)
            columns = list(df.columns)
            n_rows = len(df)
        else:
            with path.open("r", encoding="utf-8", newline="") as handle:
                reader = csv.DictReader(handle)
                columns = reader.fieldnames or []
                rows = list(reader)
                n_rows = len(rows)
    except Exception as exc:
        raise SystemExit(f"Failed to read {path}: {exc}") from exc

    def numeric_series(name: str) -> np.ndarray:
        if pd is not None and df is not None:
            series = pd.to_numeric(df[name], errors="coerce")
            return series.to_numpy(dtype=float)
        values: List[float] = []
        for row in rows:
            try:
                values.append(float(row.get(name, float("nan"))))
            except Exception:
                values.append(float("nan"))
        return np.asarray(values, dtype=float)

    cos_col = find_matching_column(columns, ["coszen_down", "coszen", "cos_theta", "costheta", "cos_theta_down"])
    theta_col = find_matching_column(columns, ["theta"])
    track_col = find_matching_column(columns, ["trackLenCCD_cm", "trackLen_cm", "trackLenCCD", "trackLen"])
    edep_col = find_matching_column(columns, ["edepCCD_GeV", "edep_GeV", "EdepCCD", "edep"])
    charge_col = find_matching_column(columns, ["charge", "cluster_charge", "charge_e"])
    nsteps_col = find_matching_column(columns, ["nStepsCCD", "n_steps", "nSteps"])

    missing: List[str] = []
    if not (cos_col or theta_col):
        missing.append("coszen_down/cos_theta/theta")
    if not track_col:
        missing.append("trackLenCCD_cm/trackLen_cm")
    if not (edep_col or charge_col):
        missing.append("edepCCD_GeV/edep_GeV/charge")
    if missing:
        raise SystemExit(
            f"Missing required columns in {path}: {', '.join(missing)}. "
            f"Columns ({len(columns)}): {columns}; rows={n_rows}"
        )

    cos_source = None
    if cos_col:
        cos_source = numeric_series(cos_col)
    elif theta_col:
        theta_vals = numeric_series(theta_col)
        theta_abs_max = np.nanmax(np.abs(theta_vals)) if np.any(np.isfinite(theta_vals)) else 0.0
        theta_rad = np.deg2rad(theta_vals) if theta_abs_max > (2.0 * math.pi) else theta_vals
        cos_source = np.cos(theta_rad)
    cos_down = np.clip(np.abs(cos_source) if cos_source is not None else np.array([], dtype=float), 0.0, 1.0)

    track = numeric_series(track_col) if track_col else np.array([], dtype=float)
    edep = numeric_series(edep_col) if edep_col else np.array([], dtype=float)
    charge = numeric_series(charge_col) if charge_col else None
    if charge is None:
        charge = edep * 1.0e9 / EV_PER_ELECTRON
    n_steps_arr = numeric_series(nsteps_col) if nsteps_col else None

    return (
        path,
        {"cos_down": cos_down, "track": track, "edep": edep, "charge": charge, "n_steps": n_steps_arr},
        n_rows,
        len(columns),
    )


@dataclass
class ModeData:
    label: str
    dir: Path
    root_path: Optional[Path]
    edep: np.ndarray
    track: np.ndarray
    cos_down: np.ndarray
    n_steps: Optional[np.ndarray]
    mask_hit: np.ndarray
    mask_through: np.ndarray
    charge: np.ndarray
    lcos: np.ndarray
    dedx_hits: np.ndarray
    dedx_through: np.ndarray
    n_thrown: int
    n_hits: int
    n_through: int
    charge_p99_quality: Optional[float]
    cutflow: Dict[str, int]
    validation: Dict[str, float]
    run_config: Dict


def discover_paper_dirs(
    tag: Optional[str],
    paper_none_arg: Optional[str],
    paper_cad_arg: Optional[str],
    cli_none: Optional[str],
    cli_cad: Optional[str],
) -> Tuple[Path, Path]:
    if paper_none_arg:
        paper_none = Path(paper_none_arg)
    elif cli_none:
        paper_none = Path(cli_none)
    elif tag:
        candidates = [Path("paper_outputs") / tag / "none", Path(f"paper_outputs/{tag}_none")]
        existing = [p for p in candidates if p.exists()]
        if not existing:
            raise SystemExit(f"No paper output directory found for 'none'. Tried: {candidates}")
        paper_none = existing[0]
    else:
        raise SystemExit("Paper output directory for 'none' not provided and tag is missing.")

    if paper_cad_arg:
        paper_cad = Path(paper_cad_arg)
    elif cli_cad:
        paper_cad = Path(cli_cad)
    elif tag:
        candidates = [Path("paper_outputs") / tag / "cad", Path(f"paper_outputs/{tag}_cad")]
        existing = [p for p in candidates if p.exists()]
        if not existing:
            raise SystemExit(f"No paper output directory found for 'cad'. Tried: {candidates}")
        paper_cad = existing[0]
    else:
        raise SystemExit("Paper output directory for 'cad' not provided and tag is missing.")

    return paper_none, paper_cad


def load_mode_data(
    label: str,
    mode_dir: Path,
    root_override: Optional[Path],
    thickness_cm: float,
    min_charge_e: float,
) -> ModeData:
    run_config = load_run_config(mode_dir)
    root_path = root_override
    if root_path is None and run_config.get("input"):
        root_path = Path(run_config["input"])
    if root_path is not None and not root_path.is_absolute():
        root_path = (REPO_ROOT / root_path).resolve()

    cutflow = load_cutflow(mode_dir)
    validation = load_validation_row(mode_dir)
    charge_p99_quality = percentile_from_pixel_metrics(mode_dir, "charge")

    edep = np.array([])
    track = np.array([])
    cos_down = np.array([])
    n_steps: Optional[np.ndarray] = None
    mask_hit = np.array([], dtype=bool)
    mask_through = np.array([], dtype=bool)
    charge = np.array([])
    lcos = np.array([])
    dedx_hits = np.array([])
    dedx_through = np.array([])
    n_thrown = 0
    n_hits = 0
    n_through = 0
    data_source = "none"
    csv_path: Optional[Path] = None
    csv_rows = 0
    csv_cols = 0

    csv_payload = load_pixel_metrics_from_csv(mode_dir)
    if root_path and root_path.exists():
        assert_not_lfs_pointer(root_path)
        events = load_events(
            root_path,
            ["EdepCCD", "trackLenCCD", "thetaPri", "muonCosTheta", "nStepsCCD"],
        )
        edep = np.asarray(events.get("EdepCCD", []), dtype=float)
        track = np.asarray(events.get("trackLenCCD", []), dtype=float)
        cos_down = compute_cos_down(events)
        n_steps = np.asarray(events["nStepsCCD"], dtype=float) if "nStepsCCD" in events else None
        n_thrown = len(edep)
        charge = edep * 1.0e9 / EV_PER_ELECTRON
        data_source = "root"
    else:
        if root_path:
            print(f"[WARN] Missing merged ROOT for {label}; comparisons limited. Looking for {root_path}")

    use_csv = (edep.size == 0 or track.size == 0 or cos_down.size == 0)
    if use_csv and csv_payload:
        csv_path, csv_arrays, csv_rows, csv_cols = csv_payload
        edep = csv_arrays["edep"]
        track = csv_arrays["track"]
        cos_down = csv_arrays["cos_down"]
        charge = csv_arrays["charge"]
        n_steps = csv_arrays.get("n_steps")
        data_source = "csv"
    elif not use_csv and csv_payload:
        csv_path, _, csv_rows, csv_cols = csv_payload
    elif use_csv and not csv_payload:
        print(f"[WARN] No pixel metrics CSV found for {label}; data arrays are empty.")

    edep = np.nan_to_num(np.asarray(edep, dtype=float), nan=0.0, posinf=0.0, neginf=0.0)
    track = np.nan_to_num(np.asarray(track, dtype=float), nan=0.0, posinf=0.0, neginf=0.0)
    cos_down = np.clip(np.nan_to_num(np.asarray(cos_down, dtype=float), nan=0.0, posinf=1.0, neginf=0.0), 0.0, 1.0)
    if n_steps is not None:
        n_steps = np.nan_to_num(np.asarray(n_steps, dtype=float), nan=0.0, posinf=0.0, neginf=0.0)

    mask_hit = compute_hit_mask(edep, track, n_steps)
    mask_through = compute_throughgoing_mask(edep, track, cos_down, thickness_cm, min_charge_e, mask_hit) if mask_hit.size else np.array([], dtype=bool)
    charge = charge if charge.size else edep * 1.0e9 / EV_PER_ELECTRON
    charge = np.nan_to_num(np.asarray(charge, dtype=float), nan=0.0, posinf=0.0, neginf=0.0)
    lcos = track * np.abs(cos_down)
    dedx_hits = edep[mask_hit] * 1000.0 / np.clip(track[mask_hit], 1.0e-6, None) if mask_hit.size else np.array([])
    dedx_through = (
        edep[mask_through] * 1000.0 / np.clip(track[mask_through], 1.0e-6, None) if mask_through.size else np.array([])
    )
    n_hits = int(np.sum(mask_hit))
    n_through = int(np.sum(mask_through))

    if not n_thrown:
        cf_thrown = cutflow.get("thrown") if cutflow else None
        if cf_thrown is None and validation:
            cf_thrown = validation.get("n_thrown") or validation.get("n_events") or validation.get("n_events_ccd")
        try:
            n_thrown = int(cf_thrown) if cf_thrown is not None else int(len(edep))
        except Exception:
            n_thrown = len(edep)
    if not n_hits:
        cf_hits = cutflow.get("hits") if cutflow else None
        if cf_hits is None and validation:
            cf_hits = validation.get("n_hits") or validation.get("n_events_ccd")
        try:
            n_hits = int(cf_hits) if cf_hits is not None else int(np.sum(mask_hit))
        except Exception:
            n_hits = int(np.sum(mask_hit))
    if not n_through:
        cf_through = cutflow.get("throughgoing") if cutflow else None
        if cf_through is None and validation:
            cf_through = validation.get("n_throughgoing")
        try:
            n_through = int(cf_through) if cf_through is not None else int(np.sum(mask_through))
        except Exception:
            n_through = int(np.sum(mask_through))

    if not cutflow and n_thrown:
        cutflow = {
            "thrown": n_thrown,
            "hits": n_hits,
            "throughgoing": n_through,
        }

    mode_tag = (label.lower().split()[0] if label else "mode")[:4]
    csv_label = str(csv_path) if csv_path else "none"
    print(f"[mode {mode_tag:<4}] csv={csv_label} rows={csv_rows} hit={n_hits} through={n_through} cols={csv_cols}")
    if mask_hit.size and n_hits == 0:
        print(f"[diag] {label}: hit mask removed all {mask_hit.size} entries from {data_source} data.")
    if mask_through.size and n_through == 0:
        print(f"[diag] {label}: throughgoing mask removed all entries (hit count {n_hits}).")

    return ModeData(
        label=label,
        dir=mode_dir,
        root_path=root_path,
        edep=edep,
        track=track,
        cos_down=cos_down,
        n_steps=n_steps,
        mask_hit=mask_hit,
        mask_through=mask_through,
        charge=charge,
        lcos=lcos,
        dedx_hits=dedx_hits,
        dedx_through=dedx_through,
        n_thrown=n_thrown,
        n_hits=n_hits,
        n_through=n_through,
        charge_p99_quality=charge_p99_quality,
        cutflow=cutflow,
        validation=validation,
        run_config=run_config,
    )


def build_comparison_plots(
    none_mode: ModeData,
    cad_mode: ModeData,
    out_dir: Path,
    thickness_cm: float,
    label_none: str,
    label_cad: str,
    bins_cos: int,
) -> List[Path]:
    plots: List[Path] = []
    norm_none_thrown = max(none_mode.n_thrown, 1)
    norm_cad_thrown = max(cad_mode.n_thrown, 1)
    norm_none_hits = max(none_mode.n_hits, 1)
    norm_cad_hits = max(cad_mode.n_hits, 1)

    def fallback_value(val: float, default: float) -> float:
        try:
            if val is None or not np.isfinite(val) or val <= 0:
                return default
        except Exception:
            return default
        return float(val)

    if none_mode.cos_down.size == 0:
        print(f"[diag] No cos(zenith) entries for {label_none}; cos-binned plots will be empty.")
    if cad_mode.cos_down.size == 0:
        print(f"[diag] No cos(zenith) entries for {label_cad}; cos-binned plots will be empty.")
    if none_mode.mask_hit.size == 0:
        print(f"[diag] {label_none} has no hit entries; hit histograms will be empty.")
    if cad_mode.mask_hit.size == 0:
        print(f"[diag] {label_cad} has no hit entries; hit histograms will be empty.")

    bins_cos_edges = np.linspace(0.0, 1.0, bins_cos + 1)
    centers_cos = centers_from_edges(bins_cos_edges)

    # Hit efficiency vs cos(zenith)
    none_tot, _ = np.histogram(none_mode.cos_down, bins=bins_cos_edges)
    cad_tot, _ = np.histogram(cad_mode.cos_down, bins=bins_cos_edges)
    none_hits_b, _ = np.histogram(none_mode.cos_down[none_mode.mask_hit], bins=bins_cos_edges)
    cad_hits_b, _ = np.histogram(cad_mode.cos_down[cad_mode.mask_hit], bins=bins_cos_edges)
    none_eff, none_eff_err = safe_binomial(none_hits_b, none_tot)
    cad_eff, cad_eff_err = safe_binomial(cad_hits_b, cad_tot)
    out_path = out_dir / "compare_hit_efficiency_vs_coszen.pdf"
    make_ratio_series_plot(
        centers_cos,
        none_eff,
        cad_eff,
        none_eff_err,
        cad_eff_err,
        xlabel="cos(zenith) down",
        ylabel="Hit efficiency",
        title="Hit efficiency vs cos(zenith)",
        out_path=out_path,
        label_none=label_none,
        label_cad=label_cad,
    )
    plots.append(out_path)

    # Throughgoing fraction per thrown
    none_through_b, _ = np.histogram(none_mode.cos_down[none_mode.mask_through], bins=bins_cos_edges)
    cad_through_b, _ = np.histogram(cad_mode.cos_down[cad_mode.mask_through], bins=bins_cos_edges)
    none_thr_frac, none_thr_err = safe_binomial(none_through_b, none_tot)
    cad_thr_frac, cad_thr_err = safe_binomial(cad_through_b, cad_tot)
    out_path = out_dir / "compare_through_fraction_per_thrown_vs_coszen.pdf"
    make_ratio_series_plot(
        centers_cos,
        none_thr_frac,
        cad_thr_frac,
        none_thr_err,
        cad_thr_err,
        xlabel="cos(zenith) down",
        ylabel="Throughgoing / thrown",
        title="Throughgoing fraction per thrown",
        out_path=out_path,
        label_none=label_none,
        label_cad=label_cad,
    )
    plots.append(out_path)

    # Throughgoing fraction of hits
    none_thr_hits_frac, none_thr_hits_err = safe_binomial(none_through_b, none_hits_b)
    cad_thr_hits_frac, cad_thr_hits_err = safe_binomial(cad_through_b, cad_hits_b)
    out_path = out_dir / "compare_through_fraction_of_hits_vs_coszen.pdf"
    make_ratio_series_plot(
        centers_cos,
        none_thr_hits_frac,
        cad_thr_hits_frac,
        none_thr_hits_err,
        cad_thr_hits_err,
        xlabel="cos(zenith) down",
        ylabel="Throughgoing / hits",
        title="Throughgoing fraction of hits vs cos(zenith)",
        out_path=out_path,
        label_none=label_none,
        label_cad=label_cad,
    )
    plots.append(out_path)

    # Edep CCD core/tail (per thrown)
    bins_edep_core = np.linspace(0.0, 0.0025, 100)
    out_path = out_dir / "compare_edep_core.pdf"
    make_hist_ratio_plot(
        none_mode.edep[none_mode.mask_hit],
        cad_mode.edep[cad_mode.mask_hit],
        bins=bins_edep_core,
        norm_none=norm_none_thrown,
        norm_cad=norm_cad_thrown,
        xlabel="Edep CCD [GeV]",
        ylabel="Rate per thrown",
        title="Edep CCD core (hits, per thrown)",
        out_path=out_path,
        label_none=label_none,
        label_cad=label_cad,
        logy=True,
    )
    plots.append(out_path)

    bins_edep_tail = np.linspace(0.0, 0.04, 160)
    out_path = out_dir / "compare_edep_tail.pdf"
    make_hist_ratio_plot(
        none_mode.edep[none_mode.mask_hit],
        cad_mode.edep[cad_mode.mask_hit],
        bins=bins_edep_tail,
        norm_none=norm_none_thrown,
        norm_cad=norm_cad_thrown,
        xlabel="Edep CCD [GeV]",
        ylabel="Rate per thrown",
        title="Edep CCD tail (hits, per thrown)",
        out_path=out_path,
        label_none=label_none,
        label_cad=label_cad,
        logy=True,
    )
    plots.append(out_path)

    # dE/dx core/tail (per hit)
    bins_dedx_core = np.linspace(0.0, 25.0, 120)
    out_path = out_dir / "compare_dedx_core.pdf"
    make_hist_ratio_plot(
        none_mode.dedx_hits,
        cad_mode.dedx_hits,
        bins=bins_dedx_core,
        norm_none=norm_none_hits,
        norm_cad=norm_cad_hits,
        xlabel="dE/dx [MeV/cm]",
        ylabel="Rate per hit",
        title="dE/dx core (hits, per hit)",
        out_path=out_path,
        label_none=label_none,
        label_cad=label_cad,
        logy=False,
    )
    plots.append(out_path)

    bins_dedx_tail = np.linspace(0.0, 400.0, 160)
    out_path = out_dir / "compare_dedx_tail.pdf"
    make_hist_ratio_plot(
        none_mode.dedx_hits,
        cad_mode.dedx_hits,
        bins=bins_dedx_tail,
        norm_none=norm_none_hits,
        norm_cad=norm_cad_hits,
        xlabel="dE/dx [MeV/cm]",
        ylabel="Rate per hit",
        title="dE/dx tail (hits, per hit)",
        out_path=out_path,
        label_none=label_none,
        label_cad=label_cad,
        logy=True,
    )
    plots.append(out_path)

    # Lcos distribution (hits, per thrown)
    lcos_pct = combined_percentile([none_mode.lcos[none_mode.mask_hit], cad_mode.lcos[cad_mode.mask_hit]], 99.5)
    default_lcos_hi = thickness_cm * 1.2 if thickness_cm > 0 else 0.1
    lcos_combined_hi = fallback_value(lcos_pct, default_lcos_hi)
    lcos_combined_hi = fallback_value(lcos_combined_hi, max(default_lcos_hi, 0.1))
    bins_lcos = np.linspace(0.0, max(lcos_combined_hi, default_lcos_hi, 0.1), 140)
    out_path = out_dir / "compare_lcos_distribution.pdf"
    make_hist_ratio_plot(
        none_mode.lcos[none_mode.mask_hit],
        cad_mode.lcos[cad_mode.mask_hit],
        bins=bins_lcos,
        norm_none=norm_none_thrown,
        norm_cad=norm_cad_thrown,
        xlabel="L*cos(zenith) [cm]",
        ylabel="Rate per thrown",
        title="Lcos distribution (hits, per thrown)",
        out_path=out_path,
        label_none=label_none,
        label_cad=label_cad,
        logy=True,
    )
    plots.append(out_path)

    # Charge distributions (hits + throughgoing)
    charge_hits_none = none_mode.charge[none_mode.mask_hit]
    charge_hits_cad = cad_mode.charge[cad_mode.mask_hit]
    charge_through_none = none_mode.charge[none_mode.mask_through]
    charge_through_cad = cad_mode.charge[cad_mode.mask_through]

    charge_core_pct = combined_percentile([charge_hits_none, charge_hits_cad], 99.0)
    charge_core_hi = fallback_value(charge_core_pct, 1.0)
    bins_charge_core = np.linspace(0.0, max(charge_core_hi, 1.0), 140)
    out_path = out_dir / "compare_charge_hits_core.pdf"
    make_hist_ratio_plot(
        charge_hits_none,
        charge_hits_cad,
        bins=bins_charge_core,
        norm_none=norm_none_thrown,
        norm_cad=norm_cad_thrown,
        xlabel="Cluster charge [e-]",
        ylabel="Rate per thrown",
        title="Charge distribution (hits, core)",
        out_path=out_path,
        label_none=label_none,
        label_cad=label_cad,
        logy=True,
    )
    plots.append(out_path)

    charge_tail_pct = combined_percentile([charge_hits_none, charge_hits_cad], 99.9)
    charge_tail_hi = fallback_value(
        charge_tail_pct,
        charge_core_hi if np.isfinite(charge_core_hi) and charge_core_hi > 0 else 1.0,
    )
    bins_charge_tail = np.linspace(0.0, max(charge_tail_hi, charge_core_hi, 1.0), 160)
    out_path = out_dir / "compare_charge_hits_tail.pdf"
    make_hist_ratio_plot(
        charge_hits_none,
        charge_hits_cad,
        bins=bins_charge_tail,
        norm_none=norm_none_thrown,
        norm_cad=norm_cad_thrown,
        xlabel="Cluster charge [e-]",
        ylabel="Rate per thrown",
        title="Charge distribution (hits, tail)",
        out_path=out_path,
        label_none=label_none,
        label_cad=label_cad,
        logy=True,
    )
    plots.append(out_path)

    if charge_through_none.size or charge_through_cad.size:
        charge_thr_core_pct = combined_percentile([charge_through_none, charge_through_cad], 99.0)
        charge_thr_core_hi = fallback_value(charge_thr_core_pct, 1.0)
        bins_charge_thr_core = np.linspace(0.0, max(charge_thr_core_hi, 1.0), 120)
        out_path = out_dir / "compare_charge_through_core.pdf"
        make_hist_ratio_plot(
            charge_through_none,
            charge_through_cad,
            bins=bins_charge_thr_core,
            norm_none=norm_none_thrown,
            norm_cad=norm_cad_thrown,
            xlabel="Cluster charge [e-]",
            ylabel="Rate per thrown",
            title="Charge distribution (throughgoing, core)",
            out_path=out_path,
            label_none=label_none,
            label_cad=label_cad,
            logy=True,
        )
        plots.append(out_path)

        charge_thr_tail_pct = combined_percentile([charge_through_none, charge_through_cad], 99.9)
        charge_thr_tail_hi = fallback_value(
            charge_thr_tail_pct,
            charge_thr_core_hi if np.isfinite(charge_thr_core_hi) and charge_thr_core_hi > 0 else 1.0,
        )
        bins_charge_thr_tail = np.linspace(0.0, max(charge_thr_tail_hi, charge_thr_core_hi, 1.0), 150)
        out_path = out_dir / "compare_charge_through_tail.pdf"
        make_hist_ratio_plot(
            charge_through_none,
            charge_through_cad,
            bins=bins_charge_thr_tail,
            norm_none=norm_none_thrown,
            norm_cad=norm_cad_thrown,
            xlabel="Cluster charge [e-]",
            ylabel="Rate per thrown",
            title="Charge distribution (throughgoing, tail)",
            out_path=out_path,
            label_none=label_none,
            label_cad=label_cad,
            logy=True,
        )
        plots.append(out_path)

    return plots


def write_comparison_summary(path: Path, metrics: List[Dict[str, float]]) -> None:
    fieldnames = ["metric", "cad_value", "cad_err", "none_value", "none_err", "ratio", "ratio_err"]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in metrics:
            writer.writerow({k: row.get(k, "") for k in fieldnames})


def build_metrics(
    none_mode: ModeData,
    cad_mode: ModeData,
    charge_p99_hits_none: float,
    charge_p99_hits_cad: float,
    charge_p99_through_none: float,
    charge_p99_through_cad: float,
) -> List[Dict[str, float]]:
    metrics: List[Dict[str, float]] = []

    def add_metric(name: str, none_val: float, none_err: float, cad_val: float, cad_err: float) -> None:
        ratio = cad_val / none_val if none_val else math.nan
        ratio_err = math.nan
        if cad_val and none_val:
            ratio_err = ratio * math.sqrt(
                (cad_err / cad_val) ** 2 + (none_err / none_val) ** 2 if cad_err or none_err else 0.0
            )
        metrics.append(
            {
                "metric": name,
                "cad_value": cad_val,
                "cad_err": cad_err,
                "none_value": none_val,
                "none_err": none_err,
                "ratio": ratio,
                "ratio_err": ratio_err,
            }
        )

    none_hit_frac, none_hit_err = proportion_with_error(none_mode.n_hits, none_mode.n_thrown)
    cad_hit_frac, cad_hit_err = proportion_with_error(cad_mode.n_hits, cad_mode.n_thrown)
    add_metric("hits_per_thrown", none_hit_frac, none_hit_err, cad_hit_frac, cad_hit_err)

    none_thr_frac, none_thr_err = proportion_with_error(none_mode.n_through, none_mode.n_thrown)
    cad_thr_frac, cad_thr_err = proportion_with_error(cad_mode.n_through, cad_mode.n_thrown)
    add_metric("through_per_thrown", none_thr_frac, none_thr_err, cad_thr_frac, cad_thr_err)

    none_thr_hit_frac, none_thr_hit_err = proportion_with_error(none_mode.n_through, max(none_mode.n_hits, 1))
    cad_thr_hit_frac, cad_thr_hit_err = proportion_with_error(cad_mode.n_through, max(cad_mode.n_hits, 1))
    add_metric("through_fraction_of_hits", none_thr_hit_frac, none_thr_hit_err, cad_thr_hit_frac, cad_thr_hit_err)

    add_metric("charge_p99_hits", charge_p99_hits_none, 0.0, charge_p99_hits_cad, 0.0)
    add_metric("charge_p99_through", charge_p99_through_none, 0.0, charge_p99_through_cad, 0.0)

    return metrics


def write_markdown_summary(
    path: Path,
    tag: str,
    config_info: Dict,
    none_mode: ModeData,
    cad_mode: ModeData,
    metrics: List[Dict[str, float]],
    compare_plots: List[Path],
) -> None:
    thrown_none = cad_mode.cutflow.get("thrown", None)  # type: ignore
    thrown_cad = cad_mode.cutflow.get("thrown", None)
    thrown_none = none_mode.cutflow.get("thrown", thrown_none)
    thrown_cad = cad_mode.cutflow.get("thrown", thrown_cad)

    lines: List[str] = []
    lines.append(f"# CAD vs none summary ({tag})")
    lines.append("")
    lines.append("## Run config")
    compare_base = compare_plots[0].parent if compare_plots else path.parent
    cad_file = config_info.get("cad_file") or cad_mode.run_config.get("run_info", {}).get("macroPath", "")
    geom_none = config_info.get("geometry", {}).get("none", "")
    geom_cad = config_info.get("geometry", {}).get("cad", "")
    lines.append(f"- tag: {tag}")
    if thrown_none is not None:
        lines.append(f"- thrown_none: {thrown_none}")
    if thrown_cad is not None:
        lines.append(f"- thrown_cad: {thrown_cad}")
    if geom_none:
        lines.append(f"- geometry_none: {geom_none}")
    if geom_cad:
        lines.append(f"- geometry_cad: {geom_cad}")
    if cad_file:
        lines.append(f"- cad_file: {cad_file}")
    lines.append(f"- compare_dir: {compare_base}")
    lines.append("")

    lines.append("## Cutflow (counts)")
    lines.append("| stage | none | cad |")
    lines.append("| --- | --- | --- |")
    stages = sorted(set(none_mode.cutflow.keys()) | set(cad_mode.cutflow.keys()))
    for stage in stages:
        lines.append(
            f"| {stage} | {none_mode.cutflow.get(stage, 0)} | {cad_mode.cutflow.get(stage, 0)} |"
        )
    lines.append("")

    lines.append("## Headline effects")
    lines.append("| metric | none | cad | ratio |")
    lines.append("| --- | --- | --- | --- |")
    for row in metrics:
        none_val = row["none_value"]
        cad_val = row["cad_value"]
        ratio = row["ratio"]
        none_err = row["none_err"]
        cad_err = row["cad_err"]
        ratio_err = row["ratio_err"]

        def fmt(val: float, err: float) -> str:
            if isinstance(val, float) and (math.isnan(val) or math.isinf(val)):
                return "nan"
            return f"{val:.6g} +/- {err:.2g}"

        lines.append(
            f"| {row['metric']} | {fmt(none_val, none_err)} | {fmt(cad_val, cad_err)} | {fmt(ratio, ratio_err)} |"
        )
    lines.append("")

    lines.append("## Key plots")
    for plot in sorted(compare_plots):
        rel = plot
        if plot.is_absolute():
            try:
                rel = plot.relative_to(path.parent.parent)
            except ValueError:
                rel = plot
        lines.append(f"- `{rel}`")
    lines.append("")

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compare CAD vs none paper outputs with ratio panels.")
    parser.add_argument("--cad", help="CAD mode output folder (from make_paper_outputs).")
    parser.add_argument("--none", dest="none", help="None mode output folder.")
    parser.add_argument("--paper-cad", dest="paper_cad", help="Explicit paper output folder for CAD mode.")
    parser.add_argument("--paper-none", dest="paper_none", help="Explicit paper output folder for none mode.")
    parser.add_argument("--tag", default=None, help="Tag name for autodiscovery (paper_outputs/<tag>_{none,cad}).")
    parser.add_argument(
        "--out",
        default=None,
        help="Output folder for comparison plots (default paper_outputs/<tag>/compare when tag is set).",
    )
    parser.add_argument("--cad-root", dest="cad_root", default=None, help="Optional override merged ROOT for CAD.")
    parser.add_argument("--none-root", dest="none_root", default=None, help="Optional override merged ROOT for none.")
    parser.add_argument("--min-charge-e", type=float, default=1.0e4, help="Min charge for throughgoing selection.")
    parser.add_argument("--thickness-microns", type=float, default=725.0, help="CCD thickness in microns.")
    parser.add_argument("--label-none", default="none", help="Legend label for none mode.")
    parser.add_argument("--label-cad", default="CAD (steel+copper)", help="Legend label for CAD mode.")
    parser.add_argument("--bins-cos", type=int, default=COS_BINS_DEFAULT, help="Bins for cos(zenith) profiles.")
    parser.add_argument("--summary", default=None, help="Optional markdown summary output path.")
    parser.add_argument("--config", default=None, help="Optional YAML config to include in summary.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    paper_none_dir, paper_cad_dir = discover_paper_dirs(
        args.tag,
        args.paper_none,
        args.paper_cad,
        args.none,
        args.cad,
    )
    tag = args.tag
    if tag is None:
        if args.out:
            tag = Path(args.out).parent.name
        else:
            tag = paper_none_dir.parent.name if paper_none_dir.name == "none" else paper_none_dir.name.replace("_none", "")

    if args.out:
        out_dir = ensure_dir(Path(args.out))
    elif tag:
        out_dir = ensure_dir(Path("paper_outputs") / tag / "compare")
    else:
        raise SystemExit("No output directory specified and tag is missing; provide --out or --tag.")
    thickness_cm = args.thickness_microns * 1.0e-4

    none_mode = load_mode_data(
        args.label_none,
        paper_none_dir,
        Path(args.none_root) if args.none_root else None,
        thickness_cm,
        args.min_charge_e,
    )
    cad_mode = load_mode_data(
        args.label_cad,
        paper_cad_dir,
        Path(args.cad_root) if args.cad_root else None,
        thickness_cm,
        args.min_charge_e,
    )

    compare_plots = build_comparison_plots(
        none_mode,
        cad_mode,
        out_dir,
        thickness_cm,
        args.label_none,
        args.label_cad,
        bins_cos=args.bins_cos,
    )

    # Tail metrics from pixel metrics if available; fallback to event charges.
    charge_p99_hits_none = (
        none_mode.charge_p99_quality
        if none_mode.charge_p99_quality is not None
        else combined_percentile([none_mode.charge[none_mode.mask_hit]], 99.0)
    )
    charge_p99_hits_cad = (
        cad_mode.charge_p99_quality
        if cad_mode.charge_p99_quality is not None
        else combined_percentile([cad_mode.charge[cad_mode.mask_hit]], 99.0)
    )
    charge_p99_through_none = combined_percentile([none_mode.charge[none_mode.mask_through]], 99.0)
    charge_p99_through_cad = combined_percentile([cad_mode.charge[cad_mode.mask_through]], 99.0)

    metrics = build_metrics(
        none_mode,
        cad_mode,
        charge_p99_hits_none,
        charge_p99_hits_cad,
        charge_p99_through_none,
        charge_p99_through_cad,
    )
    comparison_summary_path = out_dir / "comparison_summary.csv"
    write_comparison_summary(comparison_summary_path, metrics)

    summary_path: Optional[Path] = None
    if args.summary:
        summary_path = Path(args.summary)
    elif tag:
        summary_path = Path("docs") / f"cad_vs_none_summary_{tag}.md"

    if summary_path:
        config_info: Dict = {}
        if args.config and Path(args.config).exists() and yaml:
            try:
                config_info = yaml.safe_load(Path(args.config).read_text(encoding="utf-8")) or {}
            except Exception:
                config_info = {}
        write_markdown_summary(
            summary_path,
            tag or summary_path.stem,
            config_info,
            none_mode,
            cad_mode,
            metrics,
            compare_plots,
        )

    manifest_lines = [
        f"paper_none={paper_none_dir}",
        f"paper_cad={paper_cad_dir}",
        f"out_dir={out_dir}",
        f"comparison_summary={comparison_summary_path}",
        f"summary_markdown={summary_path if summary_path else 'none'}",
        "outputs:",
    ]
    for p in sorted(compare_plots + [comparison_summary_path] + ([summary_path] if summary_path else [])):
        manifest_lines.append(str(p))
    (out_dir / "compare_manifest.txt").write_text("\n".join(manifest_lines) + "\n", encoding="utf-8")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
