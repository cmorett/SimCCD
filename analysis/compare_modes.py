#!/usr/bin/env python3
"""
Compare CAD vs none ROOT outputs with uncertainty bands and ratio panels.
"""
from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

try:
    import uproot  # noqa: E402
except ImportError as exc:
    raise SystemExit("This script requires the 'uproot' package. Install via pip.") from exc

from pixelization.pixelize_helpers import EV_PER_ELECTRON  # noqa: E402


def load_events(path: Path) -> Dict[str, np.ndarray]:
    with uproot.open(path) as root:
        if "B02Evts" not in root:
            raise RuntimeError(f"B02Evts tree missing in {path}")
        tree = root["B02Evts"]
        needed = ["EdepCCD", "trackLenCCD", "thetaPri", "muonCosTheta", "nStepsCCD"]
        events: Dict[str, np.ndarray] = {}
        for name in needed:
            if name in tree:
                events[name] = tree[name].array(library="np")
        # Always include EdepCCD and trackLenCCD for size
        for name in ("EdepCCD", "trackLenCCD", "thetaPri"):
            if name not in events and name in tree:
                events[name] = tree[name].array(library="np")
        return events


def compute_cos_down(events: Dict[str, np.ndarray]) -> np.ndarray:
    if "muonCosTheta" in events:
        cos_down = np.asarray(events["muonCosTheta"])
    else:
        cos_down = np.clip(-np.cos(np.asarray(events["thetaPri"])), 0.0, 1.0)
    return np.clip(cos_down, 0.0, 1.0)


def compute_hit_mask(events: Dict[str, np.ndarray]) -> np.ndarray:
    edep = np.asarray(events["EdepCCD"])
    track_len = np.asarray(events["trackLenCCD"])
    mask = (edep > 0) | (track_len > 0)
    if "nStepsCCD" in events:
        mask |= np.asarray(events["nStepsCCD"]) > 0
    return mask


def compute_throughgoing_mask(
    events: Dict[str, np.ndarray],
    cos_down: np.ndarray,
    thickness_cm: float,
    min_charge_e: float,
) -> np.ndarray:
    track_len = np.asarray(events["trackLenCCD"])
    edep = np.asarray(events["EdepCCD"])
    hit_mask = compute_hit_mask(events)
    charge_e = edep * 1.0e9 / EV_PER_ELECTRON
    lcos = track_len * np.abs(cos_down)
    tol = 0.02 * thickness_cm
    return (
        hit_mask
        & (track_len > 0.01)
        & (charge_e >= min_charge_e)
        & (np.abs(lcos - thickness_cm) <= tol)
    )


def binomial_efficiency(n_hit: np.ndarray, n_tot: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    eff = np.zeros_like(n_tot, dtype=float)
    err = np.zeros_like(n_tot, dtype=float)
    mask = n_tot > 0
    eff[mask] = n_hit[mask] / n_tot[mask]
    err[mask] = np.sqrt(eff[mask] * (1.0 - eff[mask]) / n_tot[mask])
    return eff, err


def ratio_with_uncertainty(
    num: np.ndarray, den: np.ndarray, num_err: np.ndarray, den_err: np.ndarray
) -> Tuple[np.ndarray, np.ndarray]:
    ratio = np.full_like(num, np.nan, dtype=float)
    err = np.full_like(num, np.nan, dtype=float)
    mask = (den > 0) & (num > 0)
    ratio[mask] = num[mask] / den[mask]
    err[mask] = ratio[mask] * np.sqrt(
        (num_err[mask] / num[mask]) ** 2 + (den_err[mask] / den[mask]) ** 2
    )
    return ratio, err


def auto_range(data: np.ndarray, pct: Tuple[float, float] = (1, 99)) -> Tuple[float, float]:
    if data.size == 0:
        return (0.0, 1.0)
    lo, hi = np.percentile(data, pct)
    if not np.isfinite(lo) or not np.isfinite(hi):
        lo, hi = float(np.min(data)), float(np.max(data))
    if hi <= lo:
        hi = lo + 1.0
    return float(lo), float(hi)


def make_ratio_plot(
    x: np.ndarray,
    y_cad: np.ndarray,
    y_none: np.ndarray,
    err_cad: np.ndarray,
    err_none: np.ndarray,
    xlabel: str,
    ylabel: str,
    out_path: Path,
    title: str,
) -> None:
    ratio, ratio_err = ratio_with_uncertainty(y_cad, y_none, err_cad, err_none)
    fig, (ax, axr) = plt.subplots(
        2, 1, figsize=(7.0, 6.0), sharex=True, gridspec_kw={"height_ratios": [3, 1]}
    )
    ax.errorbar(x, y_none, yerr=err_none, fmt="o", ms=4, label="none")
    ax.errorbar(x, y_cad, yerr=err_cad, fmt="s", ms=4, label="cad")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(alpha=0.3)
    ax.legend()

    axr.axhline(1.0, color="black", linewidth=1)
    axr.errorbar(x, ratio, yerr=ratio_err, fmt="o", ms=4)
    axr.set_xlabel(xlabel)
    axr.set_ylabel("CAD/none")
    axr.grid(alpha=0.3)

    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path)
    plt.close(fig)


def make_hist_ratio_plot(
    cad_vals: np.ndarray,
    none_vals: np.ndarray,
    n_thrown_cad: int,
    n_thrown_none: int,
    bins: int,
    xlabel: str,
    ylabel: str,
    out_path: Path,
    title: str,
) -> None:
    combined = np.concatenate([cad_vals, none_vals]) if cad_vals.size and none_vals.size else (
        cad_vals if cad_vals.size else none_vals
    )
    lo, hi = auto_range(combined)
    counts_cad, edges = np.histogram(cad_vals, bins=bins, range=(lo, hi))
    counts_none, _ = np.histogram(none_vals, bins=edges)
    centers = 0.5 * (edges[:-1] + edges[1:])

    rate_cad = counts_cad / max(1, n_thrown_cad)
    rate_none = counts_none / max(1, n_thrown_none)
    err_cad = np.sqrt(counts_cad) / max(1, n_thrown_cad)
    err_none = np.sqrt(counts_none) / max(1, n_thrown_none)

    ratio, ratio_err = ratio_with_uncertainty(rate_cad, rate_none, err_cad, err_none)

    fig, (ax, axr) = plt.subplots(
        2, 1, figsize=(7.0, 6.0), sharex=True, gridspec_kw={"height_ratios": [3, 1]}
    )
    ax.step(centers, rate_none, where="mid", label="none")
    ax.step(centers, rate_cad, where="mid", label="cad")
    ax.errorbar(centers, rate_none, yerr=err_none, fmt="o", ms=3)
    ax.errorbar(centers, rate_cad, yerr=err_cad, fmt="s", ms=3)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(alpha=0.3)
    ax.legend()

    axr.axhline(1.0, color="black", linewidth=1)
    axr.errorbar(centers, ratio, yerr=ratio_err, fmt="o", ms=3)
    axr.set_xlabel(xlabel)
    axr.set_ylabel("CAD/none")
    axr.grid(alpha=0.3)

    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path)
    plt.close(fig)


def binned_mean(values: np.ndarray, bins: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    means = np.zeros(len(bins) - 1, dtype=float)
    errs = np.zeros(len(bins) - 1, dtype=float)
    for i in range(len(bins) - 1):
        mask = (values[:, 0] >= bins[i]) & (values[:, 0] < bins[i + 1])
        if not np.any(mask):
            means[i] = float("nan")
            errs[i] = float("nan")
            continue
        vals = values[mask][:, 1]
        means[i] = float(np.mean(vals))
        errs[i] = float(np.std(vals, ddof=1) / math.sqrt(vals.size)) if vals.size > 1 else 0.0
    return means, errs


def make_binned_mean_ratio_plot(
    cad_x: np.ndarray,
    cad_y: np.ndarray,
    none_x: np.ndarray,
    none_y: np.ndarray,
    bins: np.ndarray,
    xlabel: str,
    ylabel: str,
    out_path: Path,
    title: str,
) -> None:
    cad_vals = np.column_stack([cad_x, cad_y])
    none_vals = np.column_stack([none_x, none_y])
    cad_mean, cad_err = binned_mean(cad_vals, bins)
    none_mean, none_err = binned_mean(none_vals, bins)
    centers = 0.5 * (bins[:-1] + bins[1:])
    ratio, ratio_err = ratio_with_uncertainty(cad_mean, none_mean, cad_err, none_err)

    fig, (ax, axr) = plt.subplots(
        2, 1, figsize=(7.0, 6.0), sharex=True, gridspec_kw={"height_ratios": [3, 1]}
    )
    ax.errorbar(centers, none_mean, yerr=none_err, fmt="o", ms=4, label="none")
    ax.errorbar(centers, cad_mean, yerr=cad_err, fmt="s", ms=4, label="cad")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(alpha=0.3)
    ax.legend()

    axr.axhline(1.0, color="black", linewidth=1)
    axr.errorbar(centers, ratio, yerr=ratio_err, fmt="o", ms=4)
    axr.set_xlabel(xlabel)
    axr.set_ylabel("CAD/none")
    axr.grid(alpha=0.3)

    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path)
    plt.close(fig)


def compute_summary(
    events: Dict[str, np.ndarray],
    cos_down: np.ndarray,
    through_mask: np.ndarray,
) -> Dict[str, float]:
    n_events = len(events["EdepCCD"])
    hit_mask = compute_hit_mask(events)
    n_hits = int(np.sum(hit_mask))
    hit_fraction = n_hits / n_events if n_events else 0.0
    hit_err = math.sqrt(hit_fraction * (1.0 - hit_fraction) / n_events) if n_events else 0.0
    track_len = np.asarray(events["trackLenCCD"])
    edep_mev = np.asarray(events["EdepCCD"]) * 1000.0
    through_edep = edep_mev[through_mask]
    through_dedx = (
        edep_mev[through_mask] / np.clip(track_len[through_mask], 1.0e-9, None)
        if np.any(through_mask)
        else np.asarray([])
    )
    summary = {
        "n_events": float(n_events),
        "n_hits": float(n_hits),
        "hit_fraction": hit_fraction,
        "hit_fraction_err": hit_err,
        "through_edep_mean_mev": float(np.mean(through_edep)) if through_edep.size else 0.0,
        "through_edep_sem_mev": float(np.std(through_edep, ddof=1) / math.sqrt(through_edep.size))
        if through_edep.size > 1
        else 0.0,
        "through_dedx_mean_mev_per_cm": float(np.mean(through_dedx)) if through_dedx.size else 0.0,
        "through_dedx_sem_mev_per_cm": float(np.std(through_dedx, ddof=1) / math.sqrt(through_dedx.size))
        if through_dedx.size > 1
        else 0.0,
    }
    return summary


def write_comparison_summary(
    path: Path, cad: Dict[str, float], none: Dict[str, float]
) -> None:
    rows: List[Dict[str, float]] = []

    def add_metric(name: str, cad_val: float, cad_err: float, none_val: float, none_err: float) -> None:
        ratio = cad_val / none_val if none_val else float("nan")
        ratio_err = (
            ratio * math.sqrt((cad_err / cad_val) ** 2 + (none_err / none_val) ** 2)
            if cad_val and none_val
            else float("nan")
        )
        rows.append(
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

    add_metric(
        "hit_fraction",
        cad["hit_fraction"],
        cad["hit_fraction_err"],
        none["hit_fraction"],
        none["hit_fraction_err"],
    )
    add_metric(
        "through_edep_mean_mev",
        cad["through_edep_mean_mev"],
        cad["through_edep_sem_mev"],
        none["through_edep_mean_mev"],
        none["through_edep_sem_mev"],
    )
    add_metric(
        "through_dedx_mean_mev_per_cm",
        cad["through_dedx_mean_mev_per_cm"],
        cad["through_dedx_sem_mev_per_cm"],
        none["through_dedx_mean_mev_per_cm"],
        none["through_dedx_sem_mev_per_cm"],
    )

    header = "metric,cad_value,cad_err,none_value,none_err,ratio,ratio_err"
    lines = [header]
    for row in rows:
        lines.append(
            f"{row['metric']},{row['cad_value']:.6g},{row['cad_err']:.6g},"
            f"{row['none_value']:.6g},{row['none_err']:.6g},"
            f"{row['ratio']:.6g},{row['ratio_err']:.6g}"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="Compare CAD vs none ROOT outputs.")
    parser.add_argument("--cad", required=True, help="CAD merged ROOT file.")
    parser.add_argument("--none", required=True, help="None merged ROOT file.")
    parser.add_argument("--out-dir", required=True, help="Output directory for plots.")
    parser.add_argument("--min-charge-e", type=float, default=1.0e4, help="Min charge for throughgoing.")
    parser.add_argument("--thickness-microns", type=float, default=725.0, help="CCD thickness.")
    parser.add_argument("--bins", type=int, default=60, help="Bins for spectra.")
    args = parser.parse_args()

    cad_path = Path(args.cad)
    none_path = Path(args.none)
    out_dir = Path(args.out_dir)

    cad_events = load_events(cad_path)
    none_events = load_events(none_path)

    thickness_cm = args.thickness_microns * 1.0e-4
    cad_cos = compute_cos_down(cad_events)
    none_cos = compute_cos_down(none_events)

    cad_through = compute_throughgoing_mask(cad_events, cad_cos, thickness_cm, args.min_charge_e)
    none_through = compute_throughgoing_mask(none_events, none_cos, thickness_cm, args.min_charge_e)

    # Hit efficiency vs cos(zenith)
    bins_cos = np.linspace(0.0, 1.0, 11)
    cad_tot, _ = np.histogram(cad_cos, bins=bins_cos)
    cad_hit, _ = np.histogram(cad_cos[compute_hit_mask(cad_events)], bins=bins_cos)
    none_tot, _ = np.histogram(none_cos, bins=bins_cos)
    none_hit, _ = np.histogram(none_cos[compute_hit_mask(none_events)], bins=bins_cos)
    cad_eff, cad_err = binomial_efficiency(cad_hit, cad_tot)
    none_eff, none_err = binomial_efficiency(none_hit, none_tot)
    centers_cos = 0.5 * (bins_cos[:-1] + bins_cos[1:])

    make_ratio_plot(
        centers_cos,
        cad_eff,
        none_eff,
        cad_err,
        none_err,
        xlabel="cos(zenith) down",
        ylabel="Hit efficiency",
        out_path=out_dir / "cad_vs_none_hit_efficiency_vs_coszen.pdf",
        title="Hit efficiency vs cos(zenith)",
    )

    # Edep CCD throughgoing
    cad_edep_mev = np.asarray(cad_events["EdepCCD"])[cad_through] * 1000.0
    none_edep_mev = np.asarray(none_events["EdepCCD"])[none_through] * 1000.0
    make_hist_ratio_plot(
        cad_edep_mev,
        none_edep_mev,
        n_thrown_cad=len(cad_events["EdepCCD"]),
        n_thrown_none=len(none_events["EdepCCD"]),
        bins=args.bins,
        xlabel="Edep CCD [MeV]",
        ylabel="Rate per thrown",
        out_path=out_dir / "cad_vs_none_edep_ccd_throughgoing.pdf",
        title="Throughgoing Edep CCD",
    )

    # dEdx throughgoing
    cad_track = np.asarray(cad_events["trackLenCCD"])[cad_through]
    none_track = np.asarray(none_events["trackLenCCD"])[none_through]
    cad_dedx = cad_edep_mev / np.clip(cad_track, 1.0e-9, None)
    none_dedx = none_edep_mev / np.clip(none_track, 1.0e-9, None)
    make_hist_ratio_plot(
        cad_dedx,
        none_dedx,
        n_thrown_cad=len(cad_events["EdepCCD"]),
        n_thrown_none=len(none_events["EdepCCD"]),
        bins=args.bins,
        xlabel="dE/dx [MeV/cm]",
        ylabel="Rate per thrown",
        out_path=out_dir / "cad_vs_none_dedx_throughgoing.pdf",
        title="Throughgoing dE/dx",
    )

    # Charge vs coszen (throughgoing)
    cad_charge = np.asarray(cad_events["EdepCCD"])[cad_through] * 1.0e9 / EV_PER_ELECTRON
    none_charge = np.asarray(none_events["EdepCCD"])[none_through] * 1.0e9 / EV_PER_ELECTRON
    make_binned_mean_ratio_plot(
        cad_cos[cad_through],
        cad_charge,
        none_cos[none_through],
        none_charge,
        bins=bins_cos,
        xlabel="cos(zenith) down",
        ylabel="Mean charge [e-]",
        out_path=out_dir / "cad_vs_none_charge_vs_coszen_throughgoing.pdf",
        title="Throughgoing charge vs cos(zenith)",
    )

    # Summary CSV for markdown synthesis
    cad_summary = compute_summary(cad_events, cad_cos, cad_through)
    none_summary = compute_summary(none_events, none_cos, none_through)
    write_comparison_summary(out_dir / "comparison_summary.csv", cad_summary, none_summary)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
