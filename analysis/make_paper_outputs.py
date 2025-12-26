#!/usr/bin/env python3
"""
Generate paper-ready validation plots/tables and offline pixelized examples.

Usage:
  python analysis/make_paper_outputs.py --input build/vs2022/main/Release/B02ntuples.root --tag demo --examples 100 --dist-events 20000
"""

import argparse
import math
import sys
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

try:
    import uproot  # noqa: E402
except ImportError as exc:
    raise SystemExit(
        "This script requires the 'uproot' package. Please install it in your environment "
        "(pip install uproot)."
    ) from exc

try:
    import pandas as pd  # noqa: E402
except ImportError:
    pd = None

# Ensure repo root is on path so pixelization helpers import cleanly when run as a script.
REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from pixelization.pixelize_helpers import (  # noqa: E402
    CCDParams,
    DEFAULT_CCD_THICKNESS_MICRONS,
    DEFAULT_PIXEL_SIZE_MICRONS,
    EV_PER_ELECTRON,
    build_track_image,
    ensure_dir,
    image_moments,
    save_json,
)


def load_tree(path: Path):
    f = uproot.open(path)
    if "B02Evts" not in f:
        raise RuntimeError(f"B02Evts tree missing in {path}")
    return f["B02Evts"]


def arr(tree, name):
    return tree[name].array(library="np")


def get_field(events: Dict[str, np.ndarray], name: str, default=None):
    return events[name] if name in events else default


def basic_stats(values: np.ndarray) -> Dict[str, float]:
    vals = np.asarray(values)
    if vals.size == 0:
        return {k: 0.0 for k in ["mean", "median", "rms", "p05", "p16", "p50", "p84", "p95"]}
    return {
        "mean": float(np.mean(vals)),
        "median": float(np.median(vals)),
        "rms": float(np.sqrt(np.mean(vals ** 2))),
        "p05": float(np.percentile(vals, 5)),
        "p16": float(np.percentile(vals, 16)),
        "p50": float(np.percentile(vals, 50)),
        "p84": float(np.percentile(vals, 84)),
        "p95": float(np.percentile(vals, 95)),
    }


def percentile_limits(data: np.ndarray, p_lo: float = 1.0, p_hi: float = 99.0):
    arr_np = np.asarray(data)
    if arr_np.size == 0:
        return None
    return np.percentile(arr_np, [p_lo, p_hi])


def make_hist1d(ax, datasets, bins, xlabel, ylabel="Events", title=None, logy=False, percentile_clip=None):
    combined = []
    plotted = 0
    for vals, label in datasets:
        vals_np = np.asarray(vals)
        if vals_np.size == 0:
            continue
        combined.append(vals_np)
        ax.hist(vals_np, bins=bins, histtype="step", lw=2, label=label)
        plotted += 1
    if plotted == 0:
        ax.text(0.5, 0.5, "No events", transform=ax.transAxes, ha="center", va="center")
        return
    if percentile_clip and combined:
        merged = np.concatenate(combined)
        lo, hi = np.percentile(merged, percentile_clip)
        ax.set_xlim(lo, hi)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title)
    if logy:
        ax.set_yscale("log")
    if plotted > 1:
        ax.legend()


def make_2d(ax, x, y, bins, xlabel, ylabel, title=None, percentile_clip: Tuple[float, float] = (1, 99)):
    x_np = np.asarray(x)
    y_np = np.asarray(y)
    if x_np.size == 0 or y_np.size == 0:
        ax.text(0.5, 0.5, "No events", transform=ax.transAxes, ha="center", va="center")
        return
    h = ax.hist2d(x_np, y_np, bins=bins, cmap="viridis")
    plt.colorbar(h[3], ax=ax)
    if percentile_clip:
        xlim = percentile_limits(x_np, *percentile_clip)
        ylim = percentile_limits(y_np, *percentile_clip)
        if xlim is not None:
            ax.set_xlim(xlim[0], xlim[1])
        if ylim is not None:
            ax.set_ylim(ylim[0], ylim[1])
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title)


def write_csv(path: Path, rows: List[Dict[str, float]]):
    if pd:
        pd.DataFrame(rows).to_csv(path, index=False)
    else:
        import csv

        if not rows:
            return
        with path.open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)


def select_diverse_indices(idxs: np.ndarray, values: np.ndarray, n: int) -> np.ndarray:
    if n >= idxs.size:
        return idxs
    order = np.argsort(values)
    sorted_idxs = idxs[order]
    quantiles = np.linspace(0, 1, n, endpoint=False) + 0.5 / n
    picks = []
    for q in quantiles:
        pos = int(q * (sorted_idxs.size - 1))
        picks.append(sorted_idxs[pos])
    return np.array(picks, dtype=int)


def compute_pixel_metric(
    events: Dict[str, np.ndarray],
    idx: int,
    params: CCDParams,
    canvas_mode: str,
    canvas_size,
    margin_pix: int,
):
    track_len_cm = float(events["trackLenCCD"][idx])
    edep_gev = float(events["EdepCCD"][idx])
    theta = float(events["thetaPri"][idx])
    phi = float(events["phiPri"][idx])
    img = build_track_image(
        edep_gev,
        track_len_cm,
        theta,
        phi,
        (events["xEntryCCD"][idx], events["yEntryCCD"][idx], events["zEntryCCD"][idx]),
        params,
        canvas_mode=canvas_mode,
        canvas_size=canvas_size,
        margin_pix=margin_pix,
    )
    peak = float(np.max(img))
    threshold = peak * 0.02 if peak > 0 else 0.0
    metrics = image_moments(
        img,
        threshold=threshold,
        pixel_size_microns=params.pixel_size_microns,
        track_len_cm=track_len_cm,
        edge_margin=1,
    )
    # track_len_cm is 3D chord length (cm); PCA length is 2D in the image plane using pixel_size_microns.
    length_expected_cm = expected_l2d_cm_for_event(events, idx, track_len_cm, theta)
    length_expected_pix = float("nan")
    if math.isfinite(length_expected_cm) and length_expected_cm >= 0:
        length_expected_pix = length_expected_cm * 1.0e4 / params.pixel_size_microns
    length_image_cm = metrics["length_pix_img"] * params.pixel_size_microns * 1.0e-4
    metrics.update(
        {
            "event_idx": int(idx),
            "edep_GeV": edep_gev,
            "trackLen_cm": track_len_cm,
            "theta": theta,
            "phi": phi,
            "length_geom_cm": track_len_cm,
            "length_expected_cm": length_expected_cm,
            "length_expected_pix": length_expected_pix,
            "length_image_cm": length_image_cm,
            "canvas_w": int(img.shape[1]),
            "canvas_h": int(img.shape[0]),
            "peak_charge": peak,
        }
    )
    return img, metrics


def collect_pixel_metrics(
    events: Dict[str, np.ndarray],
    idxs: Sequence[int],
    params: CCDParams,
    canvas_mode: str,
    canvas_size,
    margin_pix: int,
):
    all_metrics = []
    quality_metrics = []
    for idx in idxs:
        _, metrics = compute_pixel_metric(events, idx, params, canvas_mode, canvas_size, margin_pix)
        all_metrics.append(metrics)
        if not metrics["is_truncated"]:
            quality_metrics.append(metrics)
    return all_metrics, quality_metrics


def render_examples(
    events: Dict[str, np.ndarray],
    idxs: np.ndarray,
    outdir: Path,
    params: CCDParams,
    canvas_mode: str,
    canvas_size,
    margin_pix: int,
    n_examples: int,
    rng: np.random.Generator,
):
    images_dir = ensure_dir(outdir / "images")
    chosen = []
    diverse_pool = select_diverse_indices(
        idxs, np.asarray(events["trackLenCCD"])[idxs], min(n_examples * 3, idxs.size)
    )
    for idx in rng.permutation(diverse_pool):
        img, metrics = compute_pixel_metric(events, idx, params, canvas_mode, canvas_size, margin_pix)
        if metrics["charge"] <= 0 or metrics["is_truncated"]:
            continue
        chosen.append((idx, img, metrics))
        plt.imsave(images_dir / f"event_{idx}.png", img, origin="lower", cmap="magma")
        if len(chosen) >= n_examples:
            break
    if not chosen:
        print("WARNING: No non-truncated examples available for rendering.")
        return

    display = chosen[: min(len(chosen), max(12, min(24, n_examples)))]
    n_panels = len(display)
    ncols = min(4, max(3, int(math.ceil(math.sqrt(n_panels)))))
    nrows = int(math.ceil(n_panels / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 3.0, nrows * 3.0))
    axes_flat = axes.flat if hasattr(axes, "flat") else [axes]
    for ax, (idx, img, metrics) in zip(axes_flat, display):
        ax.imshow(img, origin="lower", cmap="magma")
        ax.axis("off")
        cos_down = max(0.0, -math.cos(events["thetaPri"][idx]))
        label = (
            f"idx={idx}\ncosz={cos_down:.2f}\nEdep={events['EdepCCD'][idx]*1000:.1f} MeV\n"
            f"Lccd={events['trackLenCCD'][idx]:.2f} cm"
        )
        ax.text(0.02, 0.98, label, transform=ax.transAxes, ha="left", va="top", color="white", fontsize=8)
        bar_pix = max(10, min(int(img.shape[1] * 0.25), 60))
        bar_um = bar_pix * params.pixel_size_microns
        x0 = 4
        y0 = 4
        ax.plot([x0, x0 + bar_pix], [y0, y0], color="white", lw=2)
        ax.text(
            x0,
            y0 + max(2, bar_pix * 0.05),
            f"{bar_pix} px ~ {bar_um:.0f} um",
            color="white",
            fontsize=8,
            va="bottom",
        )
    for ax in axes_flat[n_panels:]:
        ax.axis("off")
    fig.tight_layout()
    fig.savefig(outdir / "fig_pixelized_examples.pdf", dpi=200)
    plt.close(fig)


def plot_pixel_metrics(
    metrics_all: List[Dict[str, float]],
    metrics_quality: List[Dict[str, float]],
    outdir: Path,
    quality_only: bool,
    max_canvas_pix: int,
    fail_on_warning: bool = False,
):
    def values(metrics, key):
        return np.asarray([m[key] for m in metrics if m.get(key) is not None])

    def datasets_for(key):
        all_vals = values(metrics_all, key)
        qual_vals = values(metrics_quality, key)
        if quality_only:
            return [(qual_vals, "Quality")] if qual_vals.size else [(all_vals, "All")]
        if not qual_vals.size:
            return [(all_vals, "All")]
        return [(all_vals, "All"), (qual_vals, "Quality")]

    plots = [
        ("npix", "Cluster size [pix]", "fig_cluster_size.pdf", True, (1, 99)),
        ("charge", "Total charge [e-]", "fig_cluster_charge.pdf", True, (1, 99)),
        ("sigma_trans", "Sigma_trans (PCA) [pix]", "fig_sigma_trans.pdf", False, (1, 99)),
        ("sigma_long", "Sigma_long (PCA) [pix]", "fig_sigma_long.pdf", False, (1, 99)),
        ("elongation", "Elongation (sigma_long/sigma_trans)", "fig_elongation.pdf", False, (1, 99)),
        ("sigma_x", "Sigma X axis-aligned [pix]", "fig_sigma_x_axis.pdf", False, (1, 99)),
        ("sigma_y", "Sigma Y axis-aligned [pix]", "fig_sigma_y_axis.pdf", False, (1, 99)),
    ]
    for key, xlabel, fname, logy, clip in plots:
        datasets = datasets_for(key)
        fig, ax = plt.subplots()
        make_hist1d(ax, datasets, bins=60, xlabel=xlabel, logy=logy, percentile_clip=clip)
        fig.tight_layout()
        fig.savefig(outdir / fname)
        plt.close(fig)

    # Length overlays
    def overlay_sets(all_vals, qual_vals, label_all, label_qual):
        if quality_only:
            return [(qual_vals, label_qual)] if qual_vals.size else [(all_vals, label_all)]
        if not qual_vals.size:
            return [(all_vals, label_all)]
        return [(all_vals, label_all), (qual_vals, label_qual)]

    geom_all = values(metrics_all, "length_pix_geom")
    geom_quality = values(metrics_quality, "length_pix_geom")
    fig, ax = plt.subplots()
    make_hist1d(
        ax,
        overlay_sets(geom_all, geom_quality, "Geom (all)", "Geom (quality)"),
        bins=60,
        xlabel="Length_geom [pix]",
        percentile_clip=(1, 99),
    )
    fig.tight_layout()
    fig.savefig(outdir / "fig_length_pix_geom.pdf")
    plt.close(fig)

    img_all = values(metrics_all, "length_pix_img")
    img_quality = values(metrics_quality, "length_pix_img")
    fig, ax = plt.subplots()
    make_hist1d(
        ax,
        overlay_sets(img_all, img_quality, "Image PCA (all)", "Image PCA (quality)"),
        bins=60,
        xlabel="Length_image [pix]",
        percentile_clip=(1, 99),
    )
    fig.tight_layout()
    fig.savefig(outdir / "fig_length_pix_img.pdf")
    plt.close(fig)

    # Combined definition comparison (quality-preferred)
    length_source = metrics_quality if (quality_only and metrics_quality) else (metrics_quality or metrics_all)
    geom_vals = values(length_source, "length_pix_geom")
    img_vals = values(length_source, "length_pix_img")

    def describe(name, arr):
        if arr.size == 0:
            print(f"[length] {name}: no entries")
            return {"min": 0.0, "max": 0.0, "p95": 0.0, "p995": 0.0}
        stats = {
            "min": float(np.min(arr)),
            "max": float(np.max(arr)),
            "p50": float(np.percentile(arr, 50)),
            "p90": float(np.percentile(arr, 90)),
            "p95": float(np.percentile(arr, 95)),
            "p99": float(np.percentile(arr, 99)),
            "p995": float(np.percentile(arr, 99.5)),
        }
        print(f"[length] {name}: min={stats['min']:.2f}, p50={stats['p50']:.2f}, p95={stats['p95']:.2f}, p99={stats['p99']:.2f}, p99.5={stats['p995']:.2f}, max={stats['max']:.2f}")
        return stats

    geom_stats = describe("geom_pix", geom_vals)
    img_stats = describe("img_pix", img_vals)

    hi_axis = max(geom_stats["p995"], img_stats["p995"])
    if max_canvas_pix > 0:
        hi_axis = min(hi_axis, float(max_canvas_pix))
    hi_axis = max(hi_axis, 1.0)
    lo_axis = 0.0

    overflow_geom = float(np.sum(geom_vals > hi_axis)) if geom_vals.size else 0.0
    overflow_img = float(np.sum(img_vals > hi_axis)) if img_vals.size else 0.0
    print(
        f"[length] plotting overlay to {hi_axis:.1f} px "
        f"(geom overflow {overflow_geom:.0f}, img overflow {overflow_img:.0f})"
    )

    fig, ax = plt.subplots()
    ax.hist(geom_vals, bins=100, range=(lo_axis, hi_axis), histtype="step", lw=2, label="Geometry length (3D) [pix]")
    ax.hist(img_vals, bins=100, range=(lo_axis, hi_axis), histtype="step", lw=2, label="Image PCA length (2D) [pix]")
    ax.set_xlabel("Length [pix] (legacy 3D vs 2D)")
    ax.set_ylabel("Events")
    ax.set_xlim(lo_axis, hi_axis)
    if overflow_geom > 0 or overflow_img > 0:
        ax.text(
            0.98,
            0.95,
            f"overflow: geom={int(overflow_geom)}, img={int(overflow_img)}",
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=9,
        )
    ax.legend()
    if geom_stats["p95"] > 64.0 and hi_axis <= 70.0:
        warn_msg = f"Length axis ({hi_axis:.1f} px) is too tight while geom p95={geom_stats['p95']:.1f} px"
        print(f"WARNING: {warn_msg}")
        if fail_on_warning:
            raise SystemExit(warn_msg)
    fig.tight_layout()
    fig.savefig(outdir / "fig_length_pix_legacy.pdf")
    plt.close(fig)


def plot_length_comparison(
    length_source: List[Dict[str, float]],
    events: Dict[str, np.ndarray],
    cos_down: np.ndarray,
    mask_hit: np.ndarray,
    ccd_thickness_cm: float,
    outdir: Path,
) -> Dict[str, float]:
    idxs = np.asarray([m["event_idx"] for m in length_source], dtype=int)
    length_expected_cm = np.asarray([m["length_expected_cm"] for m in length_source], dtype=float)
    length_image_cm = np.asarray([m["length_image_cm"] for m in length_source], dtype=float)

    valid = np.isfinite(length_expected_cm) & (length_expected_cm > 0) & np.isfinite(length_image_cm)
    idxs = idxs[valid]
    length_expected_cm = length_expected_cm[valid]
    length_image_cm = length_image_cm[valid]
    ratio = length_image_cm / length_expected_cm if length_expected_cm.size else np.array([])

    ratio_median = float(np.median(ratio)) if ratio.size else 0.0
    ratio_through_median = 0.0
    through_count = 0
    if idxs.size:
        if "zEntryCCD" in events and "zExitCCD" in events and ccd_thickness_cm > 0:
            z_span = np.abs(np.asarray(events["zExitCCD"])[idxs] - np.asarray(events["zEntryCCD"])[idxs])
            through_mask = z_span > 0.8 * ccd_thickness_cm
        else:
            lcos = np.asarray(events["trackLenCCD"])[idxs] * np.abs(cos_down[idxs])
            through_mask = np.abs(lcos - ccd_thickness_cm) < 0.2 * ccd_thickness_cm if ccd_thickness_cm > 0 else np.ones_like(lcos, dtype=bool)
        through_count = int(np.sum(through_mask))
        ratio_through = ratio[through_mask]
        if ratio_through.size:
            ratio_through_median = float(np.median(ratio_through))

    fig, axes = plt.subplots(1, 3, figsize=(15, 4))

    if length_expected_cm.size:
        h = axes[0].hist2d(length_expected_cm, length_image_cm, bins=60, cmap="viridis")
        plt.colorbar(h[3], ax=axes[0])
        xlim = percentile_limits(length_expected_cm, 1, 99)
        ylim = percentile_limits(length_image_cm, 1, 99)
        if xlim is not None:
            axes[0].set_xlim(xlim[0], xlim[1])
        if ylim is not None:
            axes[0].set_ylim(ylim[0], ylim[1])
        lim_max = max(axes[0].get_xlim()[1], axes[0].get_ylim()[1])
        axes[0].plot([0, lim_max], [0, lim_max], color="white", lw=1, linestyle="--")
    else:
        axes[0].text(0.5, 0.5, "No events", transform=axes[0].transAxes, ha="center", va="center")
    axes[0].set_xlabel("L2D_expected [cm]")
    axes[0].set_ylabel("L2D_image (PCA) [cm]")
    axes[0].set_title("PCA vs expected projected length")

    if ratio.size:
        ratio_xlim = percentile_limits(ratio, 1, 99)
        axes[1].hist(ratio, bins=60, histtype="step", lw=2)
        if ratio_xlim is not None:
            axes[1].set_xlim(ratio_xlim[0], ratio_xlim[1])
        axes[1].axvline(ratio_median, color="red", linestyle="--", label=f"median={ratio_median:.2f}")
        if ratio_through_median > 0:
            axes[1].axvline(
                ratio_through_median,
                color="orange",
                linestyle="--",
                label=f"through median={ratio_through_median:.2f}",
            )
        axes[1].legend()
    else:
        axes[1].text(0.5, 0.5, "No events", transform=axes[1].transAxes, ha="center", va="center")
    axes[1].set_xlabel("L2D_image / L2D_expected")
    axes[1].set_ylabel("Events")
    axes[1].set_title("Length ratio")

    l3d_hit = np.asarray(events["trackLenCCD"])[mask_hit]
    cos_hit = np.asarray(cos_down)[mask_hit]
    if l3d_hit.size:
        h = axes[2].hist2d(cos_hit, l3d_hit, bins=60, cmap="viridis")
        plt.colorbar(h[3], ax=axes[2])
        ylim = percentile_limits(l3d_hit, 1, 99)
        if ylim is not None:
            axes[2].set_ylim(ylim[0], ylim[1])
        axes[2].set_xlim(0.0, 1.0)
        if ccd_thickness_cm > 0:
            cos_vals = np.linspace(0.05, 1.0, 50)
            axes[2].plot(cos_vals, ccd_thickness_cm / cos_vals, color="white", lw=1, linestyle="--", label="thickness/|cos|")
            axes[2].legend()
    else:
        axes[2].text(0.5, 0.5, "No events", transform=axes[2].transAxes, ha="center", va="center")
    axes[2].set_xlabel("cos_{zenith}^{down}")
    axes[2].set_ylabel("L3D [cm]")
    axes[2].set_title("3D chord vs cos(zenith)")

    fig.tight_layout()
    fig.savefig(outdir / "fig_length_pix.pdf")
    plt.close(fig)

    print(
        f"[length] L2D_image/L2D_expected median={ratio_median:.3f} "
        f"(through-going median={ratio_through_median:.3f}, n_through={through_count})"
    )
    return {
        "length_ratio_median": ratio_median,
        "length_ratio_through_median": ratio_through_median,
        "length_ratio_through_count": float(through_count),
        "length_ratio_count": float(ratio.size),
    }


def impact_uniformity_stat(x, y, bins=20):
    H, _, _ = np.histogram2d(x, y, bins=bins)
    expected = np.mean(H)
    if expected <= 0:
        return 0.0
    return float(np.max(np.abs(H - expected)) / expected)


def write_units_file(path: Path, params: CCDParams):
    units_text = f"""Units and conventions:
- thetaPri/phiPri: Geant4 polar/azimuthal angles (rad); downward cos(zenith) = -cos(thetaPri).
- EevtPri: primary kinetic energy (GeV).
- muonX0/Y0/Z0 and muonXImp/YImp/ZImp: cm.
- EdepCCD: total deposited energy in CCD sensitive volume (GeV).
- trackLenCCD: 3D chord length between first and last CCD step (cm).
- xEntryCCD/xExitCCD/yEntryCCD/yExitCCD: cm; L2D_expected = sqrt((xExit-xEntry)^2 + (yExit-yEntry)^2).
- dirX/dirY/dirZ: primary momentum direction (unit vector).
- Pixelization: pixel size {params.pixel_size_microns} um, CCD thickness {params.thickness_microns} um, {params.ev_per_electron} eV per electron.
- dE/dx: computed as EdepCCD * 1000 / trackLenCCD (MeV/cm) for trackLenCCD>0. Plots optionally require trackLenCCD>0.01 cm to avoid extreme ratios from grazing hits.
- Canvas: adaptive mode sizes canvas to track length with margin; truncation flag marks charge touching the image edge.
- PCA metrics: sigma_long/sigma_trans = sqrt(eigenvalues of charge covariance); elongation = sigma_long/sigma_trans; length_pix_img = p99-p1 along PCA axis (2D image plane); length_pix_geom = trackLenCCD * 1e4 / pixel_size (3D chord in pixels); length_expected_cm/pix = entry/exit projected length.
- Quality cuts: EdepCCD>0, trackLenCCD>0, is_truncated=False."""
    path.write_text(units_text)


def dedx_mev_per_cm(edep_GeV: np.ndarray, track_len_cm: np.ndarray) -> np.ndarray:
    mask = (track_len_cm > 0) & (edep_GeV > 0)
    return edep_GeV[mask] * 1000.0 / track_len_cm[mask]


def compute_hit_mask(events: Dict[str, np.ndarray]) -> np.ndarray:
    edep = np.asarray(events["EdepCCD"])
    track_len = np.asarray(events["trackLenCCD"])
    mask = (edep > 0) | (track_len > 0)
    if "nStepsCCD" in events:
        mask |= np.asarray(events["nStepsCCD"]) > 0
    return mask


def compute_expected_l2d_cm(events: Dict[str, np.ndarray]) -> np.ndarray:
    if all(k in events for k in ("xEntryCCD", "xExitCCD", "yEntryCCD", "yExitCCD")):
        dx = np.asarray(events["xExitCCD"]) - np.asarray(events["xEntryCCD"])
        dy = np.asarray(events["yExitCCD"]) - np.asarray(events["yEntryCCD"])
        return np.sqrt(dx * dx + dy * dy)
    if "dirX" in events and "dirY" in events:
        transverse = np.sqrt(np.asarray(events["dirX"]) ** 2 + np.asarray(events["dirY"]) ** 2)
        return np.asarray(events["trackLenCCD"]) * transverse
    return np.asarray(events["trackLenCCD"]) * np.sin(np.asarray(events["thetaPri"]))


def expected_l2d_cm_for_event(
    events: Dict[str, np.ndarray], idx: int, track_len_cm: float, theta: float
) -> float:
    # Entry/exit positions are in cm; expected L2D is the x-y projected span.
    if all(k in events for k in ("xEntryCCD", "xExitCCD", "yEntryCCD", "yExitCCD")):
        dx = float(events["xExitCCD"][idx]) - float(events["xEntryCCD"][idx])
        dy = float(events["yExitCCD"][idx]) - float(events["yEntryCCD"][idx])
        return math.sqrt(dx * dx + dy * dy)
    if "dirX" in events and "dirY" in events:
        dirx = float(events["dirX"][idx])
        diry = float(events["dirY"][idx])
        return track_len_cm * math.sqrt(dirx * dirx + diry * diry)
    return track_len_cm * math.sin(theta)


def binned_efficiency(values: np.ndarray, mask: np.ndarray, bins: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    counts_all, edges = np.histogram(values, bins=bins)
    counts_hit, _ = np.histogram(values[mask], bins=edges)
    with np.errstate(divide="ignore", invalid="ignore"):
        eff = counts_hit / counts_all
    eff = np.nan_to_num(eff)
    centers = 0.5 * (edges[:-1] + edges[1:])
    return centers, eff


def binned_profile(x: np.ndarray, y: np.ndarray, bins: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    counts, edges = np.histogram(x, bins=bins)
    sums, _ = np.histogram(x, bins=edges, weights=y)
    with np.errstate(divide="ignore", invalid="ignore"):
        means = sums / counts
    centers = 0.5 * (edges[:-1] + edges[1:])
    return centers, means, counts


def main():
    parser = argparse.ArgumentParser(description="Make paper-ready outputs for SimCCD.")
    parser.add_argument("--input", required=True, help="Path to B02ntuples ROOT file")
    parser.add_argument("--output", default="paper_outputs", help="Base output directory")
    parser.add_argument("--tag", default="run", help="Tag name for outputs")
    parser.add_argument("--examples", type=int, default=200, help="Number of pixelized examples to render")
    parser.add_argument("--dist-events", default="50000", help="Events to use for pixelization metrics ('all' or int)")
    parser.add_argument("--seed", type=int, default=12345, help="Numpy seed for reproducibility")
    parser.add_argument("--pixel-size-microns", type=float, default=DEFAULT_PIXEL_SIZE_MICRONS, help="Pixel size in microns")
    parser.add_argument("--thickness-microns", type=float, default=DEFAULT_CCD_THICKNESS_MICRONS, help="CCD thickness in microns")
    parser.add_argument("--ev-per-electron", type=float, default=EV_PER_ELECTRON, help="Energy per e- (eV)")
    parser.add_argument("--canvas-mode", choices=["adaptive", "fixed"], default="adaptive", help="Canvas sizing strategy")
    parser.add_argument("--canvas-size", type=int, default=None, help="Canvas size (pix) when using fixed mode")
    parser.add_argument("--margin-pix", type=int, default=24, help="Margin in pixels around the expected track")
    parser.add_argument("--quality-only", action="store_true", help="Use only quality (non-truncated) events for plots")
    parser.add_argument("--max-canvas", type=int, default=512, help="Cap for adaptive canvas dimension (pix)")
    parser.add_argument("--fail-on-warning", action="store_true", help="Exit non-zero if internal warnings trigger")
    args = parser.parse_args()

    rng = np.random.default_rng(args.seed)
    input_path = Path(args.input)
    out_base = Path(args.output) / args.tag
    plots_dir = ensure_dir(out_base)
    tables_dir = ensure_dir(out_base / "tables")
    warnings: List[str] = []

    params = CCDParams(
        pixel_size_microns=args.pixel_size_microns,
        thickness_microns=args.thickness_microns,
        ev_per_electron=args.ev_per_electron,
        max_canvas_pix=args.max_canvas,
    )

    uproot_file = uproot.open(input_path)
    if "B02Evts" not in uproot_file:
        raise RuntimeError(f"B02Evts tree missing in {input_path}")
    tree = uproot_file["B02Evts"]
    run_info_tree = uproot_file["B02RunInfo"] if "B02RunInfo" in uproot_file else None
    events = {name: arr(tree, name) for name in tree.keys()}
    run_info = {}
    if run_info_tree:
        for name in run_info_tree.keys():
            values = run_info_tree[name].array(library="np")
            if values.size == 0:
                continue
            val = values[0]
            if isinstance(val, bytes):
                val = val.decode("utf-8")
            if hasattr(val, "item"):
                val = val.item()
            run_info[name] = val
    if not run_info:
        fallback_map = {
            "prov_seed1": "seed1",
            "prov_seed2": "seed2",
            "prov_useTimeSeed": "useTimeSeed",
            "prov_overburdenEnabled": "overburdenEnabled",
            "prov_overburdenThickness_cm": "overburdenThickness_cm",
            "prov_overburdenZTop_cm": "overburdenZTop_cm",
            "prov_ccdGammaCut_cm": "ccdGammaCut_cm",
            "prov_ccdElectronCut_cm": "ccdElectronCut_cm",
            "prov_ccdPositronCut_cm": "ccdPositronCut_cm",
            "prov_ccdMaxStep_cm": "ccdMaxStep_cm",
            "prov_ccdThickness_cm": "ccdThickness_cm",
            "prov_muonChargeRatio": "muonChargeRatio",
            "prov_gitHashCode": "gitHashCode",
            "prov_macroHashCode": "macroHashCode",
            "prov_macroPathHash": "macroPathHash",
            "prov_physicsListHash": "physicsListHash",
        }
        for src, dst in fallback_map.items():
            if src in events and len(events[src]) > 0:
                val = np.asarray(events[src])[0]
                if hasattr(val, "item"):
                    val = val.item()
                run_info[dst] = val
    n_events = len(events["EevtPri"])
    edep_ccd = np.asarray(events["EdepCCD"])
    track_len_ccd = np.asarray(events["trackLenCCD"])
    weights = None
    if "eventLivetime_s" in events:
        weights = np.asarray(events["eventLivetime_s"])
    elif "muonWeight_s" in events:
        weights = np.asarray(events["muonWeight_s"])
    total_livetime = float(np.sum(weights)) if weights is not None else 0.0
    ccd_thickness_cm = float(run_info.get("ccdThickness_cm", params.thickness_microns * 1e-4))
    lcos_median_val = None
    lcos_ratio_val = None
    mask_hit = compute_hit_mask(events)
    mask_pixel = (edep_ccd > 0) & (track_len_ccd > 0)
    n_hits = int(np.sum(mask_hit))
    l2d_expected_cm = compute_expected_l2d_cm(events)
    mask_quality = np.zeros(n_events, dtype=bool)

    if "muonCosTheta" in events:
        cos_down = np.asarray(events["muonCosTheta"])
    else:
        cos_down = np.clip(-np.cos(events["thetaPri"]), 0.0, 1.0)
    clamp_low = int(np.sum(-np.cos(events["thetaPri"]) < 0))
    clamp_high = int(np.sum(-np.cos(events["thetaPri"]) > 1))

    energy_pri = np.asarray(events["EevtPri"])
    energy_hits = energy_pri[mask_hit]

    # Generator/source plots
    fig, ax = plt.subplots()
    make_hist1d(
        ax,
        [(energy_pri, "All thrown")],
        bins=120,
        xlabel="E_{#mu} [GeV]",
        logy=True,
        title="Primary muon energy (all thrown)",
    )
    fig.tight_layout()
    fig.savefig(plots_dir / "fig_energy_spectrum.pdf")
    plt.close(fig)

    fig, ax = plt.subplots()
    make_hist1d(
        ax,
        [(energy_pri, "All thrown"), (energy_hits, "CCD hits")],
        bins=120,
        xlabel="E_{#mu} [GeV]",
        logy=True,
        title="Primary energy (all vs hits)",
    )
    fig.tight_layout()
    fig.savefig(plots_dir / "validate_energy_hits_logy.pdf")
    plt.close(fig)

    fig, ax = plt.subplots()
    make_hist1d(
        ax,
        [(energy_pri, "All thrown"), (energy_hits, "CCD hits")],
        bins=80,
        xlabel="E_{#mu} [GeV]",
        title="Primary energy (zoom)",
    )
    ax.set_xlim(1.0, 50.0)
    fig.tight_layout()
    fig.savefig(plots_dir / "validate_energy_hits_zoom.pdf")
    plt.close(fig)
    fig, ax = plt.subplots()
    make_hist1d(
        ax,
        [(cos_down, "All thrown")],
        bins=60,
        xlabel="cos_{zenith}^{down}",
        title=f"Downward cos(zenith) all thrown (clamp low={clamp_low}, high={clamp_high})",
    )
    fig.tight_layout()
    fig.savefig(plots_dir / "fig_coszenith_down.pdf")
    plt.close(fig)

    fig, ax = plt.subplots()
    ax.hist(cos_down, bins=60, histtype="step", lw=2, density=True, label="All thrown")
    ax.hist(cos_down[mask_hit], bins=60, histtype="step", lw=2, density=True, label="CCD hits")
    ax.set_xlabel("cos_{zenith}^{down}")
    ax.set_ylabel("Normalized density")
    ax.set_title("cos(zenith): all vs hits (normalized)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(plots_dir / "fig_coszen_all_vs_hits.pdf")
    plt.close(fig)

    fig, ax = plt.subplots()
    cos_bins = np.linspace(0.0, 1.0, 21)
    centers, eff = binned_efficiency(cos_down, mask_hit, cos_bins)
    ax.plot(centers, eff, marker="o", lw=2)
    ax.set_xlabel("cos_{zenith}^{down}")
    ax.set_ylabel("Hit efficiency")
    ax.set_ylim(0.0, 1.05)
    ax.set_title("CCD hit efficiency vs cos(zenith)")
    fig.tight_layout()
    fig.savefig(plots_dir / "fig_hit_efficiency_vs_coszen.pdf")
    plt.close(fig)

    muon_pdg = events.get("muonPDG")
    mu_plus_frac = None
    mu_ratio_sampled = None
    if muon_pdg is not None:
        mu_pdg_arr = np.asarray(muon_pdg)
        n_mu_plus = float(np.sum(mu_pdg_arr == -13))
        n_mu_minus = float(np.sum(mu_pdg_arr == 13))
        total_mu = n_mu_plus + n_mu_minus
        if total_mu > 0:
            mu_plus_frac = n_mu_plus / total_mu
        if n_mu_minus > 0:
            mu_ratio_sampled = n_mu_plus / n_mu_minus
        print(f"[charge] N(mu+)={n_mu_plus:.0f}, N(mu-)={n_mu_minus:.0f}, mu+/mu-={mu_ratio_sampled if mu_ratio_sampled is not None else 'nan'}")

    fig, ax = plt.subplots()
    make_2d(
        ax,
        events["muonXImp"],
        events["muonYImp"],
        bins=80,
        xlabel="x_imp [cm]",
        ylabel="y_imp [cm]",
        title="Impact plane (all thrown)",
    )
    # overlay CCD footprint +/-0.75 cm
    for x in (-0.75, 0.75):
        ax.axvline(x, color="red", linestyle="--", linewidth=1)
    for y in (-0.75, 0.75):
        ax.axhline(y, color="red", linestyle="--", linewidth=1)
    fig.tight_layout()
    fig.savefig(plots_dir / "fig_xyImpact.pdf")
    plt.close(fig)

    fig, ax = plt.subplots()
    make_hist1d(ax, [(events["muonZ0"], None)], bins=40, xlabel="z0 [cm]", title="Source height (all thrown)")
    fig.tight_layout()
    fig.savefig(plots_dir / "fig_z0.pdf")
    plt.close(fig)

    if "muonX0" in events and "muonY0" in events:
        fig, ax = plt.subplots()
        make_2d(
            ax,
            events["muonX0"],
            events["muonY0"],
            bins=80,
            xlabel="x0 [cm]",
            ylabel="y0 [cm]",
            title="Source plane (all thrown)",
        )
        fig.tight_layout()
        fig.savefig(plots_dir / "fig_xy0.pdf")
        plt.close(fig)

    # CCD summaries (hits only)
    fig, ax = plt.subplots()
    make_hist1d(
        ax,
        [(edep_ccd[mask_hit], "CCD hits")],
        bins=100,
        xlabel="Edep CCD [GeV]",
        title="CCD deposited energy (hits)",
        logy=True,
    )
    fig.tight_layout()
    fig.savefig(plots_dir / "fig_edep_ccd.pdf")
    plt.close(fig)

    fig, ax = plt.subplots()
    make_hist1d(
        ax,
        [(track_len_ccd[mask_hit], "CCD hits")],
        bins=100,
        xlabel="track length CCD [cm]",
        title="CCD track length (hits)",
    )
    fig.tight_layout()
    fig.savefig(plots_dir / "fig_trackLen_ccd.pdf")
    plt.close(fig)

    fig, ax = plt.subplots()
    valid_len = track_len_ccd[mask_hit]
    valid_edep = edep_ccd[mask_hit]
    make_2d(
        ax,
        valid_len,
        valid_edep,
        bins=60,
        xlabel="track length CCD [cm]",
        ylabel="Edep CCD [GeV]",
        title="Edep vs length (hits)",
    )
    fig.tight_layout()
    fig.savefig(plots_dir / "fig_edep_vs_trackLen.pdf")
    plt.close(fig)

    fig, ax = plt.subplots()
    make_2d(
        ax,
        cos_down[mask_hit],
        track_len_ccd[mask_hit],
        bins=60,
        xlabel="cos_{zenith}^{down}",
        ylabel="track length CCD [cm]",
        title="Incidence vs CCD length (hits)",
    )
    fig.tight_layout()
    fig.savefig(plots_dir / "fig_costheta_vs_trackLen.pdf")
    plt.close(fig)

    # L*cos(theta) sanity check vs CCD thickness
    if ccd_thickness_cm > 0.0:
        lcos = track_len_ccd[mask_hit] * np.abs(cos_down[mask_hit])
        lcos = lcos[lcos > 0]
        if lcos.size:
            lcos_median = float(np.median(lcos))
            lcos_ratio = lcos_median / ccd_thickness_cm if ccd_thickness_cm > 0 else 0.0
            lcos_median_val = lcos_median
            lcos_ratio_val = lcos_ratio
            fig, ax = plt.subplots()
            ax.hist(lcos, bins=80, histtype="step", lw=2)
            ax.axvline(ccd_thickness_cm, color="red", linestyle="--", label="CCD thickness")
            ax.set_xlabel("trackLen * |cos(theta)| [cm]")
            ax.set_ylabel("Events")
            ax.set_title("L*cos(theta) vs CCD thickness")
            ax.legend()
            fig.tight_layout()
            fig.savefig(plots_dir / "fig_Lcos_vs_thickness.pdf")
            plt.close(fig)
            if lcos_ratio < 0.5 or lcos_ratio > 1.5:
                warnings.append(
                    f"L*cos(theta) median {lcos_median:.4f} cm deviates from CCD thickness {ccd_thickness_cm:.4f} cm"
                )
        else:
            warnings.append("No valid track lengths to compute L*cos(theta).")

    # dE/dx (geometry-level, hits only)
    dedx_all = dedx_mev_per_cm(edep_ccd[mask_hit], track_len_ccd[mask_hit])
    dedx_mask = (track_len_ccd > 0.01) & mask_hit
    dedx_plot = dedx_mev_per_cm(edep_ccd[dedx_mask], track_len_ccd[dedx_mask])
    if dedx_all.size:
        dedx_p = {p: float(np.percentile(dedx_all, p)) for p in [50, 84, 95, 99, 99.9]}
        print(
            f"[dEdx] median={dedx_p[50]:.2f} MeV/cm, p84={dedx_p[84]:.2f}, p95={dedx_p[95]:.2f}, "
            f"p99={dedx_p[99]:.2f}, p99.9={dedx_p[99.9]:.2f}"
        )
        hi_zoom = min(max(dedx_p[99], 5.0), 30.0)
        hi_tail = max(dedx_p[99.9], hi_zoom)

        fig, ax = plt.subplots()
        ax.hist(dedx_all, bins=120, histtype="step", lw=2)
        ax.set_xlabel("dE/dx [MeV/cm]")
        ax.set_ylabel("Events")
        ax.set_yscale("log")
        ax.set_title("dE/dx (hits)")
        fig.tight_layout()
        fig.savefig(plots_dir / "fig_dEdx_hits.pdf")
        plt.close(fig)

        fig, ax = plt.subplots()
        ax.hist(dedx_plot, bins=100, range=(0, hi_zoom), histtype="step", lw=2, label="trackLen>0.01 cm")
        ax.set_xlabel("dE/dx [MeV/cm]")
        ax.set_ylabel("Events")
        ax.set_xlim(0, hi_zoom)
        ax.set_title("dE/dx = Edep/trackLen (hits, zoom)")
        ax.legend()
        fig.tight_layout()
        fig.savefig(plots_dir / "fig_dEdx.pdf")
        plt.close(fig)

        fig, ax = plt.subplots()
        ax.hist(dedx_plot, bins=100, range=(0, hi_tail), histtype="step", lw=2, label="trackLen>0.01 cm")
        ax.set_xlabel("dE/dx [MeV/cm]")
        ax.set_ylabel("Events")
        ax.set_xlim(0, hi_tail)
        ax.set_yscale("log")
        ax.set_title("dE/dx tail (hits, log-y)")
        ax.legend()
        fig.tight_layout()
        fig.savefig(plots_dir / "fig_dEdx_tail.pdf")
        plt.close(fig)

        if dedx_p[50] < 0.5 or dedx_p[50] > 10.0:
            warnings.append(f"dE/dx median {dedx_p[50]:.2f} MeV/cm outside expected muon range (0.5-10).")

    # Edep vs length with dE/dx guide (hits only)
    if mask_hit.any():
        fig, ax = plt.subplots()
        make_2d(
            ax,
            track_len_ccd[mask_hit],
            edep_ccd[mask_hit],
            bins=60,
            xlabel="L3D (trackLenCCD) [cm]",
            ylabel="Edep CCD [GeV]",
            title="Edep vs L3D (hits)",
        )
        if dedx_all.size:
            dedx_mean = float(np.mean(dedx_all))
            x_line = np.linspace(0.0, float(np.percentile(track_len_ccd[mask_hit], 99)), 50)
            y_line = dedx_mean * x_line / 1000.0
            ax.plot(x_line, y_line, color="white", lw=1, linestyle="--", label=f"mean dE/dx={dedx_mean:.2f} MeV/cm")
            ax.legend()
        fig.tight_layout()
        fig.savefig(plots_dir / "fig_edep_vs_len_hits.pdf")
        plt.close(fig)

    # Pixelization
    valid_idxs = np.nonzero(mask_pixel)[0]
    metrics_all: List[Dict[str, float]] = []
    metrics_quality: List[Dict[str, float]] = []
    length_stats: Dict[str, float] = {}
    hit_fraction = float(n_hits) / float(n_events) if n_events else 0.0
    hit_rate_hz = float(n_hits) / total_livetime if total_livetime > 0 else 0.0
    geom_intersect_fraction = None
    if "geomIntersectsCCD" in events:
        geom_mask = np.asarray(events["geomIntersectsCCD"]) > 0.5
        geom_intersect_fraction = float(np.sum(geom_mask)) / float(n_events) if n_events else 0.0
    cfg_emin = events.get("cfg_EminGeV_eff")
    cfg_emax = events.get("cfg_EmaxGeV_eff")
    cfg_theta = events.get("cfg_thetaMax_deg")
    cfg_plane = (
        events.get("cfg_sourcePlaneZ_cm"),
        events.get("cfg_sourcePlaneLx_cm"),
        events.get("cfg_sourcePlaneLy_cm"),
    )
    if total_livetime > 0.0:
        print(f"Hit fraction={hit_fraction:.4f}, livetime={total_livetime:.2f} s, CCD hit rate={hit_rate_hz:.3f} Hz")
    else:
        print(f"Hit fraction={hit_fraction:.4f}")
    if cfg_emin is not None and cfg_emax is not None:
        print(f"Effective energy range: Emin={float(np.mean(cfg_emin)):.3g} GeV, Emax={float(np.mean(cfg_emax)):.3g} GeV")
    if valid_idxs.size == 0:
        print("WARNING: No events with CCD deposits; skipping pixelization outputs.")
    else:
        # Examples from quality candidates
        n_examples = min(args.examples, valid_idxs.size)
        render_examples(events, valid_idxs, plots_dir, params, args.canvas_mode, args.canvas_size, args.margin_pix, n_examples, rng)

        # Metrics subset
        if args.dist_events == "all":
            metrics_idxs = valid_idxs
        else:
            n_dist = min(int(args.dist_events), valid_idxs.size)
            metrics_idxs = rng.choice(valid_idxs, size=n_dist, replace=False)
        metrics_all, metrics_quality = collect_pixel_metrics(
            events, metrics_idxs, params, args.canvas_mode, args.canvas_size, args.margin_pix
        )
        if metrics_quality:
            # mask_quality covers the subset of events used for pixel metrics.
            quality_idxs = np.asarray([m["event_idx"] for m in metrics_quality], dtype=int)
            mask_quality[quality_idxs] = True
        plot_pixel_metrics(
            metrics_all,
            metrics_quality,
            plots_dir,
            args.quality_only,
            max_canvas_pix=args.max_canvas,
            fail_on_warning=args.fail_on_warning,
        )

        length_source = metrics_quality if (args.quality_only and metrics_quality) else (metrics_quality or metrics_all)
        if length_source:
            length_stats = plot_length_comparison(length_source, events, cos_down, mask_hit, ccd_thickness_cm, plots_dir)

        metrics_source = metrics_quality if metrics_quality else metrics_all
        if metrics_source:
            metric_idxs = np.asarray([m["event_idx"] for m in metrics_source], dtype=int)
            cos_metric = cos_down[metric_idxs]
            sigma_trans = np.asarray([m["sigma_trans"] for m in metrics_source], dtype=float)
            charge = np.asarray([m["charge"] for m in metrics_source], dtype=float)

            fig, ax = plt.subplots()
            h = ax.hist2d(cos_metric, sigma_trans, bins=[40, 60], cmap="viridis")
            plt.colorbar(h[3], ax=ax)
            centers, means, counts = binned_profile(cos_metric, sigma_trans, np.linspace(0.0, 1.0, 21))
            valid = counts > 0
            if np.any(valid):
                ax.plot(centers[valid], means[valid], color="white", lw=1.5, label="Mean")
                ax.legend()
            ax.set_xlabel("cos_{zenith}^{down}")
            ax.set_ylabel("Sigma_trans (PCA) [pix]")
            ax.set_title("Cluster width vs cos(zenith)")
            fig.tight_layout()
            fig.savefig(plots_dir / "fig_width_vs_coszen.pdf")
            plt.close(fig)

            fig, ax = plt.subplots()
            h = ax.hist2d(cos_metric, charge, bins=[40, 60], cmap="viridis")
            plt.colorbar(h[3], ax=ax)
            centers, means, counts = binned_profile(cos_metric, charge, np.linspace(0.0, 1.0, 21))
            valid = counts > 0
            if np.any(valid):
                ax.plot(centers[valid], means[valid], color="white", lw=1.5, label="Mean")
                ax.legend()
            ax.set_xlabel("cos_{zenith}^{down}")
            ax.set_ylabel("Cluster charge [e-]")
            ax.set_title("Cluster charge vs cos(zenith)")
            fig.tight_layout()
            fig.savefig(plots_dir / "fig_charge_vs_coszen.pdf")
            plt.close(fig)

    dedx_mean_val = float(np.mean(dedx_all)) if dedx_all.size else 0.0
    dedx_rms_val = float(np.sqrt(np.mean(dedx_all ** 2))) if dedx_all.size else 0.0
    ratio_through_median = float(length_stats.get("length_ratio_through_median", 0.0))
    ratio_through_count = int(length_stats.get("length_ratio_through_count", 0.0))
    print(
        "[diag] "
        f"N_thrown={n_events}, N_hits={n_hits}, hit_fraction={hit_fraction:.4f}, "
        f"dEdx_mean={dedx_mean_val:.2f} MeV/cm, dEdx_rms={dedx_rms_val:.2f} MeV/cm, "
        f"L2D_ratio_through_median={ratio_through_median:.3f} (n_through={ratio_through_count})"
    )

    # Tables
    truncated_count = max(0, len(metrics_all) - len(metrics_quality))
    validation_summary = {
        "n_events": int(n_events),
        "n_thrown": int(n_events),
        "n_events_ccd": int(n_hits),
        "n_hits": int(n_hits),
        "n_pixel_metrics_all": int(len(metrics_all)),
        "n_pixel_metrics_quality": int(len(metrics_quality)),
        "n_truncated_flagged": int(truncated_count),
        "truncated_fraction": float(truncated_count / len(metrics_all)) if metrics_all else 0.0,
        "cos_clamp_low": clamp_low,
        "cos_clamp_high": clamp_high,
        "hit_fraction": hit_fraction,
        "total_livetime_s": total_livetime,
        "hit_rate_hz": hit_rate_hz,
    }
    if weights is not None and weights.size > 0:
        validation_summary["event_weight_mean_s"] = float(np.mean(weights))
    if geom_intersect_fraction is not None:
        validation_summary["geom_intersect_fraction"] = geom_intersect_fraction
    if cfg_emin is not None:
        validation_summary["cfg_EminGeV_eff"] = float(np.mean(cfg_emin))
    if cfg_emax is not None:
        validation_summary["cfg_EmaxGeV_eff"] = float(np.mean(cfg_emax))
    if cfg_theta is not None:
        validation_summary["cfg_thetaMax_deg"] = float(np.mean(cfg_theta))
    if cfg_plane[0] is not None:
        validation_summary["cfg_sourcePlaneZ_cm"] = float(np.mean(cfg_plane[0]))
    if cfg_plane[1] is not None:
        validation_summary["cfg_sourcePlaneLx_cm"] = float(np.mean(cfg_plane[1]))
    if cfg_plane[2] is not None:
        validation_summary["cfg_sourcePlaneLy_cm"] = float(np.mean(cfg_plane[2]))
    validation_summary.update({f"energy_{k}": v for k, v in basic_stats(events["EevtPri"]).items()})
    validation_summary.update({f"EdepCCD_{k}": v for k, v in basic_stats(events["EdepCCD"]).items()})
    validation_summary.update({f"trackLenCCD_{k}": v for k, v in basic_stats(events["trackLenCCD"]).items()})
    if l2d_expected_cm.size:
        validation_summary.update(
            {f"L2D_expected_hits_{k}": v for k, v in basic_stats(l2d_expected_cm[mask_hit]).items()}
        )
    validation_summary["impact_uniformity"] = impact_uniformity_stat(events["muonXImp"], events["muonYImp"])
    if "EdepOther" in events:
        validation_summary.update({f"EdepOther_{k}": v for k, v in basic_stats(events["EdepOther"]).items()})
    if "muonModeCode" in events and len(events["muonModeCode"]) > 0:
        validation_summary["muonModeCode"] = int(np.median(np.asarray(events["muonModeCode"])))
    if dedx_all.size:
        validation_summary.update({f"dEdx_{k}": v for k, v in basic_stats(dedx_all).items()})
    if mu_plus_frac is not None:
        validation_summary["mu_plus_fraction"] = float(mu_plus_frac)
    if mu_ratio_sampled is not None:
        validation_summary["mu_plus_to_minus_sampled"] = float(mu_ratio_sampled)
    if lcos_median_val is not None:
        validation_summary["Lcos_median_cm"] = float(lcos_median_val)
    if lcos_ratio_val is not None:
        validation_summary["Lcos_ratio_to_thickness"] = float(lcos_ratio_val)
    if metrics_quality:
        validation_summary.update({f"sigma_trans_{k}": v for k, v in basic_stats(np.asarray([m["sigma_trans"] for m in metrics_quality])).items()})
        validation_summary.update({f"length_pix_geom_{k}": v for k, v in basic_stats(np.asarray([m["length_pix_geom"] for m in metrics_quality])).items()})
    if length_stats:
        validation_summary.update(length_stats)
    if warnings:
        validation_summary["warnings"] = "; ".join(warnings)
    for k in [
        "gitHash",
        "gitDirty",
        "macroPath",
        "macroHash",
        "provenanceTag",
        "physicsList",
        "useTimeSeed",
        "seed1",
        "seed2",
        "muonMode",
        "fluxModel",
        "muonChargeMode",
        "muonChargeRatio",
        "overburdenEnabled",
        "overburdenThickness_cm",
        "overburdenZTop_cm",
        "overburdenMaterial",
        "ccdGammaCut_cm",
        "ccdElectronCut_cm",
        "ccdPositronCut_cm",
        "ccdMaxStep_cm",
        "ccdThickness_cm",
        "gitHashCode",
        "macroHashCode",
        "macroPathHash",
        "physicsListHash",
    ]:
        if k in run_info:
            validation_summary[k] = run_info[k]
    write_csv(tables_dir / "validation_summary.csv", [validation_summary])

    if metrics_all:
        write_csv(tables_dir / "pixel_metrics_all.csv", metrics_all)
    if metrics_quality:
        write_csv(tables_dir / "pixel_metrics_quality.csv", metrics_quality)

    write_units_file(tables_dir / "units_and_conventions.txt", params)

    if warnings:
        print("[WARN] The following validation warnings were raised:")
        for w in warnings:
            print(f"  - {w}")
        if args.fail_on_warning:
            raise SystemExit("fail_on_warning set and warnings present.")

    # Save config
    config_payload = {
        "input": str(input_path),
        "output": str(out_base),
        "examples": args.examples,
        "dist_events": args.dist_events,
        "seed": args.seed,
        "quality_only": args.quality_only,
        "canvas_mode": args.canvas_mode,
        "canvas_size": args.canvas_size,
        "margin_pix": args.margin_pix,
        "fail_on_warning": args.fail_on_warning,
        "ccd_params": vars(params),
    }
    if run_info:
        config_payload["run_info"] = run_info
    save_json(config_payload, out_base / "run_config.json")


if __name__ == "__main__":
    main()
