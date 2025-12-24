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
    img = build_track_image(
        events["EdepCCD"][idx],
        events["trackLenCCD"][idx],
        events["thetaPri"][idx],
        events["phiPri"][idx],
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
        track_len_cm=float(events["trackLenCCD"][idx]),
        edge_margin=1,
    )
    metrics.update(
        {
            "event_idx": int(idx),
            "edep_GeV": float(events["EdepCCD"][idx]),
            "trackLen_cm": float(events["trackLenCCD"][idx]),
            "theta": float(events["thetaPri"][idx]),
            "phi": float(events["phiPri"][idx]),
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
    ax.hist(geom_vals, bins=100, range=(lo_axis, hi_axis), histtype="step", lw=2, label="Geometry-based length [pix]")
    ax.hist(img_vals, bins=100, range=(lo_axis, hi_axis), histtype="step", lw=2, label="Image PCA length [pix]")
    ax.set_xlabel("Length [pix] (geom vs PCA)")
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
    fig.savefig(outdir / "fig_length_pix.pdf")
    plt.close(fig)


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
- trackLenCCD: chord length between first and last CCD step (cm).
- dirX/dirY/dirZ: primary momentum direction (unit vector).
- Pixelization: pixel size {params.pixel_size_microns} um, CCD thickness {params.thickness_microns} um, {params.ev_per_electron} eV per electron.
- dE/dx: computed as EdepCCD * 1000 / trackLenCCD (MeV/cm) for trackLenCCD>0. Plots optionally require trackLenCCD>0.01 cm to avoid extreme ratios from grazing hits.
- Canvas: adaptive mode sizes canvas to track length with margin; truncation flag marks charge touching the image edge.
- PCA metrics: sigma_long/sigma_trans = sqrt(eigenvalues of charge covariance); elongation = sigma_long/sigma_trans; length_pix_img = p99-p1 along PCA axis; length_pix_geom = trackLenCCD * 1e4 / pixel_size.
- Quality cuts: EdepCCD>0, trackLenCCD>0, is_truncated=False."""
    path.write_text(units_text)


def dedx_mev_per_cm(edep_GeV: np.ndarray, track_len_cm: np.ndarray) -> np.ndarray:
    mask = (track_len_cm > 0) & (edep_GeV > 0)
    return edep_GeV[mask] * 1000.0 / track_len_cm[mask]


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

    params = CCDParams(
        pixel_size_microns=args.pixel_size_microns,
        thickness_microns=args.thickness_microns,
        ev_per_electron=args.ev_per_electron,
        max_canvas_pix=args.max_canvas,
    )

    tree = load_tree(input_path)
    events = {name: arr(tree, name) for name in tree.keys()}
    n_events = len(events["EevtPri"])

    # Generator/source plots
    fig, ax = plt.subplots()
    make_hist1d(ax, [(events["EevtPri"], None)], bins=120, xlabel="E_{#mu} [GeV]", logy=True, title="Primary muon energy")
    fig.tight_layout()
    fig.savefig(plots_dir / "fig_energy_spectrum.pdf")
    plt.close(fig)

    cos_down = np.clip(-np.cos(events["thetaPri"]), 0.0, 1.0)
    clamp_low = int(np.sum(-np.cos(events["thetaPri"]) < 0))
    clamp_high = int(np.sum(-np.cos(events["thetaPri"]) > 1))
    fig, ax = plt.subplots()
    make_hist1d(
        ax,
        [(cos_down, None)],
        bins=60,
        xlabel="cos_{zenith}^{down}",
        title=f"Downward cos(zenith) (clamp low={clamp_low}, high={clamp_high})",
    )
    fig.tight_layout()
    fig.savefig(plots_dir / "fig_coszenith_down.pdf")
    plt.close(fig)

    fig, ax = plt.subplots()
    make_2d(ax, events["muonXImp"], events["muonYImp"], bins=80, xlabel="x_imp [cm]", ylabel="y_imp [cm]", title="Impact plane")
    fig.tight_layout()
    fig.savefig(plots_dir / "fig_xyImpact.pdf")
    plt.close(fig)

    fig, ax = plt.subplots()
    make_hist1d(ax, [(events["muonZ0"], None)], bins=40, xlabel="z0 [cm]", title="Source height")
    fig.tight_layout()
    fig.savefig(plots_dir / "fig_z0.pdf")
    plt.close(fig)

    # CCD summaries
    fig, ax = plt.subplots()
    make_hist1d(ax, [(events["EdepCCD"], None)], bins=100, xlabel="Edep CCD [GeV]", title="CCD deposited energy", logy=True)
    fig.tight_layout()
    fig.savefig(plots_dir / "fig_edep_ccd.pdf")
    plt.close(fig)

    fig, ax = plt.subplots()
    make_hist1d(ax, [(events["trackLenCCD"], None)], bins=100, xlabel="track length CCD [cm]", title="CCD track length")
    fig.tight_layout()
    fig.savefig(plots_dir / "fig_trackLen_ccd.pdf")
    plt.close(fig)

    fig, ax = plt.subplots()
    valid_len = events["trackLenCCD"]
    valid_edep = events["EdepCCD"]
    make_2d(
        ax,
        valid_len,
        valid_edep,
        bins=60,
        xlabel="track length CCD [cm]",
        ylabel="Edep CCD [GeV]",
        title="Edep vs length",
    )
    fig.tight_layout()
    fig.savefig(plots_dir / "fig_edep_vs_trackLen.pdf")
    plt.close(fig)

    fig, ax = plt.subplots()
    make_2d(
        ax,
        cos_down,
        events["trackLenCCD"],
        bins=60,
        xlabel="cos_{zenith}^{down}",
        ylabel="track length CCD [cm]",
        title="Incidence vs CCD length",
    )
    fig.tight_layout()
    fig.savefig(plots_dir / "fig_costheta_vs_trackLen.pdf")
    plt.close(fig)

    # dE/dx (geometry-level)
    dedx_all = dedx_mev_per_cm(events["EdepCCD"], events["trackLenCCD"])
    dedx_mask = events["trackLenCCD"] > 0.01
    dedx_plot = dedx_mev_per_cm(events["EdepCCD"][dedx_mask], events["trackLenCCD"][dedx_mask])
    if dedx_all.size:
        dedx_p = {p: float(np.percentile(dedx_all, p)) for p in [50, 84, 95, 99, 99.9]}
        print(
            f"[dEdx] median={dedx_p[50]:.2f} MeV/cm, p84={dedx_p[84]:.2f}, p95={dedx_p[95]:.2f}, "
            f"p99={dedx_p[99]:.2f}, p99.9={dedx_p[99.9]:.2f}"
        )
        hi_zoom = min(max(dedx_p[99], 5.0), 30.0)
        hi_tail = max(dedx_p[99.9], hi_zoom)

        fig, ax = plt.subplots()
        ax.hist(dedx_plot, bins=100, range=(0, hi_zoom), histtype="step", lw=2, label="trackLen>0.01 cm")
        ax.set_xlabel("dE/dx [MeV/cm]")
        ax.set_ylabel("Events")
        ax.set_xlim(0, hi_zoom)
        ax.set_title("dE/dx = Edep/trackLen (zoom)")
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
        ax.set_title("dE/dx tail (log-y)")
        ax.legend()
        fig.tight_layout()
        fig.savefig(plots_dir / "fig_dEdx_tail.pdf")
        plt.close(fig)

    # Pixelization
    valid_mask = (events["EdepCCD"] > 0) & (events["trackLenCCD"] > 0)
    valid_idxs = np.nonzero(valid_mask)[0]
    metrics_all: List[Dict[str, float]] = []
    metrics_quality: List[Dict[str, float]] = []
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
        plot_pixel_metrics(
            metrics_all,
            metrics_quality,
            plots_dir,
            args.quality_only,
            max_canvas_pix=args.max_canvas,
            fail_on_warning=args.fail_on_warning,
        )

    # Tables
    truncated_count = max(0, len(metrics_all) - len(metrics_quality))
    validation_summary = {
        "n_events": int(n_events),
        "n_events_ccd": int(valid_mask.sum()),
        "n_pixel_metrics_all": int(len(metrics_all)),
        "n_pixel_metrics_quality": int(len(metrics_quality)),
        "n_truncated_flagged": int(truncated_count),
        "truncated_fraction": float(truncated_count / len(metrics_all)) if metrics_all else 0.0,
        "cos_clamp_low": clamp_low,
        "cos_clamp_high": clamp_high,
    }
    validation_summary.update({f"energy_{k}": v for k, v in basic_stats(events["EevtPri"]).items()})
    validation_summary.update({f"EdepCCD_{k}": v for k, v in basic_stats(events["EdepCCD"]).items()})
    validation_summary.update({f"trackLenCCD_{k}": v for k, v in basic_stats(events["trackLenCCD"]).items()})
    validation_summary["impact_uniformity"] = impact_uniformity_stat(events["muonXImp"], events["muonYImp"])
    if dedx_all.size:
        validation_summary.update({f"dEdx_{k}": v for k, v in basic_stats(dedx_all).items()})
    if metrics_quality:
        validation_summary.update({f"sigma_trans_{k}": v for k, v in basic_stats(np.asarray([m["sigma_trans"] for m in metrics_quality])).items()})
        validation_summary.update({f"length_pix_geom_{k}": v for k, v in basic_stats(np.asarray([m["length_pix_geom"] for m in metrics_quality])).items()})
    write_csv(tables_dir / "validation_summary.csv", [validation_summary])

    if metrics_all:
        write_csv(tables_dir / "pixel_metrics_all.csv", metrics_all)
    if metrics_quality:
        write_csv(tables_dir / "pixel_metrics_quality.csv", metrics_quality)

    write_units_file(tables_dir / "units_and_conventions.txt", params)

    # Save config
    save_json(
        {
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
        },
        out_base / "run_config.json",
    )


if __name__ == "__main__":
    main()
