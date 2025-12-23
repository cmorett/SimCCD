#!/usr/bin/env python3
"""
Generate paper-ready validation plots/tables and offline pixelized examples.

Usage:
  python analysis/make_paper_outputs.py --input build/vs2022/main/Release/B02ntuples.root --tag demo --examples 100 --dist-events 20000
"""

import argparse
import json
import math
import os
from pathlib import Path
from typing import Dict, List, Tuple

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

from pixelization.pixelize_helpers import (  # noqa: E402
    CCDParams,
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
        return {k: 0.0 for k in ["mean", "median", "rms", "p05", "p95"]}
    return {
        "mean": float(np.mean(vals)),
        "median": float(np.median(vals)),
        "rms": float(np.sqrt(np.mean(vals ** 2))),
        "p05": float(np.percentile(vals, 5)),
        "p95": float(np.percentile(vals, 95)),
    }


def make_hist(ax, data, bins, xlabel, ylabel="Events", title=None, logy=False, annotate=None):
    ax.hist(data, bins=bins, histtype="step", lw=2)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title)
    if logy:
        ax.set_yscale("log")
    if annotate:
        ax.text(0.95, 0.95, annotate, transform=ax.transAxes, ha="right", va="top")


def make_2d(ax, x, y, bins, xlabel, ylabel, title=None):
    h = ax.hist2d(x, y, bins=bins, cmap="viridis")
    plt.colorbar(h[3], ax=ax)
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


def render_examples(events: Dict[str, np.ndarray], idxs: np.ndarray, outdir: Path, params: CCDParams):
    imgs = []
    images_dir = ensure_dir(outdir / "images")
    for i, idx in enumerate(idxs):
        img = build_track_image(
            events["EdepCCD"][idx],
            events["trackLenCCD"][idx],
            events["thetaPri"][idx],
            events["phiPri"][idx],
            (events["xEntryCCD"][idx], events["yEntryCCD"][idx], events["zEntryCCD"][idx]),
            params,
        )
        imgs.append(img)
        plt.imsave(images_dir / f"event_{idx}.png", img, origin="lower", cmap="magma")
    # grid plot
    ncols = 5
    nrows = math.ceil(len(imgs) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 2.5, nrows * 2.5))
    for ax, img in zip(axes.flat, imgs):
        ax.imshow(img, origin="lower", cmap="magma")
        ax.axis("off")
    for ax in axes.flat[len(imgs) :]:
        ax.axis("off")
    fig.tight_layout()
    fig.savefig(outdir / "fig_pixelized_examples.pdf")
    plt.close(fig)


def pixel_metrics(events: Dict[str, np.ndarray], idxs: np.ndarray, params: CCDParams):
    metrics = []
    for idx in idxs:
        img = build_track_image(
            events["EdepCCD"][idx],
            events["trackLenCCD"][idx],
            events["thetaPri"][idx],
            events["phiPri"][idx],
            (events["xEntryCCD"][idx], events["yEntryCCD"][idx], events["zEntryCCD"][idx]),
            params,
        )
        metrics.append(image_moments(img, threshold=np.max(img) * 0.05 if np.max(img) > 0 else 0.0))
    return metrics


def plot_metrics(metrics: List[Dict[str, float]], outdir: Path):
    keys = {
        "npix": ("Cluster size [pix]", "fig_cluster_size.pdf"),
        "charge": ("Total charge [e-]", "fig_cluster_charge.pdf"),
        "sigma_x": ("Sigma X [pix]", "fig_sigma_x.pdf"),
        "sigma_y": ("Sigma Y [pix]", "fig_sigma_y.pdf"),
        "length_pix": ("Length [pix]", "fig_length_pix.pdf"),
    }
    for key, (xlabel, fname) in keys.items():
        vals = [m[key] for m in metrics if m[key] is not None]
        if not vals:
            continue
        fig, ax = plt.subplots()
        ax.hist(vals, bins=60, histtype="step", lw=2)
        ax.set_xlabel(xlabel)
        ax.set_ylabel("Events")
        fig.tight_layout()
        fig.savefig(outdir / fname)
        plt.close(fig)


def impact_uniformity_stat(x, y, bins=20):
    H, _, _ = np.histogram2d(x, y, bins=bins)
    expected = np.mean(H)
    if expected <= 0:
        return 0.0
    return float(np.max(np.abs(H - expected)) / expected)


def main():
    parser = argparse.ArgumentParser(description="Make paper-ready outputs for SimCCD.")
    parser.add_argument("--input", required=True, help="Path to B02ntuples ROOT file")
    parser.add_argument("--output", default="paper_outputs", help="Base output directory")
    parser.add_argument("--tag", default="run", help="Tag name for outputs")
    parser.add_argument("--examples", type=int, default=200, help="Number of pixelized examples to render")
    parser.add_argument("--dist-events", default="50000", help="Events to use for pixelization metrics ('all' or int)")
    parser.add_argument("--seed", type=int, default=12345, help="Numpy seed for reproducibility")
    args = parser.parse_args()

    np.random.seed(args.seed)
    input_path = Path(args.input)
    out_base = Path(args.output) / args.tag
    plots_dir = ensure_dir(out_base)
    tables_dir = ensure_dir(out_base / "tables")

    tree = load_tree(input_path)
    events = {name: arr(tree, name) for name in tree.keys()}
    n_events = len(events["EevtPri"])

    # Generator/source plots
    fig, ax = plt.subplots()
    make_hist(ax, events["EevtPri"], bins=120, xlabel="E_{#mu} [GeV]", logy=True, title="Primary muon energy")
    fig.tight_layout()
    fig.savefig(plots_dir / "fig_energy_spectrum.pdf")
    plt.close(fig)

    cos_down = np.clip(-np.cos(events["thetaPri"]), 0.0, 1.0)
    clamp_low = int(np.sum(-np.cos(events["thetaPri"]) < 0))
    clamp_high = int(np.sum(-np.cos(events["thetaPri"]) > 1))
    fig, ax = plt.subplots()
    make_hist(
        ax,
        cos_down,
        bins=60,
        xlabel="cos_{zenith}^{↓}",
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
    make_hist(ax, events["muonZ0"], bins=40, xlabel="z0 [cm]", title="Source height")
    fig.tight_layout()
    fig.savefig(plots_dir / "fig_z0.pdf")
    plt.close(fig)

    # CCD summaries
    fig, ax = plt.subplots()
    make_hist(ax, events["EdepCCD"], bins=100, xlabel="Edep CCD [GeV]", title="CCD deposited energy")
    fig.tight_layout()
    fig.savefig(plots_dir / "fig_edep_ccd.pdf")
    plt.close(fig)

    fig, ax = plt.subplots()
    make_hist(ax, events["trackLenCCD"], bins=100, xlabel="track length CCD [cm]", title="CCD track length")
    fig.tight_layout()
    fig.savefig(plots_dir / "fig_trackLen_ccd.pdf")
    plt.close(fig)

    fig, ax = plt.subplots()
    make_2d(
        ax,
        events["trackLenCCD"],
        events["EdepCCD"],
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
        xlabel="cos_{zenith}^{↓}",
        ylabel="track length CCD [cm]",
        title="Incidence vs CCD length",
    )
    fig.tight_layout()
    fig.savefig(plots_dir / "fig_costheta_vs_trackLen.pdf")
    plt.close(fig)

    # Pixelization
    params = CCDParams()
    valid_mask = (events["EdepCCD"] > 0) & (events["trackLenCCD"] > 0)
    valid_idxs = np.nonzero(valid_mask)[0]
    if valid_idxs.size == 0:
        print("WARNING: No events with CCD deposits; skipping pixelization outputs.")
        metrics = []
    else:
        # Examples
        n_examples = min(args.examples, valid_idxs.size)
        example_idxs = np.random.choice(valid_idxs, size=n_examples, replace=False)
        render_examples(events, example_idxs, plots_dir, params)

        # Metrics subset
        if args.dist_events == "all":
            metrics_idxs = valid_idxs
        else:
            n_dist = min(int(args.dist_events), valid_idxs.size)
            metrics_idxs = np.random.choice(valid_idxs, size=n_dist, replace=False)
        metrics = pixel_metrics(events, metrics_idxs, params)
        plot_metrics(metrics, plots_dir)

    # Tables
    validation_summary = {
        "n_events": int(n_events),
        "cos_clamp_low": clamp_low,
        "cos_clamp_high": clamp_high,
    }
    validation_summary.update({f"energy_{k}": v for k, v in basic_stats(events["EevtPri"]).items()})
    validation_summary.update({f"EdepCCD_{k}": v for k, v in basic_stats(events["EdepCCD"]).items()})
    validation_summary.update({f"trackLenCCD_{k}": v for k, v in basic_stats(events["trackLenCCD"]).items()})
    validation_summary["impact_uniformity"] = impact_uniformity_stat(events["muonXImp"], events["muonYImp"])
    write_csv(tables_dir / "validation_summary.csv", [validation_summary])

    if metrics:
        write_csv(tables_dir / "pixel_metrics.csv", metrics)

    units_text = """Units and conventions:
- thetaPri: Geant4 polar angle (rad); downward cos(zenith) = -cos(thetaPri).
- phiPri: Geant4 azimuthal angle (rad).
- EevtPri: primary kinetic energy (GeV).
- muonX0/Y0/Z0: primary vertex position (cm).
- muonXImp/YImp/ZImp: impact-plane sampling point (cm).
- EdepCCD: total deposited energy in CCD sensitive volume (GeV).
- trackLenCCD: chord length between first and last CCD step (cm).
- dirX/dirY/dirZ: primary momentum direction (unit vector).
- Pixelization: pixel size 15 microns, CCD thickness 725 microns, 3.7 eV per electron.
"""
    (tables_dir / "units_and_conventions.txt").write_text(units_text)

    # Save config
    save_json(
        {
            "input": str(input_path),
            "output": str(out_base),
            "examples": args.examples,
            "dist_events": args.dist_events,
            "seed": args.seed,
            "ccd_params": vars(params),
        },
        out_base / "run_config.json",
    )


if __name__ == "__main__":
    main()
