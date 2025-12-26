#!/usr/bin/env python3
"""
Quick validation for Tier B cosmic muon sampling.

Generates basic histograms for cos(theta), energy, and source-plane (x0,y0),
and reports the CCD hit fraction and effective livetime if available.
"""

import argparse
import sys
from pathlib import Path
from typing import Dict

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

try:
    import uproot  # noqa: E402
except ImportError as exc:
    raise SystemExit("This script requires uproot (pip install uproot).") from exc

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from pixelization.pixelize_helpers import ensure_dir  # noqa: E402


def load_events(path: Path) -> Dict[str, np.ndarray]:
    tree = uproot.open(path)["B02Evts"]
    return {name: tree[name].array(library="np") for name in tree.keys()}


def load_run_info(path: Path) -> Dict[str, object]:
    f = uproot.open(path)
    if "B02RunInfo" not in f:
        return {}
    tree = f["B02RunInfo"]
    result: Dict[str, object] = {}
    for name in tree.keys():
        vals = tree[name].array(library="np")
        if vals.size == 0:
            continue
        v = vals[0]
        if isinstance(v, bytes):
            v = v.decode("utf-8")
        if hasattr(v, "item"):
            v = v.item()
        result[name] = v
    return result


def dedx_mev_per_cm(edep_GeV: np.ndarray, track_len_cm: np.ndarray) -> np.ndarray:
    mask = (track_len_cm > 0) & (edep_GeV > 0)
    return edep_GeV[mask] * 1000.0 / track_len_cm[mask]


def main():
    parser = argparse.ArgumentParser(description="Validate Tier B muon sampling.")
    parser.add_argument("--input", required=True, help="Path to B02ntuples.root")
    parser.add_argument("--outdir", default="paper_outputs/tierB_validation", help="Output directory for plots")
    parser.add_argument("--max-events", type=int, default=None, help="Optional cap on events to speed up plotting")
    args = parser.parse_args()

    input_path = Path(args.input)
    out_dir = ensure_dir(Path(args.outdir))
    events = load_events(input_path)
    run_info = load_run_info(input_path)
    ccd_thickness_cm = float(run_info.get("ccdThickness_cm", 0.0))
    warnings = []
    if args.max_events:
        for key, val in events.items():
            events[key] = val[: args.max_events]

    n_events = len(events["EevtPri"])
    edep = events.get("EdepCCD", np.array([]))
    tracklen = events.get("trackLenCCD", np.array([]))
    valid_mask = (edep > 0) & (tracklen > 0)
    hit_fraction = float(np.sum(valid_mask)) / float(n_events) if n_events else 0.0
    geom_mask = None
    geom_fraction = None
    if "geomIntersectsCCD" in events:
        geom_mask = np.asarray(events["geomIntersectsCCD"]) > 0.5
        geom_fraction = float(np.sum(geom_mask)) / float(n_events) if n_events else 0.0

    weights = None
    if "eventLivetime_s" in events:
        weights = np.asarray(events["eventLivetime_s"])
    elif "muonWeight_s" in events:
        weights = np.asarray(events["muonWeight_s"])
    livetime = float(np.sum(weights)) if weights is not None else 0.0
    hit_rate = float(np.sum(valid_mask)) / livetime if livetime > 0 else 0.0

    cos_down = None
    if "muonCosTheta" in events:
        cos_down = np.asarray(events["muonCosTheta"])
    else:
        cos_down = np.clip(-np.cos(events["thetaPri"]), 0.0, 1.0)

    energy = np.asarray(events["EevtPri"])
    x0 = events.get("muonX0")
    y0 = events.get("muonY0")

    print(f"Events: {n_events}, CCD hits: {int(np.sum(valid_mask))}, hit fraction={hit_fraction:.4f}")
    if livetime > 0:
        print(f"Effective livetime: {livetime:.2f} s, hit rate={hit_rate:.3f} Hz")

    fig, ax = plt.subplots()
    ax.hist(cos_down, bins=80, range=(0, 1), histtype="step", lw=2)
    ax.set_xlabel("cos(theta_down)")
    ax.set_ylabel("Events")
    ax.set_title("Cosine zenith (downward)")
    fig.tight_layout()
    fig.savefig(out_dir / "cos_theta_down.pdf")
    plt.close(fig)

    fig, ax = plt.subplots()
    ax.hist(energy, bins=120, histtype="step", lw=2, log=True)
    ax.set_xlabel("E_mu [GeV]")
    ax.set_ylabel("Events")
    ax.set_title("Muon energy spectrum")
    fig.tight_layout()
    fig.savefig(out_dir / "energy_spectrum.pdf")
    plt.close(fig)

    fig, ax = plt.subplots()
    h = ax.hist2d(cos_down, np.log10(np.clip(energy, 1e-3, None)), bins=60, cmap="viridis")
    plt.colorbar(h[3], ax=ax)
    ax.set_xlabel("cos(theta_down)")
    ax.set_ylabel("log10(E_mu/GeV)")
    ax.set_title("Energy vs cos(theta)")
    fig.tight_layout()
    fig.savefig(out_dir / "energy_vs_costheta.pdf")
    plt.close(fig)

    if x0 is not None and y0 is not None:
        fig, ax = plt.subplots()
        h = ax.hist2d(x0, y0, bins=80, cmap="viridis")
        plt.colorbar(h[3], ax=ax)
        ax.set_xlabel("x0 [cm]")
        ax.set_ylabel("y0 [cm]")
        ax.set_title("Source plane (x0,y0)")
        fig.tight_layout()
        fig.savefig(out_dir / "source_plane_xy0.pdf")
        plt.close(fig)

    if "muonXImp" in events and "muonYImp" in events:
        fig, ax = plt.subplots()
        h = ax.hist2d(events["muonXImp"], events["muonYImp"], bins=80, cmap="viridis")
        plt.colorbar(h[3], ax=ax)
        ax.set_xlabel("xImp [cm]")
        ax.set_ylabel("yImp [cm]")
        ax.set_title("Predicted impact at z=0")
        for x in (-0.75, 0.75):
            ax.axvline(x, color="red", linestyle="--", linewidth=1)
        for y in (-0.75, 0.75):
            ax.axhline(y, color="red", linestyle="--", linewidth=1)
        fig.tight_layout()
        fig.savefig(out_dir / "impact_plane_xy.pdf")
        plt.close(fig)

    # L*cos(theta) sanity check
    if ccd_thickness_cm > 0.0 and "trackLenCCD" in events:
        lcos = np.asarray(events["trackLenCCD"]) * np.abs(cos_down)
        lcos = lcos[lcos > 0]
        if lcos.size:
            lcos_median = float(np.median(lcos))
            ratio = lcos_median / ccd_thickness_cm if ccd_thickness_cm > 0 else 0.0
            fig, ax = plt.subplots()
            ax.hist(lcos, bins=80, histtype="step", lw=2)
            ax.axvline(ccd_thickness_cm, color="red", linestyle="--", label="CCD thickness")
            ax.set_xlabel("trackLen * |cos(theta)| [cm]")
            ax.set_ylabel("Events")
            ax.legend()
            ax.set_title("L*cos(theta) vs thickness")
            fig.tight_layout()
            fig.savefig(out_dir / "fig_Lcos_vs_thickness.pdf")
            plt.close(fig)
            if ratio < 0.5 or ratio > 1.5:
                warnings.append(f"Lcos median {lcos_median:.4f} cm vs thickness {ccd_thickness_cm:.4f} cm")
        else:
            warnings.append("No track lengths to evaluate Lcos.")

    # dE/dx check
    if "EdepCCD" in events and "trackLenCCD" in events:
        dedx = dedx_mev_per_cm(events["EdepCCD"], events["trackLenCCD"])
        if dedx.size:
            median_dedx = float(np.median(dedx))
            if median_dedx < 0.5 or median_dedx > 10.0:
                warnings.append(f"dE/dx median {median_dedx:.2f} MeV/cm outside 0.5-10 range.")

    summary = {
        "n_events": int(n_events),
        "n_hits_ccd": int(np.sum(valid_mask)),
        "hit_fraction": hit_fraction,
        "livetime_s": livetime,
        "hit_rate_hz": hit_rate,
    }
    if geom_fraction is not None:
        summary["geom_intersect_fraction"] = geom_fraction
    if warnings:
        print("Warnings:")
        for w in warnings:
            print(f"  - {w}")
        summary["warnings"] = "; ".join(warnings)
    (out_dir / "summary.txt").write_text(
        "\n".join(f"{k}: {v}" for k, v in summary.items()), encoding="utf-8"
    )


if __name__ == "__main__":
    main()
