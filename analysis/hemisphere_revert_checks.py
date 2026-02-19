#!/usr/bin/env python3
"""
Validation checks for the reverted hemisphere + Smith-Duller muon generator.
Produces plots and a JSON report with quantitative pass/fail checks.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Dict, Iterable, Tuple

import matplotlib.pyplot as plt
import numpy as np
import uproot


CCD_SIDE_CM = 5.461
EXPECTED_CCD_THICKNESS_CM = 0.07
SMITH_LOGE_MIN = -1.0
SMITH_LOGE_MAX = 5.0

EVENT_COLUMNS = [
    "muonCosTheta",
    "phiPri",
    "thetaPri",
    "EevtPri",
    "muonEnergySampledGeV",
    "muonX0",
    "muonY0",
    "muonZ0",
    "muonXImp",
    "muonYImp",
    "muonZImp",
    "dirX",
    "dirY",
    "dirZ",
    "trackLenCCD",
    "nStepsCCD",
    "EdepCCD",
    "muonModeCode",
    "isTargeted",
    "cfg_thetaMax_deg",
    "muonWeight_s",
    "eventLivetime_s",
    "prov_ccdThickness_cm",
]


def _as_float_array(arrays: Dict[str, np.ndarray], key: str) -> np.ndarray:
    if key not in arrays:
        return np.array([], dtype=float)
    return np.asarray(arrays[key], dtype=float)


def _as_int_array(arrays: Dict[str, np.ndarray], key: str) -> np.ndarray:
    if key not in arrays:
        return np.array([], dtype=int)
    return np.asarray(arrays[key], dtype=int)


def _finite(x: np.ndarray) -> np.ndarray:
    return x[np.isfinite(x)]


def load_events(path: Path) -> Dict[str, np.ndarray]:
    with uproot.open(path) as root:
        if "B02Evts" not in root:
            raise RuntimeError(f"Tree B02Evts not found in {path}")
        tree = root["B02Evts"]
        available = set(tree.keys())
        wanted = [c for c in EVENT_COLUMNS if c in available]
        arrays = tree.arrays(wanted, library="np")
    return dict(arrays)


def efficiency(mask: np.ndarray, total: int) -> Tuple[float, int, int]:
    hits = int(mask.sum())
    return (hits / float(total) if total else 0.0), hits, total


def ks_uniform_stat(u: np.ndarray) -> float:
    if u.size == 0:
        return float("nan")
    u_sorted = np.sort(np.clip(u, 0.0, 1.0))
    n = u_sorted.size
    grid = (np.arange(1, n + 1, dtype=float)) / n
    d_plus = np.max(grid - u_sorted)
    d_minus = np.max(u_sorted - (np.arange(0, n, dtype=float) / n))
    return float(max(d_plus, d_minus))


def ks_cos2_sin(mu: np.ndarray) -> float:
    mu = _finite(mu)
    if mu.size == 0:
        return float("nan")
    mu_sorted = np.sort(np.clip(mu, 0.0, 1.0))
    n = mu_sorted.size
    emp = (np.arange(1, n + 1, dtype=float)) / n
    expected = mu_sorted**3  # CDF for p(mu)=3*mu^2, mu in [0,1]
    return float(np.max(np.abs(emp - expected)))


def uniform_binned_stats(values: np.ndarray, low: float, high: float, bins: int) -> Dict[str, float]:
    values = _finite(values)
    if values.size == 0:
        return {"chi2_ndf": float("nan"), "max_rel_dev": float("nan")}
    counts, _ = np.histogram(values, bins=bins, range=(low, high))
    expected = values.size / bins
    if expected <= 0:
        return {"chi2_ndf": float("nan"), "max_rel_dev": float("nan")}
    chi2 = np.sum((counts - expected) ** 2 / expected)
    rel = np.abs(counts - expected) / expected
    return {
        "chi2_ndf": float(chi2 / max(1, bins - 1)),
        "max_rel_dev": float(np.max(rel)),
    }


def smith_duller_weights(loge_mev: np.ndarray, cos_theta: float) -> np.ndarray:
    cos_theta = float(np.clip(cos_theta, 1.0e-4, 1.0))
    e_mev = np.power(10.0, loge_mev)

    # Same constants as generator implementation.
    Au = 2.0e9
    gu = 2.645
    ru = 0.76
    au = 2.5
    y0u = 1000.0
    bmu = 0.80
    cu = 299792458.0e2
    mmu = 105.7 / (cu**2)
    t0mu = 2.2e-6
    r0u = 0.00129
    Bmu = bmu * mmu * y0u * cu / (t0mu * r0u)
    lpu = 120.0
    bu = 0.771
    mpu = 139.6 / (cu**2)
    t0pu = 2.6e-8
    jpu = mpu * y0u * cu / (t0pu * r0u)

    epu = (e_mev + au * y0u * (1.0 / cos_theta - 0.100)) / ru
    epu = np.clip(epu, 1.0e-12, None)
    base = 0.100 * cos_theta * (1.0 - (au * (y0u / cos_theta - 100.0) / (ru * epu)))
    base = np.clip(base, 1.0e-30, None)
    expo = Bmu / ((ru * epu + 100.0 * au) * cos_theta)
    pmu = np.power(base, expo)
    denom = np.clip(epu * cos_theta + bu * jpu, 1.0e-30, None)
    weights = Au * np.power(epu, -gu) * pmu * lpu * bu * jpu / denom
    weights = np.where(np.isfinite(weights) & (weights > 0.0), weights, 0.0)
    return weights


def smith_pit(loge_obs: np.ndarray, cos_obs: np.ndarray, n_cos_bins: int = 30) -> np.ndarray:
    valid = np.isfinite(loge_obs) & np.isfinite(cos_obs)
    loge_obs = loge_obs[valid]
    cos_obs = np.clip(cos_obs[valid], 1.0e-4, 1.0)
    if loge_obs.size == 0:
        return np.array([], dtype=float)

    grid = np.linspace(SMITH_LOGE_MIN, SMITH_LOGE_MAX, 1800)
    quantiles = np.linspace(0.0, 1.0, n_cos_bins + 1)
    edges = np.quantile(cos_obs, quantiles)
    edges = np.unique(edges)
    if edges.size < 2:
        edges = np.array([np.min(cos_obs), np.max(cos_obs) + 1e-6], dtype=float)

    idx = np.clip(np.digitize(cos_obs, edges) - 1, 0, edges.size - 2)
    pit = np.zeros_like(loge_obs, dtype=float)

    for bin_id in np.unique(idx):
        mask = idx == bin_id
        cmean = float(np.mean(cos_obs[mask]))
        w = smith_duller_weights(grid, cmean)
        cdf = np.cumsum(w)
        if cdf[-1] <= 0.0:
            pit[mask] = 0.5
            continue
        cdf /= cdf[-1]
        pit[mask] = np.interp(loge_obs[mask], grid, cdf, left=0.0, right=1.0)
    return pit


def tangent_plane_residuals(arrays: Dict[str, np.ndarray]) -> np.ndarray:
    x0 = _as_float_array(arrays, "muonX0")
    y0 = _as_float_array(arrays, "muonY0")
    z0 = _as_float_array(arrays, "muonZ0")
    x_imp = _as_float_array(arrays, "muonXImp")
    y_imp = _as_float_array(arrays, "muonYImp")
    z_imp = _as_float_array(arrays, "muonZImp")
    dx = _as_float_array(arrays, "dirX")
    dy = _as_float_array(arrays, "dirY")
    dz = _as_float_array(arrays, "dirZ")

    n = min(len(x0), len(y0), len(z0), len(x_imp), len(y_imp), len(z_imp), len(dx), len(dy), len(dz))
    if n == 0:
        return np.array([], dtype=float)
    x0, y0, z0 = x0[:n], y0[:n], z0[:n]
    x_imp, y_imp, z_imp = x_imp[:n], y_imp[:n], z_imp[:n]
    dx, dy, dz = dx[:n], dy[:n], dz[:n]

    valid = np.isfinite(x0) & np.isfinite(y0) & np.isfinite(z0)
    valid &= np.isfinite(x_imp) & np.isfinite(y_imp) & np.isfinite(z_imp)
    valid &= np.isfinite(dx) & np.isfinite(dy) & np.isfinite(dz)
    valid &= np.abs(dz) > 1.0e-12
    if not np.any(valid):
        return np.array([], dtype=float)

    t = (z_imp[valid] - z0[valid]) / dz[valid]
    x_calc = x0[valid] + t * dx[valid]
    y_calc = y0[valid] + t * dy[valid]
    res = np.sqrt((x_calc - x_imp[valid]) ** 2 + (y_calc - y_imp[valid]) ** 2)
    return res[np.isfinite(res)]


def plot_impact(x: np.ndarray, y: np.ndarray, side_cm: float, out: Path, title: str) -> None:
    half = 0.5 * side_cm
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(x, y, s=6, alpha=0.35, edgecolor="none")
    ax.add_patch(
        plt.Rectangle(
            (-half, -half),
            side_cm,
            side_cm,
            fill=False,
            lw=2.0,
            color="tab:red",
            label="CCD active area",
        )
    )
    ax.set_xlabel("x_imp [cm]")
    ax.set_ylabel("y_imp [cm]")
    ax.set_title(title)
    ax.set_xlim(-half * 1.1, half * 1.1)
    ax.set_ylim(-half * 1.1, half * 1.1)
    ax.set_aspect("equal", adjustable="box")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out, dpi=220)
    plt.close(fig)


def plot_angles(mu: np.ndarray, phi: np.ndarray, out: Path, label: str) -> None:
    mu = _finite(mu)
    phi = _finite(phi)
    centers = np.linspace(0.01, 0.99, 60)
    expected = 3.0 * centers**2

    fig, axes = plt.subplots(1, 3, figsize=(14, 4))
    if mu.size:
        hist, edges = np.histogram(mu, bins=np.linspace(0.0, 1.0, 61), density=True)
        x = 0.5 * (edges[:-1] + edges[1:])
        axes[0].step(x, hist, where="mid", color="tab:blue", label="sample")
    axes[0].plot(centers, expected, color="tab:red", label="3 mu^2")
    axes[0].set_xlabel("mu = cos(theta)")
    axes[0].set_ylabel("PDF")
    axes[0].set_title(f"{label}: cos(theta)")
    axes[0].legend()

    if mu.size:
        theta_deg = np.degrees(np.arccos(np.clip(mu, 0.0, 1.0)))
        axes[1].hist(theta_deg, bins=45, color="tab:purple", alpha=0.75)
    axes[1].set_xlabel("theta [deg]")
    axes[1].set_ylabel("counts")
    axes[1].set_title(f"{label}: theta")

    if phi.size:
        phi_wrapped = np.mod(phi, 2.0 * math.pi)
        axes[2].hist(phi_wrapped, bins=36, range=(0.0, 2.0 * math.pi), color="tab:green", alpha=0.75)
    axes[2].set_xlabel("phi [rad]")
    axes[2].set_ylabel("counts")
    axes[2].set_title(f"{label}: phi")

    fig.tight_layout()
    fig.savefig(out, dpi=220)
    plt.close(fig)


def plot_energy(loge: np.ndarray, out: Path, label: str) -> None:
    loge = _finite(loge)
    fig, ax = plt.subplots(figsize=(6, 4))
    if loge.size:
        ax.hist(loge, bins=np.linspace(SMITH_LOGE_MIN, SMITH_LOGE_MAX, 60), color="tab:orange", alpha=0.8)
    ax.set_xlabel("log10(E [MeV])")
    ax.set_ylabel("counts")
    ax.set_title(f"{label}: sampled energy")
    fig.tight_layout()
    fig.savefig(out, dpi=220)
    plt.close(fig)


def plot_tracklen(track_len: np.ndarray, out: Path, label: str) -> None:
    vals = _finite(track_len)
    vals = vals[vals > 0.0]
    fig, ax = plt.subplots(figsize=(6, 4))
    if vals.size:
        ax.hist(vals, bins=50, color="tab:cyan", alpha=0.8)
    ax.set_xlabel("Track length in CCD [cm]")
    ax.set_ylabel("counts")
    ax.set_title(f"{label}: CCD track length")
    fig.tight_layout()
    fig.savefig(out, dpi=220)
    plt.close(fig)


def summarize_dataset(label: str, arrays: Dict[str, np.ndarray], out_dir: Path) -> Dict[str, float]:
    out_dir.mkdir(parents=True, exist_ok=True)

    mu = np.clip(_as_float_array(arrays, "muonCosTheta"), 0.0, 1.0)
    phi = _as_float_array(arrays, "phiPri")
    energy_gev = _as_float_array(arrays, "muonEnergySampledGeV")
    if energy_gev.size == 0:
        energy_gev = _as_float_array(arrays, "EevtPri")
    loge = np.log10(np.clip(energy_gev * 1.0e3, 1.0e-12, None))

    x_imp = _as_float_array(arrays, "muonXImp")
    y_imp = _as_float_array(arrays, "muonYImp")
    track_len = _as_float_array(arrays, "trackLenCCD")
    edep_gev = _as_float_array(arrays, "EdepCCD")

    muon_mode = _as_int_array(arrays, "muonModeCode")
    targeted = _as_int_array(arrays, "isTargeted")
    theta_max_cfg = _as_float_array(arrays, "cfg_thetaMax_deg")
    ccd_thickness = _as_float_array(arrays, "prov_ccdThickness_cm")
    mu_w = _as_float_array(arrays, "muonWeight_s")
    livetime = _as_float_array(arrays, "eventLivetime_s")

    total = int(len(track_len)) if len(track_len) else int(len(mu))
    hits_mask = track_len > 0.0
    eff, hits, _ = efficiency(hits_mask, total)

    pit = smith_pit(loge, mu)
    phi_wrapped = np.mod(phi, 2.0 * math.pi)
    phi_stats = uniform_binned_stats(phi_wrapped, 0.0, 2.0 * math.pi, bins=36)
    x_stats = uniform_binned_stats(x_imp, -0.5 * CCD_SIDE_CM, 0.5 * CCD_SIDE_CM, bins=24)
    y_stats = uniform_binned_stats(y_imp, -0.5 * CCD_SIDE_CM, 0.5 * CCD_SIDE_CM, bins=24)

    tangent_res = tangent_plane_residuals(arrays)
    lcos = track_len[hits_mask] * mu[hits_mask] if mu.size == track_len.size else np.array([], dtype=float)
    d_edx = np.array([], dtype=float)
    if edep_gev.size == track_len.size and np.any(hits_mask):
        d_edx = (edep_gev[hits_mask] * 1.0e3) / np.clip(track_len[hits_mask], 1.0e-12, None)

    result: Dict[str, float] = {
        "n_events": float(total),
        "n_hits": float(hits),
        "efficiency": float(eff),
        "cos_theta_mean": float(np.mean(mu)) if mu.size else float("nan"),
        "cos_theta_std": float(np.std(mu)) if mu.size else float("nan"),
        "cos_theta_ks_vs_3mu2": ks_cos2_sin(mu),
        "phi_chi2_ndf": phi_stats["chi2_ndf"],
        "phi_max_rel_dev": phi_stats["max_rel_dev"],
        "ximp_chi2_ndf": x_stats["chi2_ndf"],
        "yimp_chi2_ndf": y_stats["chi2_ndf"],
        "energy_mean_GeV": float(np.mean(energy_gev)) if energy_gev.size else float("nan"),
        "energy_std_GeV": float(np.std(energy_gev)) if energy_gev.size else float("nan"),
        "energy_median_GeV": float(np.median(energy_gev)) if energy_gev.size else float("nan"),
        "smith_pit_ks_uniform": ks_uniform_stat(pit),
        "smith_pit_mean": float(np.mean(pit)) if pit.size else float("nan"),
        "smith_pit_std": float(np.std(pit)) if pit.size else float("nan"),
        "tangent_residual_rms_cm": float(np.sqrt(np.mean(tangent_res**2))) if tangent_res.size else float("nan"),
        "tangent_residual_p95_cm": float(np.percentile(tangent_res, 95)) if tangent_res.size else float("nan"),
        "trackLen_mean_cm_hits": float(np.mean(track_len[hits_mask])) if np.any(hits_mask) else float("nan"),
        "trackLen_median_cm_hits": float(np.median(track_len[hits_mask])) if np.any(hits_mask) else float("nan"),
        "lcos_median_cm_hits": float(np.median(lcos)) if lcos.size else float("nan"),
        "lcos_p16_cm_hits": float(np.percentile(lcos, 16)) if lcos.size else float("nan"),
        "lcos_p84_cm_hits": float(np.percentile(lcos, 84)) if lcos.size else float("nan"),
        "dEdx_median_MeV_per_cm_hits": float(np.median(d_edx)) if d_edx.size else float("nan"),
        "dEdx_p16_MeV_per_cm_hits": float(np.percentile(d_edx, 16)) if d_edx.size else float("nan"),
        "dEdx_p84_MeV_per_cm_hits": float(np.percentile(d_edx, 84)) if d_edx.size else float("nan"),
        "mode_code_median": float(np.median(muon_mode)) if muon_mode.size else float("nan"),
        "targeted_fraction": float(np.mean(targeted > 0)) if targeted.size else float("nan"),
        "cfg_thetaMax_deg_median": float(np.median(theta_max_cfg)) if theta_max_cfg.size else float("nan"),
        "ccd_thickness_median_cm": float(np.median(ccd_thickness)) if ccd_thickness.size else float("nan"),
        "muon_weight_mean_s": float(np.mean(mu_w)) if mu_w.size else float("nan"),
        "event_livetime_mean_s": float(np.mean(livetime)) if livetime.size else float("nan"),
        "pass_mode_forced_footprint": float(np.all(muon_mode == 0)) if muon_mode.size else 0.0,
        "pass_targeted_impact": float(np.mean(targeted > 0) > 0.99) if targeted.size else 0.0,
        "pass_cos_theta_shape": float(ks_cos2_sin(mu) < 0.04) if mu.size else 0.0,
        "pass_phi_uniformity": float(phi_stats["chi2_ndf"] < 2.5) if np.isfinite(phi_stats["chi2_ndf"]) else 0.0,
        "pass_smith_pit": float(ks_uniform_stat(pit) < 0.05) if pit.size else 0.0,
        "pass_not_fixed_energy": float(np.std(energy_gev) > 0.05) if energy_gev.size else 0.0,
        "pass_tangent_projection": float((np.sqrt(np.mean(tangent_res**2)) < 1.0e-9) if tangent_res.size else 0.0),
        "pass_ccd_thickness_0p07cm": float(abs(np.median(ccd_thickness) - EXPECTED_CCD_THICKNESS_CM) < 1.0e-6)
        if ccd_thickness.size
        else 0.0,
    }

    plot_impact(x_imp, y_imp, side_cm=CCD_SIDE_CM, out=out_dir / f"{label}_impact.png", title=f"{label}: impact plane")
    plot_angles(mu, phi, out=out_dir / f"{label}_angles.png", label=label)
    plot_energy(loge, out=out_dir / f"{label}_energy.png", label=label)
    plot_tracklen(track_len, out=out_dir / f"{label}_tracklen.png", label=label)

    return result


def geometry_consistency(datasets: Dict[str, Dict[str, float]]) -> Dict[str, float]:
    thickness_vals = []
    lcos_vals = []
    for stats in datasets.values():
        t = stats.get("ccd_thickness_median_cm", float("nan"))
        if np.isfinite(t):
            thickness_vals.append(float(t))
        l = stats.get("lcos_median_cm_hits", float("nan"))
        if np.isfinite(l):
            lcos_vals.append(float(l))

    out: Dict[str, float] = {}
    if thickness_vals:
        out["ccd_thickness_min_cm"] = float(np.min(thickness_vals))
        out["ccd_thickness_max_cm"] = float(np.max(thickness_vals))
        out["ccd_thickness_spread_cm"] = float(np.max(thickness_vals) - np.min(thickness_vals))
        out["pass_thickness_equal_across_representations"] = float((np.max(thickness_vals) - np.min(thickness_vals)) < 1.0e-6)
        out["pass_thickness_expected_0p07cm"] = float(np.max(np.abs(np.asarray(thickness_vals) - EXPECTED_CCD_THICKNESS_CM)) < 1.0e-6)
    else:
        out["pass_thickness_equal_across_representations"] = 0.0
        out["pass_thickness_expected_0p07cm"] = 0.0

    if lcos_vals:
        out["lcos_median_min_cm"] = float(np.min(lcos_vals))
        out["lcos_median_max_cm"] = float(np.max(lcos_vals))
    return out


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Hemisphere/Smith-Duller validation checks.")
    parser.add_argument("--tag", default="hemisphere_revert", help="Run tag used under outputs/<tag>/...")
    parser.add_argument("--run-dir", type=Path, help="Base output directory (defaults to outputs/<tag>).")
    parser.add_argument("--root-none", type=Path, help="Merged ROOT for geometry none.")
    parser.add_argument("--root-cad", type=Path, help="Merged ROOT for tessellated CAD.")
    parser.add_argument("--root-default", type=Path, help="Optional ROOT from default-settings probe run.")
    parser.add_argument("--root-primitive", type=Path, help="Optional ROOT from primitive-geometry probe run.")
    parser.add_argument("--out-dir", type=Path, help="Output directory for plots and report.")
    return parser.parse_args()


def _dataset_entries(args: argparse.Namespace) -> Iterable[Tuple[str, Path]]:
    base = args.run_dir if args.run_dir else Path("outputs") / args.tag
    none_root = args.root_none or (base / "none" / "merged.root")
    cad_root = args.root_cad or (base / "cad" / "merged.root")
    yield "none", none_root
    yield "tessellated", cad_root

    default_root = args.root_default or (base / "default_probe.root")
    primitive_root = args.root_primitive or (base / "primitive_probe.root")
    if default_root.exists():
        yield "default_probe", default_root
    if primitive_root.exists():
        yield "primitive_probe", primitive_root


def main() -> int:
    args = parse_args()
    base = args.run_dir if args.run_dir else Path("outputs") / args.tag
    out_dir = args.out_dir or (base / "sanity")
    out_dir.mkdir(parents=True, exist_ok=True)

    summaries: Dict[str, Dict[str, float]] = {}
    for label, path in _dataset_entries(args):
        if not path.exists():
            if label in ("none", "tessellated"):
                raise FileNotFoundError(f"Required ROOT file not found: {path}")
            continue
        arrays = load_events(path)
        summaries[label] = summarize_dataset(label, arrays, out_dir / label)

    global_checks = geometry_consistency(summaries)
    report = {"datasets": summaries, "global": global_checks}
    report_path = out_dir / "sanity_summary.json"
    report_path.write_text(json.dumps(report, indent=2), encoding="utf-8")

    print(f"[done] wrote validation report: {report_path}")
    for label, stats in summaries.items():
        eff = stats.get("efficiency", float("nan"))
        hits = int(stats.get("n_hits", 0))
        total = int(stats.get("n_events", 0))
        smith_ok = int(stats.get("pass_smith_pit", 0.0))
        ang_ok = int(stats.get("pass_cos_theta_shape", 0.0))
        print(
            f"{label}: efficiency={eff:.4f} ({hits}/{total}), "
            f"smith_ok={smith_ok}, angular_ok={ang_ok}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
