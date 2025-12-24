import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional, Tuple, Union

import numpy as np

DEFAULT_PIXEL_SIZE_MICRONS = 15.0
DEFAULT_CCD_THICKNESS_MICRONS = 725.0
EV_PER_ELECTRON = 3.7


@dataclass
class CCDParams:
    pixel_size_microns: float = DEFAULT_PIXEL_SIZE_MICRONS
    thickness_microns: float = DEFAULT_CCD_THICKNESS_MICRONS
    ev_per_electron: float = EV_PER_ELECTRON
    # Active silicon area (default: 0.9 cm x 0.6 cm from CONNIE-like CCD)
    width_microns: float = 9000.0
    height_microns: float = 6000.0
    # Diffusion tuned to give sigma ~6-8 px at full depth for clear but moderate thickening
    diffusion_alpha: float = 1400.0  # micron^2 scale in sqrt(|alpha * ln(1 - beta z)|)
    diffusion_beta: float = 1.0 / (DEFAULT_CCD_THICKNESS_MICRONS + 1.0)
    # Adaptive canvas guardrail
    max_canvas_pix: int = 512


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def diffusion_sigma(z_microns: float, params: CCDParams) -> float:
    z = max(0.0, min(z_microns, params.thickness_microns))
    # simple monotonic diffusion curve similar to funciones_Sim_ab_initio.diffution_curve
    inside = max(1e-6, 1.0 - params.diffusion_beta * z)
    return math.sqrt(abs(params.diffusion_alpha * math.log(inside)))


def geant_to_electrons(edep_GeV: float, params: CCDParams) -> float:
    # GeV -> eV -> e-
    return max(0.0, edep_GeV * 1.0e9 / params.ev_per_electron)


def _canvas_dimensions(
    track_len_cm: float,
    params: CCDParams,
    canvas_mode: str,
    canvas_size: Union[int, Tuple[int, int], None],
    margin_pix: int,
) -> Tuple[int, int, float, float, float]:
    """
    Returns (img_w, img_h, scale_x, scale_y, length_pix_geom).
    scale_x/y convert native CCD pixel units into canvas pixels.
    """
    pix = params.pixel_size_microns
    length_pix_geom = max(0.0, track_len_cm * 1.0e4 / pix)
    if canvas_mode == "adaptive":
        # Square canvas sized to expected track length plus margin.
        min_dim = max(32, 2 * margin_pix + 8)
        target = int(math.ceil(length_pix_geom + 2 * margin_pix))
        dim = int(min(params.max_canvas_pix, max(min_dim, target)))
        return dim, dim, 1.0, 1.0, length_pix_geom

    width_pix = params.width_microns / pix
    height_pix = params.height_microns / pix
    if canvas_size is None:
        img_w = int(round(width_pix)) + 2 * margin_pix
        img_h = int(round(height_pix)) + 2 * margin_pix
    elif isinstance(canvas_size, (tuple, list)):
        img_w, img_h = int(canvas_size[0]), int(canvas_size[1])
    else:
        img_w = img_h = int(canvas_size)

    scale_x = width_pix / img_w if img_w > 0 else 1.0
    scale_y = height_pix / img_h if img_h > 0 else 1.0
    return img_w, img_h, scale_x, scale_y, length_pix_geom


def build_track_image(
    edep_GeV: float,
    track_len_cm: float,
    theta: float,
    phi: float,
    entry_cm: Tuple[float, float, float],
    params: CCDParams,
    canvas_mode: str = "adaptive",
    canvas_size: Union[int, Tuple[int, int], None] = None,
    margin_pix: int = 20,
    center_origin: bool = False,
) -> np.ndarray:
    """
    Render a simple 2D pixelized image of a muon track.

    - edep_GeV: total energy deposited in CCD (GeV)
    - track_len_cm: chord length through CCD (cm)
    - theta/phi: primary direction (rad)
    - entry_cm: entry position in cm (x,y,z)
    - canvas_mode: "adaptive" sizes canvas to track length; "fixed" maps full CCD area.
    """
    img_w, img_h, scale_x, scale_y, _ = _canvas_dimensions(
        track_len_cm, params, canvas_mode, canvas_size, margin_pix
    )
    img = np.zeros((img_h, img_w), dtype=np.float32)
    if track_len_cm <= 0.0 or edep_GeV <= 0.0:
        return img

    pix = params.pixel_size_microns
    electrons = geant_to_electrons(edep_GeV, params)

    dx_native = track_len_cm * math.sin(theta) * math.cos(phi) * 1.0e4 / pix
    dy_native = track_len_cm * math.sin(theta) * math.sin(phi) * 1.0e4 / pix

    if canvas_mode == "adaptive":
        cx = (img_w - 1) / 2.0
        cy = (img_h - 1) / 2.0
        disp_canvas = np.array([dx_native / scale_x, dy_native / scale_y])
        start = np.array([cx, cy]) - 0.5 * disp_canvas
    else:
        width_pix = params.width_microns / pix
        height_pix = params.height_microns / pix
        entry_native = np.array(
            [
                entry_cm[0] * 1.0e4 / pix + width_pix / 2.0 + margin_pix,
                entry_cm[1] * 1.0e4 / pix + height_pix / 2.0 + margin_pix,
            ]
        )
        entry_canvas = entry_native / np.array([scale_x, scale_y])
        if center_origin:
            cx = (img_w - 1) / 2.0
            cy = (img_h - 1) / 2.0
            entry_canvas = entry_canvas - np.array([img_w / 2.0, img_h / 2.0]) + np.array([cx, cy])
        disp_canvas = np.array([dx_native / scale_x, dy_native / scale_y])
        start = entry_canvas

    end = start + disp_canvas

    track_len_microns = track_len_cm * 1.0e4
    num_samples = max(5, int(math.ceil(track_len_microns / (0.5 * pix))))
    xs = np.linspace(start[0], end[0], num_samples)
    ys = np.linspace(start[1], end[1], num_samples)

    charge_per_sample = electrons / num_samples
    # Depth mapping based on the actual z-span of the track (clamped to the CCD).
    z0 = entry_cm[2] * 1.0e4  # cm -> microns
    z_path = track_len_cm * 1.0e4 * math.cos(theta)
    z1 = z0 + z_path
    z0_c = max(0.0, min(params.thickness_microns, z0))
    z1_c = max(0.0, min(params.thickness_microns, z1))
    z_span = z1_c - z0_c

    for i, (x, y) in enumerate(zip(xs, ys)):
        t = i / max(1, num_samples - 1)
        if abs(z_span) > 1e-6:
            depth_microns = z0_c + t * z_span
        else:
            depth_microns = t * params.thickness_microns
        depth_sigma = diffusion_sigma(depth_microns, params) / pix
        # Add mild path-length broadening so thickness grows along the chord even at shallow angles.
        path_pix = track_len_microns / pix
        sigma_native = math.sqrt(max(0.0, depth_sigma ** 2 + (0.05 * t * path_pix) ** 2))
        sigma_native = min(max(0.25, sigma_native), 6.0)
        sigma_x = sigma_native / scale_x
        sigma_y = sigma_native / scale_y
        radius_x = 2.5 * sigma_x
        radius_y = 2.5 * sigma_y
        xmin = int(max(0, math.floor(x - radius_x)))
        xmax = int(min(img_w - 1, math.ceil(x + radius_x)))
        ymin = int(max(0, math.floor(y - radius_y)))
        ymax = int(min(img_h - 1, math.ceil(y + radius_y)))
        if xmin > xmax or ymin > ymax:
            continue
        xs_grid = np.arange(xmin, xmax + 1)
        ys_grid = np.arange(ymin, ymax + 1)
        xv, yv = np.meshgrid(xs_grid, ys_grid, indexing="xy")
        dist2 = ((xv - x) / sigma_x) ** 2 + ((yv - y) / sigma_y) ** 2
        weights = np.exp(-0.5 * dist2)
        weights_sum = weights.sum()
        if weights_sum <= 0:
            continue
        img[ymin : ymax + 1, xmin : xmax + 1] += charge_per_sample * weights / weights_sum
    return img


def _weighted_percentile(values: np.ndarray, weights: np.ndarray, percentile: float) -> float:
    sorter = np.argsort(values)
    values_sorted = values[sorter]
    weights_sorted = weights[sorter]
    cumulative = np.cumsum(weights_sorted)
    if cumulative[-1] <= 0:
        return float(values_sorted[0])
    cutoff = percentile / 100.0 * cumulative[-1]
    idx = np.searchsorted(cumulative, cutoff)
    idx = min(max(0, idx), len(values_sorted) - 1)
    return float(values_sorted[idx])


def image_moments(
    img: np.ndarray,
    threshold: float = 0.0,
    pixel_size_microns: float = DEFAULT_PIXEL_SIZE_MICRONS,
    track_len_cm: Optional[float] = None,
    edge_margin: int = 1,
    p_span: Tuple[float, float] = (1.0, 99.0),
) -> Dict[str, float]:
    """
    Extract per-image metrics, including PCA-based, and an edge-touch flag.
    """
    mask = img > threshold
    h, w = img.shape
    if not mask.any():
        return {
            "charge": 0.0,
            "npix": 0,
            "sigma_x": 0.0,
            "sigma_y": 0.0,
            "sigma_long": 0.0,
            "sigma_trans": 0.0,
            "elongation": 0.0,
            "phi_pca": 0.0,
            "length_pix_img": 0.0,
            "length_pix_geom": float(track_len_cm * 1.0e4 / pixel_size_microns) if track_len_cm else 0.0,
            "length_bbox_pix": 0.0,
            "is_truncated": False,
        }

    coords = np.argwhere(mask)  # (y, x)
    charges = img[mask]
    total = float(charges.sum())
    mean_yx = (coords * charges[:, None]).sum(axis=0) / total
    dyx = coords - mean_yx
    cov = (charges[:, None, None] * np.einsum("ni,nj->nij", dyx, dyx)).sum(axis=0) / total
    sigma_x = math.sqrt(max(0.0, cov[1, 1]))
    sigma_y = math.sqrt(max(0.0, cov[0, 0]))

    eigvals, eigvecs = np.linalg.eigh(cov)
    order = np.argsort(eigvals)[::-1]
    eigvals = eigvals[order]
    eigvecs = eigvecs[:, order]
    sigma_long = math.sqrt(max(0.0, eigvals[0]))
    sigma_trans = math.sqrt(max(0.0, eigvals[1]))
    elongation = sigma_long / max(sigma_trans, 1e-6)
    principal_vec = eigvecs[:, 0]  # y, x
    phi_pca = math.atan2(principal_vec[0], principal_vec[1])

    projections = dyx @ principal_vec
    p_lo, p_hi = p_span
    lo = _weighted_percentile(projections, charges, p_lo)
    hi = _weighted_percentile(projections, charges, p_hi)
    length_pix_img = max(0.0, hi - lo)

    y_extent = coords[:, 0].max() - coords[:, 0].min() + 1
    x_extent = coords[:, 1].max() - coords[:, 1].min() + 1
    edge_touch = bool(
        (coords[:, 0].min() <= edge_margin)
        or (coords[:, 0].max() >= h - 1 - edge_margin)
        or (coords[:, 1].min() <= edge_margin)
        or (coords[:, 1].max() >= w - 1 - edge_margin)
    )

    return {
        "charge": total,
        "npix": int(mask.sum()),
        "sigma_x": sigma_x,
        "sigma_y": sigma_y,
        "sigma_long": sigma_long,
        "sigma_trans": sigma_trans,
        "elongation": elongation,
        "phi_pca": phi_pca,
        "length_pix_img": length_pix_img,
        "length_pix_geom": float(track_len_cm * 1.0e4 / pixel_size_microns) if track_len_cm else 0.0,
        "length_bbox_pix": float(max(x_extent, y_extent)),
        "is_truncated": edge_touch,
    }


def save_json(obj: Dict, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2))
