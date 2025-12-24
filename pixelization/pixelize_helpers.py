import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Tuple, Union

import numpy as np

DEFAULT_PIXEL_SIZE_MICRONS = 15.0
DEFAULT_CCD_THICKNESS_MICRONS = 725.0
EV_PER_ELECTRON = 3.7


@dataclass
class CCDParams:
    pixel_size_microns: float = DEFAULT_PIXEL_SIZE_MICRONS
    thickness_microns: float = DEFAULT_CCD_THICKNESS_MICRONS
    # Active silicon area (default: 0.9 cm x 0.6 cm from CONNIE-like CCD)
    width_microns: float = 9000.0
    height_microns: float = 6000.0
    # Diffusion tuned to give sigma ~6-8 px at full depth for clear but moderate thickening
    diffusion_alpha: float = 1400.0  # micron^2 scale in sqrt(|alpha * ln(1 - beta z)|)
    diffusion_beta: float = 1.0 / (DEFAULT_CCD_THICKNESS_MICRONS + 1.0)


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def diffusion_sigma(z_microns: float, params: CCDParams) -> float:
    z = max(0.0, min(z_microns, params.thickness_microns))
    # simple monotonic diffusion curve similar to funciones_Sim_ab_initio.diffution_curve
    inside = max(1e-6, 1.0 - params.diffusion_beta * z)
    return math.sqrt(abs(params.diffusion_alpha * math.log(inside)))


def geant_to_electrons(edep_GeV: float) -> float:
    # GeV -> eV -> e-
    return max(0.0, edep_GeV * 1.0e9 / EV_PER_ELECTRON)


def build_track_image(
    edep_GeV: float,
    track_len_cm: float,
    theta: float,
    phi: float,
    entry_cm: Tuple[float, float, float],
    params: CCDParams,
    image_size: Union[int, Tuple[int, int], None] = None,
    center_origin: bool = False,
    margin_pix: int = 20,
) -> np.ndarray:
    """
    Render a simple 2D pixelized image of a muon track.

    - edep_GeV: total energy deposited in CCD (GeV)
    - track_len_cm: chord length through CCD (cm)
    - theta/phi: primary direction (rad)
    - entry_cm: entry position in cm (x,y,z)
    """
    pix = params.pixel_size_microns
    electrons = geant_to_electrons(edep_GeV)

    # Map real CCD coordinates (cm) into a downsampled canvas that spans the
    # whole active area. We keep the physical pixel size for diffusion and only
    # downscale when placing into the output image.
    width_pix = params.width_microns / pix
    height_pix = params.height_microns / pix
    if image_size is None:
        img_w = int(round(width_pix)) + 2 * margin_pix
        img_h = int(round(height_pix)) + 2 * margin_pix
    elif isinstance(image_size, Iterable):
        img_w, img_h = (int(image_size[0]), int(image_size[1]))
    else:
        img_w = img_h = int(image_size)

    img = np.zeros((img_h, img_w), dtype=np.float32)
    if track_len_cm <= 0.0 or edep_GeV <= 0.0:
        return img

    scale_x = width_pix / img_w
    scale_y = height_pix / img_h

    # Entry point in canvas coordinates (0,0 lower-left)
    entry_native = np.array(
        [
            entry_cm[0] * 1.0e4 / pix + width_pix / 2.0 + margin_pix,
            entry_cm[1] * 1.0e4 / pix + height_pix / 2.0 + margin_pix,
        ]
    )
    entry_canvas = entry_native / np.array([scale_x, scale_y])
    if center_origin:
        # Keep legacy behavior of recentering to avoid clipping when using small crops.
        cx = (img_w - 1) / 2.0
        cy = (img_h - 1) / 2.0
        entry_canvas = entry_canvas - np.array([img_w / 2.0, img_h / 2.0]) + np.array([cx, cy])

    # Displacement from entry to exit in native pixels, then downscaled.
    dx_native = track_len_cm * math.sin(theta) * math.cos(phi) * 1.0e4 / pix
    dy_native = track_len_cm * math.sin(theta) * math.sin(phi) * 1.0e4 / pix
    disp_canvas = np.array([dx_native / scale_x, dy_native / scale_y])
    start = entry_canvas
    end = entry_canvas + disp_canvas

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


def image_moments(img: np.ndarray, threshold: float = 0.0) -> Dict[str, float]:
    mask = img > threshold
    if not mask.any():
        return {"charge": 0.0, "npix": 0, "sigma_x": 0.0, "sigma_y": 0.0, "length_pix": 0.0}
    coords = np.argwhere(mask)
    charges = img[mask]
    total = charges.sum()
    mean_yx = (coords * charges[:, None]).sum(axis=0) / total
    dyx = coords - mean_yx
    cov = (charges[:, None, None] * np.einsum("ni,nj->nij", dyx, dyx)).sum(axis=0) / total
    sigma_x = math.sqrt(max(0.0, cov[1, 1]))
    sigma_y = math.sqrt(max(0.0, cov[0, 0]))
    # crude length proxy = max extent in mask
    y_extent = coords[:, 0].max() - coords[:, 0].min() + 1
    x_extent = coords[:, 1].max() - coords[:, 1].min() + 1
    length_pix = float(max(x_extent, y_extent))
    return {
        "charge": float(total),
        "npix": int(mask.sum()),
        "sigma_x": sigma_x,
        "sigma_y": sigma_y,
        "length_pix": length_pix,
    }


def save_json(obj: Dict, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2))
