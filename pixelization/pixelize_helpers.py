import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Tuple

import numpy as np

DEFAULT_PIXEL_SIZE_MICRONS = 15.0
DEFAULT_CCD_THICKNESS_MICRONS = 725.0
EV_PER_ELECTRON = 3.7


@dataclass
class CCDParams:
    pixel_size_microns: float = DEFAULT_PIXEL_SIZE_MICRONS
    thickness_microns: float = DEFAULT_CCD_THICKNESS_MICRONS
    diffusion_alpha: float = 12.0  # rough scale in micron^2
    diffusion_beta: float = 1.0 / (DEFAULT_CCD_THICKNESS_MICRONS + 1.0)


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def diffusion_sigma(z_microns: float, params: CCDParams) -> float:
    z = max(0.0, min(z_microns, params.thickness_microns))
    # simple monotonic diffusion curve similar to funciones_Sim_ab_initio.diffution_curve
    inside = max(1e-8, 1.0 - params.diffusion_beta * z)
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
    image_size: int = 64,
    center_origin: bool = True,
) -> np.ndarray:
    """
    Render a simple 2D pixelized image of a muon track.

    - edep_GeV: total energy deposited in CCD (GeV)
    - track_len_cm: chord length through CCD (cm)
    - theta/phi: primary direction (rad)
    - entry_cm: entry position in cm (x,y,z)
    """
    img = np.zeros((image_size, image_size), dtype=np.float32)
    if track_len_cm <= 0.0 or edep_GeV <= 0.0:
        return img

    pix = params.pixel_size_microns
    thickness_cm = params.thickness_microns * 1e-4  # micron -> cm
    electrons = geant_to_electrons(edep_GeV)

    # Projected XY displacement (cm) using chord through CCD thickness if available.
    dz = thickness_cm
    dx = track_len_cm * math.sin(theta) * math.cos(phi)
    dy = track_len_cm * math.sin(theta) * math.sin(phi)
    # start at center
    cx = cy = (image_size - 1) / 2.0 if center_origin else 0.0
    start = np.array([cx - dx / (2 * pix * 1e-4), cy - dy / (2 * pix * 1e-4)])
    end = np.array([cx + dx / (2 * pix * 1e-4), cy + dy / (2 * pix * 1e-4)])

    # Line parameterization
    num_samples = max(2, int(np.hypot(*(end - start)) * 2))
    xs = np.linspace(start[0], end[0], num_samples)
    ys = np.linspace(start[1], end[1], num_samples)

    charge_per_sample = electrons / num_samples
    # smear transversely with a Gaussian whose sigma grows with depth
    for i, (x, y) in enumerate(zip(xs, ys)):
        depth_microns = (i / max(1, num_samples - 1)) * params.thickness_microns
        sigma_pix = diffusion_sigma(depth_microns, params) / pix
        sigma_pix = max(0.5, min(5.0, sigma_pix))
        # influence region ~3 sigma
        xmin = int(max(0, math.floor(x - 3 * sigma_pix)))
        xmax = int(min(image_size - 1, math.ceil(x + 3 * sigma_pix)))
        ymin = int(max(0, math.floor(y - 3 * sigma_pix)))
        ymax = int(min(image_size - 1, math.ceil(y + 3 * sigma_pix)))
        if xmin > xmax or ymin > ymax:
            continue
        xs_grid = np.arange(xmin, xmax + 1)
        ys_grid = np.arange(ymin, ymax + 1)
        xv, yv = np.meshgrid(xs_grid, ys_grid, indexing="xy")
        dist2 = (xv - x) ** 2 + (yv - y) ** 2
        weights = np.exp(-0.5 * dist2 / (sigma_pix ** 2))
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
