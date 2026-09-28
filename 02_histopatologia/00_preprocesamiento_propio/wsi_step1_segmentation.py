#!/usr/bin/env python
"""Step 1 - Tissue segmentation of whole-slide images (WSI).

For every slide in ``--wsi-dir`` the script:

1. reads a low-magnification thumbnail (default 4x) with OpenSlide;
2. builds a texture map from the local RGB variance, enhances it with CLAHE and
   smooths it with a Gaussian filter;
3. thresholds the map with the triangle method, fills holes and removes small
   connected components;
4. removes pen marks (blue, green, black ink) and out-of-focus regions.

Outputs, written to ``<out-dir>/<slide_id>/``:

* ``<slide_id>_tissue_mask.png``    binary tissue mask (0 / 255);
* ``<slide_id>_thumbnail_4x.png``   RGB thumbnail;
* ``<slide_id>_overlay_qc.png``     thumbnail with the mask overlaid, for visual QC.

The slide identifier is the file name without extension, so slides must be named
after the identifier that the following steps expect. Step 2
(``wsi_step2_patch_extraction.py``) must be run with the same ``--out-dir``.

Usage
-----
    python wsi_step1_segmentation.py --wsi-dir /path/to/slides --out-dir /path/to/processed
"""

import argparse
import logging
import sys
from pathlib import Path

import cv2
import numpy as np
import openslide
from PIL import Image

# ----------------------------------------------------------------------------
# Parameters
# ----------------------------------------------------------------------------
TARGET_MAG = 4            # magnification of the thumbnail used for segmentation
VARIANCE_WIN = 5          # window (px) of the local-variance filter
CLAHE_CLIP = 2.0          # CLAHE clip limit
CLAHE_GRID = (8, 8)       # CLAHE tile grid
GAUSS_SIGMA = 2           # Gaussian smoothing (px)
MIN_TISSUE_AREA = 500     # minimum connected-component area (px) kept in the mask
MAX_DIM = 4096            # thumbnails larger than this are downscaled for segmentation
BLUR_BLOCK = 64           # window (px) of the local Laplacian-variance blur detector
BLUR_THRESHOLD = 50       # local Laplacian variance below this value = out of focus

log = logging.getLogger("wsi_step1")


def get_level_for_magnification(slide, target_mag):
    """Return the pyramid level whose downsample factor is closest to ``objective / target_mag``."""
    try:
        objective = float(slide.properties[openslide.PROPERTY_NAME_OBJECTIVE_POWER])
    except (KeyError, ValueError):
        log.warning("Objective power not available; assuming 40x")
        objective = 40.0
    downsample_target = objective / target_mag
    diffs = [abs(ds - downsample_target) for ds in slide.level_downsamples]
    best_level = int(np.argmin(diffs))
    log.debug("  objective=%sx | level=%d | downsample=%.2f",
              objective, best_level, slide.level_downsamples[best_level])
    return best_level


def rgb_variance_mask(rgb, win=VARIANCE_WIN):
    """Local RGB variance (mean over channels), rescaled to 0-255 (uint8)."""
    rgb_f = rgb.astype(np.float32)
    var_maps = []
    for c in range(3):
        ch = rgb_f[:, :, c]
        mean_sq = cv2.blur(ch ** 2, (win, win))
        sq_mean = cv2.blur(ch, (win, win)) ** 2
        var_maps.append(np.maximum(mean_sq - sq_mean, 0))
    variance = np.mean(var_maps, axis=0)
    vmin, vmax = variance.min(), variance.max()
    if vmax - vmin < 1e-6:
        return np.zeros(variance.shape, dtype=np.uint8)
    return ((variance - vmin) / (vmax - vmin) * 255).astype(np.uint8)


def apply_clahe(gray):
    """Contrast-limited adaptive histogram equalisation."""
    clahe = cv2.createCLAHE(clipLimit=CLAHE_CLIP, tileGridSize=CLAHE_GRID)
    return clahe.apply(gray)


def triangle_threshold(gray):
    """Binarise with the triangle method."""
    _, mask = cv2.threshold(gray, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_TRIANGLE)
    return mask


def clean_mask_morphology(mask, min_area=MIN_TISSUE_AREA):
    """Fill holes (flood fill from the corner) and drop components smaller than ``min_area``."""
    h, w = mask.shape
    filled = mask.copy()
    flood = np.zeros((h + 2, w + 2), dtype=np.uint8)
    cv2.floodFill(filled, flood, (0, 0), 255)
    filled_inv = cv2.bitwise_not(filled)
    mask_filled = mask | filled_inv
    num_labels, labels, stats, _ = cv2.connectedComponentsWithStats(mask_filled, connectivity=8)
    clean = np.zeros_like(mask_filled)
    for lbl in range(1, num_labels):
        if stats[lbl, cv2.CC_STAT_AREA] >= min_area:
            clean[labels == lbl] = 255
    return clean


def remove_artifacts(mask, rgb):
    """Remove pen marks (blue / green / black) and blurred regions from the tissue mask."""
    hsv = cv2.cvtColor(rgb, cv2.COLOR_RGB2HSV)
    hue, sat, val = hsv[:, :, 0], hsv[:, :, 1], hsv[:, :, 2]
    pen_blue = ((hue >= 100) & (hue <= 130) & (sat > 50)).astype(np.uint8) * 255
    pen_green = ((hue >= 40) & (hue <= 80) & (sat > 50)).astype(np.uint8) * 255
    pen_black = (val < 30).astype(np.uint8) * 255
    pen_mask = cv2.bitwise_or(pen_blue, cv2.bitwise_or(pen_green, pen_black))
    kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (15, 15))
    pen_mask = cv2.dilate(pen_mask, kernel, iterations=1)

    # Vectorised blur detection: local variance of the Laplacian, E[x^2] - E[x]^2
    gray = cv2.cvtColor(rgb, cv2.COLOR_RGB2GRAY).astype(np.float32)
    lap = cv2.Laplacian(gray, cv2.CV_32F)
    mean_lap2 = cv2.blur(lap ** 2, (BLUR_BLOCK, BLUR_BLOCK))
    mean_lap = cv2.blur(lap, (BLUR_BLOCK, BLUR_BLOCK))
    local_var = mean_lap2 - mean_lap ** 2
    blur_mask = (local_var < BLUR_THRESHOLD).astype(np.uint8) * 255

    artifact_mask = cv2.bitwise_or(pen_mask, blur_mask)
    return cv2.bitwise_and(mask, cv2.bitwise_not(artifact_mask))


def segment_wsi(svs_path, out_dir):
    """Segment one slide and write mask, thumbnail and QC overlay. Returns True on success."""
    slide_id = svs_path.stem
    log.info("Processing: %s", slide_id)
    try:
        slide = openslide.OpenSlide(str(svs_path))
    except Exception as exc:  # OpenSlide raises several unrelated exception types
        log.error("Could not open %s: %s", svs_path.name, exc)
        return False

    try:
        level = get_level_for_magnification(slide, TARGET_MAG)
        dims = slide.level_dimensions[level]
        rgb = np.array(slide.read_region((0, 0), level, dims).convert("RGB"))
        log.info("  Thumbnail %dx: %dx%d px", TARGET_MAG, dims[0], dims[1])

        # Cap the working size for very large slides
        h_orig, w_orig = rgb.shape[:2]
        scale = min(MAX_DIM / max(h_orig, w_orig), 1.0)
        if scale < 1.0:
            new_w, new_h = int(w_orig * scale), int(h_orig * scale)
            rgb = cv2.resize(rgb, (new_w, new_h), interpolation=cv2.INTER_AREA)
            log.info("  Resized to %dx%d (scale=%.2f)", new_w, new_h, scale)

        var_map = rgb_variance_mask(rgb)
        enhanced = apply_clahe(var_map)
        blurred = cv2.GaussianBlur(enhanced, (0, 0), GAUSS_SIGMA)
        mask_raw = triangle_threshold(blurred)
        mask_morph = clean_mask_morphology(mask_raw)
        mask_final = remove_artifacts(mask_morph, rgb)

        # Bring the mask (and the RGB image) back to the full thumbnail size
        if scale < 1.0:
            mask_final = cv2.resize(mask_final, (w_orig, h_orig), interpolation=cv2.INTER_NEAREST)
            rgb = np.array(slide.read_region((0, 0), level, dims).convert("RGB"))

        tissue_pct = 100 * mask_final.sum() / (255 * mask_final.size)
        log.info("  Tissue detected: %.1f%%", tissue_pct)
        if tissue_pct < 1.0:
            log.warning("  Less than 1%% tissue detected in %s", slide_id)

        slide_dir = out_dir / slide_id
        slide_dir.mkdir(exist_ok=True)
        cv2.imwrite(str(slide_dir / f"{slide_id}_tissue_mask.png"), mask_final)
        Image.fromarray(rgb).save(str(slide_dir / f"{slide_id}_thumbnail_4x.png"))

        overlay = rgb.copy()
        overlay[mask_final == 255] = (
            overlay[mask_final == 255] * 0.6 + np.array([0, 200, 0]) * 0.4
        ).astype(np.uint8)
        Image.fromarray(overlay).save(str(slide_dir / f"{slide_id}_overlay_qc.png"))
        log.info("  Saved to %s/", slide_dir)
        return True
    except Exception as exc:
        log.error("Error in %s: %s", slide_id, exc, exc_info=True)
        return False
    finally:
        slide.close()


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--wsi-dir", type=Path, required=True, help="directory containing the slides")
    parser.add_argument("--out-dir", type=Path, required=True, help="output directory (one sub-directory per slide)")
    parser.add_argument("--pattern", default="*.svs", help="glob pattern of the slide files (default: %(default)s)")
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[logging.StreamHandler(sys.stdout),
                  logging.FileHandler(args.out_dir / "segmentation.log")],
    )

    slides = sorted(args.wsi_dir.glob(args.pattern))
    log.info("Slides found: %d", len(slides))

    ok, failed = [], []
    for i, svs in enumerate(slides, 1):
        log.info("[%d/%d]", i, len(slides))
        (ok if segment_wsi(svs, args.out_dir) else failed).append(svs.stem)

    log.info("Summary: %d OK / %d failed", len(ok), len(failed))
    if failed:
        log.warning("Failed: %s", failed)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
