#!/usr/bin/env python
"""Step 2 - Patch extraction and colorimetric descriptors for whole-slide images.

Run after ``wsi_step1_segmentation.py`` with the same ``--out-dir``. For every slide
the script tiles the level-0 image on a regular, non-overlapping grid, keeps the
tiles that are sufficiently covered by tissue (according to the step-1 mask) and
that pass a brightness / saturation check, and stores for each kept tile:

* a 1000x1000 px JPEG, named ``<k>.jpg`` with ``k`` = 0, 1, 2, ... (patch index);
* one row in ``patch_info.csv`` with its colorimetric descriptors.

Tile size depends on the objective power of the slide: 1000x1000 px at 40x and
500x500 px at 20x, so that every tile covers the same tissue area as a 500x500 px
patch at 20x. Descriptors are computed on a 500x500 px version of the tile; the
JPEG is resized to 1000x1000 px so that 16 patches build a 4000x4000 px mosaic.

``patch_info.csv`` layout (consumed by the mosaic-building notebook)::

    ID, cell index, <14 colorimetric descriptors>, tissue_cov

The notebook treats every column after the first two as a descriptor, so
``tissue_cov`` (last column) is part of the descriptor vector used to compare patches.

Memory note: the tissue mask is resized to the level-0 dimensions of the slide
(one byte per pixel), which can require several GB for large 40x slides.

Usage
-----
    python wsi_step2_patch_extraction.py --wsi-dir /path/to/slides --out-dir /path/to/processed
"""

import argparse
import logging
import sys
from pathlib import Path

import cv2
import numpy as np
import openslide
import pandas as pd
from PIL import Image

# ----------------------------------------------------------------------------
# Parameters
# ----------------------------------------------------------------------------
PATCH_FINAL_SIZE = 500     # size (px) at which descriptors are computed (20x-equivalent)
PATCH_MOSAIC_SIZE = 1000   # size (px) of the saved JPEG (4 x 4 patches = 4000 x 4000 mosaic)
MIN_TISSUE_COV = 0.70      # minimum fraction of the tile covered by tissue
MAX_BRIGHTNESS = 230       # tiles with mean CIELAB L above this value are discarded (background)
MIN_SATURATION = 10        # tiles with mean HSV saturation below this value are discarded
EPSILON = 1e-6             # avoids log(0) in optical-density computations
JPEG_QUALITY = 95

FEATURE_COLS = ["R_mean", "G_mean", "B_mean", "R_std", "G_std", "B_std",
                "H_mean", "S_mean", "V_mean", "OD_R", "OD_G", "OD_B",
                "brightness", "contrast"]

log = logging.getLogger("wsi_step2")


def infer_magnification(slide):
    """Objective power of the slide; falls back to a size-based guess if the metadata are missing."""
    try:
        return int(float(slide.properties[openslide.PROPERTY_NAME_OBJECTIVE_POWER]))
    except (KeyError, ValueError):
        w, h = slide.level_dimensions[0]
        return 20 if max(w, h) < 30000 else 40


def load_mask_resized(mask_path, target_shape):
    """Load the step-1 mask and resize it (nearest neighbour) to ``target_shape`` = (height, width)."""
    mask = cv2.imread(str(mask_path), cv2.IMREAD_GRAYSCALE)
    return cv2.resize(mask, (target_shape[1], target_shape[0]), interpolation=cv2.INTER_NEAREST)


def cell_index(patch_rgb):
    """Cellularity proxy: mean RGB optical density of the patch (higher = more stain)."""
    rgb_norm = patch_rgb.astype(np.float32) / 255.0
    od = -np.log(rgb_norm + EPSILON)
    return float(od.mean())


def colorimetric_features(patch_rgb):
    """14 colorimetric descriptors, in the order given by ``FEATURE_COLS``."""
    r, g, b = patch_rgb[:, :, 0], patch_rgb[:, :, 1], patch_rgb[:, :, 2]
    feats = [r.mean(), g.mean(), b.mean(), r.std(), g.std(), b.std()]
    hsv = cv2.cvtColor(patch_rgb, cv2.COLOR_RGB2HSV)
    feats += [hsv[:, :, 0].mean(), hsv[:, :, 1].mean(), hsv[:, :, 2].mean()]
    rgb_n = patch_rgb.astype(np.float32) / 255.0
    feats += [(-np.log(rgb_n[:, :, c] + EPSILON)).mean() for c in range(3)]
    lab = cv2.cvtColor(patch_rgb, cv2.COLOR_RGB2LAB)
    feats += [float(lab[:, :, 0].mean()), float(lab[:, :, 0].std())]
    return feats


def extract_patches(svs_path, out_dir):
    """Extract patches from one slide. Returns True if at least one patch was kept."""
    slide_id = svs_path.stem
    slide_dir = out_dir / slide_id
    mask_path = slide_dir / f"{slide_id}_tissue_mask.png"
    csv_path = slide_dir / "patch_info.csv"

    if not mask_path.exists():
        log.error("  Tissue mask not found: %s (run step 1 first)", slide_id)
        return False
    if csv_path.exists():           # regenerate from scratch
        csv_path.unlink()

    log.info("Extracting: %s", slide_id)
    slide = openslide.OpenSlide(str(svs_path))
    try:
        magnification = infer_magnification(slide)
        w0, h0 = slide.level_dimensions[0]

        # 40x: 1000x1000 px tiles; 20x: 500x500 px tiles (same tissue area, see module docstring)
        extract_size = 1000 if magnification == 40 else 500
        stride = extract_size
        log.info("  %dx | tile size %dx%d", magnification, extract_size, extract_size)

        mask_l0 = load_mask_resized(mask_path, target_shape=(h0, w0))

        records, kept, total = [], 0, 0
        for y in range(0, h0 - extract_size + 1, stride):
            for x in range(0, w0 - extract_size + 1, stride):
                total += 1
                mask_patch = mask_l0[y:y + extract_size, x:x + extract_size]
                tissue_cov = (mask_patch > 0).sum() / (extract_size ** 2)
                if tissue_cov < MIN_TISSUE_COV:
                    continue

                region = slide.read_region((x, y), 0, (extract_size, extract_size))
                patch_raw = np.array(region.convert("RGB"))
                patch_500 = cv2.resize(patch_raw, (PATCH_FINAL_SIZE, PATCH_FINAL_SIZE),
                                       interpolation=cv2.INTER_AREA)

                # Quality filter (background / unstained tiles)
                hsv_p = cv2.cvtColor(patch_500, cv2.COLOR_RGB2HSV)
                lab_p = cv2.cvtColor(patch_500, cv2.COLOR_RGB2LAB)
                if lab_p[:, :, 0].mean() > MAX_BRIGHTNESS or hsv_p[:, :, 1].mean() < MIN_SATURATION:
                    continue

                patch_1000 = cv2.resize(patch_raw, (PATCH_MOSAIC_SIZE, PATCH_MOSAIC_SIZE),
                                        interpolation=cv2.INTER_AREA)
                Image.fromarray(patch_1000).save(str(slide_dir / f"{kept}.jpg"), quality=JPEG_QUALITY)

                record = {"ID": kept, "cell index": cell_index(patch_500),
                          "tissue_cov": round(tissue_cov, 4)}
                record.update(dict(zip(FEATURE_COLS, colorimetric_features(patch_500))))
                records.append(record)
                kept += 1

        log.info("  Candidates: %d | Kept: %d", total, kept)
        if kept == 0:
            log.warning("  No valid patches: %s", slide_id)
            return False

        columns = ["ID", "cell index"] + FEATURE_COLS + ["tissue_cov"]
        pd.DataFrame(records)[columns].to_csv(str(csv_path), index=False)
        log.info("  OK: %d patches saved to %s/", kept, slide_dir)
        return True
    except Exception as exc:
        log.error("  Error in %s: %s", slide_id, exc, exc_info=True)
        return False
    finally:
        slide.close()


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--wsi-dir", type=Path, required=True, help="directory containing the slides")
    parser.add_argument("--out-dir", type=Path, required=True,
                        help="directory produced by step 1 (masks); patches are written next to them")
    parser.add_argument("--pattern", default="*.svs", help="glob pattern of the slide files (default: %(default)s)")
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[logging.StreamHandler(sys.stdout),
                  logging.FileHandler(args.out_dir / "patch_extraction.log")],
    )

    slides = sorted(args.wsi_dir.glob(args.pattern))
    log.info("Slides to process: %d", len(slides))

    ok, failed = [], []
    for i, svs in enumerate(slides, 1):
        log.info("[%d/%d]", i, len(slides))
        (ok if extract_patches(svs, args.out_dir) else failed).append(svs.stem)

    log.info("Summary: %d OK / %d failed", len(ok), len(failed))
    total_jpgs = sum(len(list((args.out_dir / s).glob("*.jpg"))) for s in ok)
    log.info("Total JPEG patches generated: %d", total_jpgs)
    if failed:
        log.warning("Failed: %s", failed)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
