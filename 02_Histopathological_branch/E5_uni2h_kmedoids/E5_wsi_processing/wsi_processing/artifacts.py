"""Artifact filtering and patch extraction.

Two stages are provided:

* **Macro level** — :func:`filter_tissue_artifacts` removes pen markings from a
  raw tissue mask and smooths its boundaries with morphological operations.
* **Micro level** — :func:`extract_valid_patches` walks the slide on a regular
  grid at a target magnification, keeps the tiles that are sufficiently covered
  by tissue and pass the blur and contrast checks, and flags them on a
  low-resolution QC mask. :func:`artifact_filtering` is a thin wrapper that
  returns only that mask, without holding the patch pixels in memory.
"""

import os
from typing import Dict, Iterator, List, Sequence, Tuple

import large_image
import numpy as np
from skimage import color, exposure, filters, morphology
from tqdm import tqdm

from .utils import drop_alpha_channel, get_wsi_metadata

# RGB thresholds for the pen markers commonly used by pathologists on slides
PENS_RGB: Dict[str, List[Tuple[int, int, int]]] = {
    "red": [
        (120, 80, 90), (110, 20, 30), (185, 65, 105),
        (195, 85, 125), (220, 115, 145), (125, 40, 70),
        (200, 120, 150), (100, 50, 65), (85, 25, 45)
    ],
    "green": [
        (150, 160, 140), (70, 110, 110), (45, 115, 100),
        (30, 75, 60), (195, 220, 210), (225, 230, 225),
        (170, 210, 200), (20, 30, 20), (50, 60, 40),
        (30, 50, 35), (65, 70, 60), (100, 110, 105),
        (165, 180, 180), (140, 140, 150), (185, 195, 195)
    ],
    "blue": [
        (60, 120, 190), (120, 170, 200), (120, 170, 200),
        (175, 210, 230), (145, 210, 210), (37, 95, 160),
        (30, 65, 130), (130, 155, 180), (40, 35, 85),
        (30, 20, 65), (90, 90, 140), (60, 60, 120),
        (110, 110, 175)
    ],
    "black": [
        (70, 70, 70)
    ],
}

DEFAULT_PEN_COLORS: Tuple[str, ...] = ("red", "green", "blue", "black")


def _tqdm_position() -> int:
    """Progress-bar row, so parallel workers do not overwrite each other's bars."""
    return int(os.environ.get("TQDM_POSITION", "0"))


def marker_detection(image: np.ndarray, pen_color: str) -> np.ndarray:
    """Detect strokes of one pen marker colour on a slide thumbnail.

    Each colour is defined by a list of RGB thresholds; a pixel is flagged if it
    satisfies any of them. The comparison direction differs per colour: for red,
    the red channel must be high while green and blue are low, and so on.

    Args:
        image: Thumbnail in RGB format (H, W, 3).
        pen_color: One of ``'red'``, ``'green'``, ``'blue'``, ``'black'``.

    Returns:
        Boolean mask (H, W), True on the detected strokes.

    Raises:
        ValueError: If the colour is not supported.
    """
    if pen_color not in PENS_RGB:
        raise ValueError(f"Unsupported pen color for detection: {pen_color}")

    r, g, b = image[:, :, 0], image[:, :, 1], image[:, :, 2]
    thresholds = PENS_RGB[pen_color]
    mask = np.zeros_like(r, dtype=bool)

    if pen_color == "red":
        for t in thresholds:
            mask |= (r > t[0]) & (g < t[1]) & (b < t[2])
    elif pen_color == "green":
        for t in thresholds:
            mask |= (r < t[0]) & (g > t[1]) & (b > t[2])
    elif pen_color == "blue":
        for t in thresholds:
            mask |= (r < t[0]) & (g < t[1]) & (b > t[2])
    elif pen_color == "black":
        t = thresholds[0]
        mask = (r < t[0]) & (g < t[1]) & (b < t[2])

    return mask


def filter_tissue_artifacts(
    raw_tissue_mask: np.ndarray,
    thumb_rgb: np.ndarray,
    disk_size_opening: int = 1,
    disk_size_erosion: int = 1,
    filter_markers: bool = True,
    colors_to_filter: Sequence[str] = DEFAULT_PEN_COLORS,
) -> np.ndarray:
    """Clean a raw tissue mask by excluding pen markings and smoothing it.

    Args:
        raw_tissue_mask: Raw boolean tissue mask (H, W).
        thumb_rgb: Slide thumbnail in RGB format (H, W, 3).
        disk_size_opening: Radius of the structuring element used for binary
            opening, which removes small isolated objects.
        disk_size_erosion: Radius of the structuring element used for binary
            erosion, which trims boundary fragments.
        filter_markers: Whether to detect and exclude pen strokes.
        colors_to_filter: Pen colours to look for.

    Returns:
        The cleaned tissue mask (H, W).
    """
    if filter_markers:
        thumb_rgb = drop_alpha_channel(thumb_rgb)

        # Accumulate the strokes of every requested colour
        pen_mask = np.zeros(raw_tissue_mask.shape[:2], dtype=bool)
        for color_name in colors_to_filter:
            pen_mask |= marker_detection(thumb_rgb, color_name)

        clean_mask = raw_tissue_mask & (~pen_mask)
    else:
        clean_mask = raw_tissue_mask.copy()

    if disk_size_opening > 0:
        clean_mask = morphology.opening(clean_mask, morphology.disk(disk_size_opening))

    if disk_size_erosion > 0:
        clean_mask = morphology.erosion(clean_mask, morphology.disk(disk_size_erosion))

    return clean_mask


def _patch_geometry(
    tile_source: large_image.tilesource.FileTileSource,
    patch_width: int,
    patch_height: int,
    target_magnification: float,
    overlap: float,
) -> Tuple[int, int, int, int, int, int]:
    """Native-resolution patch size and strides for the requested magnification.

    Returns:
        ``(native_w, native_h, w_native, h_native, stride_x, stride_y)``.

    Raises:
        ValueError: If the native magnification cannot be determined.
    """
    meta = get_wsi_metadata(tile_source)
    native_w, native_h = meta["width"], meta["height"]
    native_mag = meta["magnification"]

    if not native_mag:
        raise ValueError("Could not determine WSI native magnification.")

    scale_factor = native_mag / target_magnification
    w_native = int(patch_width * scale_factor)
    h_native = int(patch_height * scale_factor)

    stride_x = int(w_native * (1.0 - overlap))
    stride_y = int(h_native * (1.0 - overlap))
    return native_w, native_h, w_native, h_native, stride_x, stride_y


def _iter_valid_patches(
    tile_source: large_image.tilesource.FileTileSource,
    raw_tissue_mask: np.ndarray,
    rx: float,
    ry: float,
    patch_width: int,
    patch_height: int,
    target_magnification: float,
    overlap: float,
    tissue_threshold: float,
    sharpness_threshold: float,
    contrast_threshold: float,
    show_progress: bool,
) -> Iterator[Tuple[int, int, np.ndarray, Tuple[int, int, int, int]]]:
    """Walk the slide grid and yield the tiles that pass every quality check.

    A tile is read from disk only after its tissue coverage has been verified on
    the low-resolution mask, which is what keeps the I/O cost manageable. Tiles
    are then rejected if they are low contrast or blurred (variance of the
    Laplacian below ``sharpness_threshold``).

    Yields:
        ``(x, y, region, (x_start, y_start, x_end, y_end))`` — the native
        coordinates of the tile, its RGB pixels, and its bounding box in
        low-resolution mask space.
    """
    native_w, native_h, w_native, h_native, stride_x, stride_y = _patch_geometry(
        tile_source, patch_width, patch_height, target_magnification, overlap
    )

    x_coords = list(range(0, native_w - w_native, stride_x))
    y_coords = list(range(0, native_h - h_native, stride_y))
    grid_coords = [(x, y) for y in y_coords for x in x_coords]

    iterator = grid_coords
    if show_progress:
        iterator = tqdm(grid_coords, desc="Fine QC Filtering", leave=False, position=_tqdm_position())

    for x, y in iterator:
        # Map the native bounding box into low-resolution mask space
        x_start, y_start = int(x * rx), int(y * ry)
        x_end, y_end = int((x + w_native) * rx), int((y + h_native) * ry)

        mask_patch = raw_tissue_mask[y_start:y_end, x_start:x_end]
        if mask_patch.size == 0:
            continue

        # Tissue coverage check on the mask, before any disk read
        tissue_ratio = np.sum(mask_patch) / mask_patch.size
        if tissue_ratio < tissue_threshold:
            continue

        try:
            region, _ = tile_source.getRegion(
                region=dict(left=x, top=y, width=w_native, height=h_native),
                scale=dict(magnification=target_magnification),
                format=large_image.constants.TILE_FORMAT_NUMPY,
            )
        except Exception:
            # Occasional read errors on damaged tiles: skip the tile
            continue

        region = drop_alpha_channel(region)
        if region.shape != (patch_height, patch_width, 3):
            continue

        if exposure.is_low_contrast(region, fraction_threshold=contrast_threshold):
            continue

        gray_patch = color.rgb2gray(region)
        if filters.laplace(gray_patch).var() < sharpness_threshold:
            continue

        yield x, y, region, (x_start, y_start, x_end, y_end)


def extract_valid_patches(
    tile_source: large_image.tilesource.FileTileSource,
    raw_tissue_mask: np.ndarray,
    rx: float,
    ry: float,
    slide_id: str,
    patch_width: int = 256,
    patch_height: int = 256,
    target_magnification: float = 20.0,
    overlap: float = 0.0,
    tissue_threshold: float = 0.8,
    sharpness_threshold: float = 0.0005,
    contrast_threshold: float = 0.05,
) -> Tuple[np.ndarray, np.ndarray, List[str], np.ndarray]:
    """Extract the valid patches of a slide, with their coordinates and QC mask.

    Every returned patch is held in memory, so a large slide can require several
    GB. Use :func:`artifact_filtering` when only the mask is needed.

    Args:
        tile_source: Open WSI source.
        raw_tissue_mask: Low-resolution binary tissue mask (H_mask, W_mask).
        rx: Scale factor on the X axis (``width_mask / width_native``).
        ry: Scale factor on the Y axis (``height_mask / height_native``).
        slide_id: Filename or identifier of the slide.
        patch_width: Patch width in pixels at the target magnification.
        patch_height: Patch height in pixels at the target magnification.
        target_magnification: Objective power at which patches are read.
        overlap: Overlap fraction between adjacent tiles, in [0.0, 1.0).
        tissue_threshold: Minimum tissue ratio required in the mask, in [0.0, 1.0].
        sharpness_threshold: Minimum variance of the Laplacian (blur rejection).
        contrast_threshold: Low-contrast tolerance passed to ``is_low_contrast``.

    Returns:
        Tuple of ``(coordinates, patches, slide_ids, quality_mask)``:
        coordinates of shape (P, 2) as ``[x, y]`` at native resolution (int32),
        patches of shape (P, patch_height, patch_width, 3) (uint8), the slide
        identifier repeated P times, and the refined QC mask.
    """
    quality_mask = np.zeros_like(raw_tissue_mask, dtype=bool)
    coords_list, patches_list = [], []

    for x, y, region, (x_start, y_start, x_end, y_end) in _iter_valid_patches(
        tile_source, raw_tissue_mask, rx, ry, patch_width, patch_height,
        target_magnification, overlap, tissue_threshold, sharpness_threshold,
        contrast_threshold, show_progress=True,
    ):
        coords_list.append([x, y])
        patches_list.append(region)
        quality_mask[y_start:y_end, x_start:x_end] = True

    if coords_list:
        coordinates = np.array(coords_list, dtype=np.int32)
        patches = np.array(patches_list, dtype=np.uint8)
    else:
        coordinates = np.empty((0, 2), dtype=np.int32)
        patches = np.empty((0, patch_height, patch_width, 3), dtype=np.uint8)

    slide_ids_list = [slide_id] * len(coords_list)
    return coordinates, patches, slide_ids_list, quality_mask


def artifact_filtering(
    tile_source: large_image.tilesource.FileTileSource,
    raw_tissue_mask: np.ndarray,
    rx: float,
    ry: float,
    patch_width: int = 256,
    patch_height: int = 256,
    target_magnification: float = 20.0,
    overlap: float = 0.0,
    tissue_threshold: float = 0.8,
    sharpness_threshold: float = 0.0005,
    contrast_threshold: float = 0.05,
) -> np.ndarray:
    """Return only the refined QC mask, without keeping the patch pixels.

    Same filtering as :func:`extract_valid_patches`; use this one when the
    patches themselves are not needed, as it does not accumulate them in memory.

    Args:
        See :func:`extract_valid_patches`.

    Returns:
        Refined QC mask, same shape as ``raw_tissue_mask``.
    """
    quality_mask = np.zeros_like(raw_tissue_mask, dtype=bool)

    for _, _, _, (x_start, y_start, x_end, y_end) in _iter_valid_patches(
        tile_source, raw_tissue_mask, rx, ry, patch_width, patch_height,
        target_magnification, overlap, tissue_threshold, sharpness_threshold,
        contrast_threshold, show_progress=False,
    ):
        quality_mask[y_start:y_end, x_start:x_end] = True

    return quality_mask
