"""WSI reading and scaling utilities.

Wraps the physical image loader (``large_image``) so that the rest of the
package is independent of the medical-imaging I/O backend.
"""

from pathlib import Path
from typing import Any, Dict, Tuple, Union

import large_image
import numpy as np

# Micrometers-per-pixel ranges used to infer the objective power when the
# slide metadata does not report it (0.25 MPP ~ 40x, 0.50 ~ 20x, 1.0 ~ 10x).
MPP_TO_MAGNIFICATION: Tuple[Tuple[float, float, float], ...] = (
    (0.20, 0.30, 40.0),
    (0.40, 0.60, 20.0),
    (0.80, 1.20, 10.0),
)

DEFAULT_MAGNIFICATION = 20.0
DEFAULT_THUMBNAIL_SIZE = 2048


def open_wsi(image_path: Union[str, Path]) -> large_image.tilesource.FileTileSource:
    """Open a whole-slide image safely.

    Args:
        image_path: Path to the WSI file (e.g. ``.svs``, ``.tiff``).

    Returns:
        The ``large_image`` tile source.

    Raises:
        FileNotFoundError: If the file does not exist.
        RuntimeError: If ``large_image`` cannot parse the image format.
    """
    path_obj = Path(image_path)
    if not path_obj.exists():
        raise FileNotFoundError(f"WSI file not found at: {image_path}")

    try:
        return large_image.getTileSource(str(path_obj))
    except Exception as exc:
        raise RuntimeError(f"Error opening WSI with large_image at {image_path}: {exc}") from exc


def _infer_magnification_from_mpp(mpp_x: float) -> Union[float, None]:
    """Map a micrometers-per-pixel value to the closest standard objective power."""
    for low, high, magnification in MPP_TO_MAGNIFICATION:
        if low <= mpp_x <= high:
            return magnification
    return None


def _magnification_from_internal_metadata(
    tile_source: large_image.tilesource.FileTileSource,
) -> Union[float, None]:
    """Read the objective power from the OpenSlide properties, if exposed."""
    try:
        internal_meta = tile_source.getInternalMetadata()
        openslide_meta = internal_meta.get("openslide", {})
        obj_power = openslide_meta.get("openslide.objective-power", None)
        return float(obj_power) if obj_power is not None else None
    except Exception:
        return None


def get_wsi_metadata(tile_source: large_image.tilesource.FileTileSource) -> Dict[str, Any]:
    """Extract standardised metadata from an open WSI.

    Missing values are filled in by successive fallbacks: the objective power is
    inferred from the pixel spacing, then read from the OpenSlide properties, and
    finally defaults to 20x with a warning. The pixel spacing is in turn derived
    from the magnification when absent. Any value produced by a fallback is an
    assumption, so it is reported on stdout.

    Args:
        tile_source: Open WSI object.

    Returns:
        Dictionary with ``width`` and ``height`` (native pixels), ``magnification``
        (nominal objective power) and ``mpp_x`` / ``mpp_y`` (micrometers per pixel).
    """
    metadata = tile_source.getMetadata()

    width = metadata.get("sizeX", 0)
    height = metadata.get("sizeY", 0)
    magnification = metadata.get("magnification", None)

    # large_image reports spacing in millimetres; convert to micrometers
    raw_mm_x = metadata.get("mm_x", None)
    raw_mm_y = metadata.get("mm_y", None)
    mpp_x = raw_mm_x * 1000.0 if raw_mm_x is not None else None
    mpp_y = raw_mm_y * 1000.0 if raw_mm_y is not None else None

    if magnification is None and mpp_x is not None:
        magnification = _infer_magnification_from_mpp(mpp_x)

    if magnification is None:
        magnification = _magnification_from_internal_metadata(tile_source)

    if magnification is None:
        magnification = DEFAULT_MAGNIFICATION
        print(f"WARNING: Magnification missing in WSI metadata. Assuming default {DEFAULT_MAGNIFICATION}X.")

    if mpp_x is None:
        mpp_x = 0.50 if magnification == 20.0 else 0.25
        print(f"WARNING: mm_x missing in WSI metadata. Assuming {mpp_x} MPP based on {magnification}X.")
    if mpp_y is None:
        mpp_y = 0.50 if magnification == 20.0 else 0.25
        print(f"WARNING: mm_y missing in WSI metadata. Assuming {mpp_y} MPP based on {magnification}X.")

    return {
        "width": width,
        "height": height,
        "magnification": magnification,
        "mpp_x": mpp_x,
        "mpp_y": mpp_y,
    }


def drop_alpha_channel(image: np.ndarray) -> np.ndarray:
    """Return the image without its alpha channel (RGBA to RGB); other inputs pass through."""
    if image.ndim == 3 and image.shape[-1] == 4:
        return image[:, :, :3]
    return image


def get_wsi_thumbnail(
    tile_source: large_image.tilesource.FileTileSource,
    target_size: int = DEFAULT_THUMBNAIL_SIZE,
) -> Tuple[np.ndarray, float, float]:
    """Extract an RGB thumbnail and the scale factors that map native to thumbnail space.

    The thumbnail is constrained so that its largest side fits inside a
    ``target_size`` box, preserving the aspect ratio.

    Args:
        tile_source: Open WSI object.
        target_size: Maximum size of the largest side of the thumbnail.

    Returns:
        Tuple of (thumbnail in RGB, ``rx``, ``ry``), where ``rx = W_thumb / W_native``
        and ``ry = H_thumb / H_native``.
    """
    thumb_data, _ = tile_source.getThumbnail(
        width=target_size,
        height=target_size,
        format=large_image.constants.TILE_FORMAT_NUMPY,
    )
    thumb_rgb = drop_alpha_channel(thumb_data)

    native_meta = get_wsi_metadata(tile_source)
    native_w = native_meta["width"]
    native_h = native_meta["height"]
    thumb_h, thumb_w = thumb_rgb.shape[:2]

    rx = thumb_w / native_w if native_w > 0 else 0.0
    ry = thumb_h / native_h if native_h > 0 else 0.0

    return thumb_rgb, rx, ry
