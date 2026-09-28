"""Tissue segmentation for H&E whole-slide images.

Isolates histological tissue from the bright glass background in low-resolution
slide images.
"""

import numpy as np
from skimage import color, filters

from .utils import drop_alpha_channel

THRESHOLD_METHODS = {
    "triangle": filters.threshold_triangle,
    "otsu": filters.threshold_otsu,
}


def get_raw_tissue_mask(
    thumb_rgb: np.ndarray,
    sigma: float = 1.0,
    threshold_method: str = "triangle",
) -> np.ndarray:
    """Generate a raw binary mask of tissue from an RGB thumbnail.

    Assumes a bright (white) background and darker tissue regions, as produced by
    H&E staining. The image is converted to grayscale, smoothed and inverted, so
    that tissue becomes the bright foreground the thresholding operates on.

    Args:
        thumb_rgb: WSI thumbnail in RGB format (H, W, 3).
        sigma: Standard deviation of the Gaussian smoothing, which suppresses
            high-frequency noise and granularity.
        threshold_method: Adaptive thresholding method, ``"triangle"`` or ``"otsu"``.

    Returns:
        Boolean mask (H, W), True where tissue is preliminarily detected.

    Raises:
        ValueError: If the thresholding method is not supported.
    """
    method_lower = threshold_method.lower()
    if method_lower not in THRESHOLD_METHODS:
        raise ValueError(f"Unsupported thresholding method: {threshold_method}")

    thumb_rgb = drop_alpha_channel(thumb_rgb)

    # Inverting makes tissue (dark) bright and background (white) dark, so the
    # thresholding algorithms treat tissue as the foreground.
    gray_image = color.rgb2gray(thumb_rgb)
    gray_smoothed = filters.gaussian(gray_image, sigma=sigma)
    complement_image = 1.0 - gray_smoothed

    threshold = THRESHOLD_METHODS[method_lower](complement_image)
    return complement_image > threshold
