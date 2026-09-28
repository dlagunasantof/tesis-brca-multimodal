"""WSI preprocessing and feature-extraction pipeline.

Modules:
    segmentation: Raw tissue isolation from slide thumbnails.
    artifacts:    Pen-mark removal, morphological cleaning and patch-level QC.
    utils:        Slide opening, metadata parsing and scaling utilities.
    storage:      QC reports and case-level HDF5 persistence.
    models:       UNI2-h foundation-model feature extractor.

Typical use, per slide::

    from wsi_processing import (open_wsi, get_wsi_thumbnail, get_raw_tissue_mask,
                                filter_tissue_artifacts, extract_valid_patches)

    tile_source = open_wsi(path)
    thumb, rx, ry = get_wsi_thumbnail(tile_source)
    raw_mask = get_raw_tissue_mask(thumb)
    clean_mask = filter_tissue_artifacts(raw_mask, thumb)
    coords, patches, slide_ids, qc_mask = extract_valid_patches(
        tile_source, clean_mask, rx, ry, slide_id=path.name)

``models`` requires torch and timm, which the other modules do not, so it is
imported lazily: ``UNI2FeatureExtractor`` is resolved on first access.
"""

from typing import TYPE_CHECKING

from .artifacts import (
    artifact_filtering,
    extract_valid_patches,
    filter_tissue_artifacts,
    marker_detection,
)
from .segmentation import get_raw_tissue_mask
from .storage import (
    save_case_data_to_h5,
    save_embeddings_to_h5,
    save_qc_report,
)
from .utils import (
    get_wsi_metadata,
    get_wsi_thumbnail,
    open_wsi,
)

if TYPE_CHECKING:                       # for type checkers and IDEs only
    from .models import UNI2FeatureExtractor

__all__ = [
    "get_raw_tissue_mask",
    "marker_detection",
    "filter_tissue_artifacts",
    "artifact_filtering",
    "extract_valid_patches",
    "open_wsi",
    "get_wsi_metadata",
    "get_wsi_thumbnail",
    "save_qc_report",
    "save_case_data_to_h5",
    "save_embeddings_to_h5",
    "UNI2FeatureExtractor",
]


def __getattr__(name: str):
    """Import the torch-dependent model wrapper only when it is requested."""
    if name == "UNI2FeatureExtractor":
        from .models import UNI2FeatureExtractor as _UNI2FeatureExtractor
        return _UNI2FeatureExtractor
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
