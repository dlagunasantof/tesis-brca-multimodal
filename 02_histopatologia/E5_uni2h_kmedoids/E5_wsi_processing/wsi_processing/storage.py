"""Storage of WSI outputs and QC visualisation.

Persists quality-control reports as images and stores patches, coordinates and
embeddings in HDF5 files consolidated at the case (patient) level, so that a
patient with several slides ends up in a single file.
"""

import os
from pathlib import Path
from typing import List, Optional, Union

import h5py
import matplotlib
import numpy as np

matplotlib.use("Agg")           # headless backend: no display needed on a server
import matplotlib.pyplot as plt  # noqa: E402  (import after the backend is set)
from tqdm import tqdm            # noqa: E402

CHUNK_SIZE = 500                 # rows written per HDF5 write operation
QC_DPI = 150


def _tqdm_position() -> int:
    """Progress-bar row, so parallel workers do not overwrite each other's bars."""
    return int(os.environ.get("TQDM_POSITION", "0"))


def _slide_stem(wsi_filename: str) -> str:
    """Filename without extension, handling double extensions such as ``.ome.tiff``."""
    stem = Path(wsi_filename).stem
    if stem.endswith(".ome"):
        stem = Path(stem).stem
    return stem


def save_qc_report(
    thumb_rgb: np.ndarray,
    coarse_mask: np.ndarray,
    fine_mask: np.ndarray,
    output_dir: Union[str, Path],
    wsi_filename: str,
) -> Path:
    """Save a four-panel visual QC report for one slide.

    The panels are the original thumbnail, the coarse tissue mask (after pen
    removal), the fine QC mask (after the blur and contrast checks) and a green
    overlay of the accepted regions on the thumbnail. Reading them left to right
    shows what each filtering stage discarded.

    The file is named after the slide, with the extension replaced by ``.png``.

    Args:
        thumb_rgb: Original WSI thumbnail in RGB.
        coarse_mask: Tissue mask after macro filtering.
        fine_mask: Refined mask after micro quality filtering.
        output_dir: Destination directory, created if missing.
        wsi_filename: Original filename of the slide (e.g. ``case_name.svs``).

    Returns:
        Path of the written report.
    """
    out_dir_path = Path(output_dir)
    out_dir_path.mkdir(parents=True, exist_ok=True)
    report_path = out_dir_path / f"{_slide_stem(wsi_filename)}.png"

    if thumb_rgb.ndim == 3 and thumb_rgb.shape[-1] == 4:
        thumb_rgb = thumb_rgb[:, :, :3]

    fig, axes = plt.subplots(1, 4, figsize=(24, 6))

    axes[0].imshow(thumb_rgb)
    axes[0].set_title("1. Original Thumbnail")

    axes[1].imshow(coarse_mask, cmap="gray")
    axes[1].set_title("2. Coarse Tissue Mask")

    axes[2].imshow(fine_mask, cmap="gray")
    axes[2].set_title("3. Fine Quality Mask")

    # Green overlay marking only the accepted patch regions
    overlay = thumb_rgb.copy()
    selected = fine_mask.astype(bool)
    overlay[selected] = (
        overlay[selected] * 0.5 + np.array([0, 255, 0], dtype=np.uint8) * 0.5
    ).astype(np.uint8)

    axes[3].imshow(overlay)
    axes[3].set_title("4. Accepted Patch Regions (Overlay)")

    for ax in axes:
        ax.axis("off")

    plt.tight_layout()
    plt.savefig(str(report_path), dpi=QC_DPI)
    plt.close(fig)
    return report_path


def _write_in_chunks(dataset, data: np.ndarray, offset: int, desc: str) -> None:
    """Write ``data`` into ``dataset`` starting at ``offset``, in chunks of CHUNK_SIZE rows."""
    total = len(data)
    for start in tqdm(range(0, total, CHUNK_SIZE), desc=desc, unit="chunk",
                      leave=False, position=_tqdm_position()):
        end = min(start + CHUNK_SIZE, total)
        dataset[offset + start: offset + end] = data[start:end]


def _validate_lengths(
    coordinates: np.ndarray,
    slide_ids: List[str],
    patches: Optional[np.ndarray],
    embeddings: Optional[np.ndarray],
) -> int:
    """Check that every array carries one entry per patch; return that count."""
    num_new = coordinates.shape[0]
    if len(slide_ids) != num_new:
        raise ValueError(f"Mismatch between number of slide_ids ({len(slide_ids)}) and coordinates ({num_new})")
    if patches is not None and len(patches) != num_new:
        raise ValueError(f"Mismatch between number of patches ({len(patches)}) and coordinates ({num_new})")
    if embeddings is not None and embeddings.shape[0] != num_new:
        raise ValueError(f"Mismatch between number of embeddings ({embeddings.shape[0]}) and coordinates ({num_new})")
    return num_new


def save_case_data_to_h5(
    output_path: Union[str, Path],
    coordinates: np.ndarray,
    slide_ids: List[str],
    patches: Optional[np.ndarray] = None,
    embeddings: Optional[np.ndarray] = None,
) -> None:
    """Store or append patches, coordinates and embeddings in a case-level HDF5 file.

    If the file already exists — because another slide of the same case was
    processed earlier — the datasets are resized and the new rows appended, so
    every slide of a patient accumulates in one file. Patch images are stored
    with gzip compression to limit disk usage.

    Rows are aligned across datasets: row *i* of ``coords``, ``slide_ids``,
    ``patches`` and ``embeddings`` always refer to the same patch.

    Args:
        output_path: Target path of the case HDF5 file.
        coordinates: Coordinate matrix (P, 2) of int32, as ``[x, y]``.
        slide_ids: Slide filename of each patch, length P.
        patches: Patch images (P, patch_h, patch_w, 3) of uint8.
        embeddings: Feature matrix (P, embed_dim) of float32.

    Raises:
        ValueError: If the arrays do not agree on the number of patches, or if
            an existing file lacks a dataset that the new data provides.
    """
    dest_path = Path(output_path)
    dest_path.parent.mkdir(parents=True, exist_ok=True)

    num_new = _validate_lengths(coordinates, slide_ids, patches, embeddings)
    slide_ids_bytes = [sid.encode("utf-8") for sid in slide_ids]

    with h5py.File(str(dest_path), "a") as h5f:
        appending = "coords" in h5f
        current_size = h5f["coords"].shape[0] if appending else 0
        new_size = current_size + num_new

        if appending:
            # Refuse to append data the file has no dataset for: silently dropping
            # it would break the row alignment between datasets.
            if patches is not None and "patches" not in h5f:
                raise ValueError(f"{dest_path} has no 'patches' dataset; cannot append patches to it.")
            if embeddings is not None and "embeddings" not in h5f:
                raise ValueError(f"{dest_path} has no 'embeddings' dataset; cannot append embeddings to it.")
            if "patches" in h5f and patches is None:
                raise ValueError(f"{dest_path} holds patches; new data must provide patches as well.")
            if "embeddings" in h5f and embeddings is None:
                raise ValueError(f"{dest_path} holds embeddings; new data must provide embeddings as well.")

            h5f["coords"].resize(new_size, axis=0)
            h5f["slide_ids"].resize(new_size, axis=0)
            h5f["coords"][current_size:new_size] = coordinates
            h5f["slide_ids"][current_size:new_size] = slide_ids_bytes
        else:
            h5f.create_dataset("coords", shape=(num_new, 2), maxshape=(None, 2), dtype="int32")
            h5f.create_dataset("slide_ids", shape=(num_new,), maxshape=(None,),
                               dtype=h5py.string_dtype())
            h5f["coords"][:] = coordinates
            h5f["slide_ids"][:] = slide_ids_bytes

        if patches is not None:
            if appending:
                h5f["patches"].resize(new_size, axis=0)
            else:
                patch_h, patch_w, channels = patches.shape[1:]
                h5f.create_dataset(
                    "patches",
                    shape=(num_new, patch_h, patch_w, channels),
                    maxshape=(None, patch_h, patch_w, channels),
                    dtype="uint8",
                    compression="gzip",
                    compression_opts=4,
                )
            _write_in_chunks(h5f["patches"], patches, current_size, "Writing patches to H5")

        if embeddings is not None:
            if appending:
                h5f["embeddings"].resize(new_size, axis=0)
            else:
                h5f.create_dataset(
                    "embeddings",
                    shape=(num_new, embeddings.shape[1]),
                    maxshape=(None, embeddings.shape[1]),
                    dtype="float32",
                )
            _write_in_chunks(h5f["embeddings"], embeddings, current_size, "Writing embeddings to H5")


def save_embeddings_to_h5(
    output_path: Union[str, Path],
    embeddings: np.ndarray,
    coordinates: np.ndarray,
    patches: np.ndarray,
    slide_ids: List[str],
) -> None:
    """Backwards-compatible wrapper around :func:`save_case_data_to_h5`."""
    save_case_data_to_h5(
        output_path=output_path,
        coordinates=coordinates,
        slide_ids=slide_ids,
        patches=patches,
        embeddings=embeddings,
    )
