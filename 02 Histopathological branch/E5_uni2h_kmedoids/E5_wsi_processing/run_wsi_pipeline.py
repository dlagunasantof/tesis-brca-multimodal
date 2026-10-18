#!/usr/bin/env python
"""End-to-end WSI preprocessing and UNI2-h feature extraction (Approach E5).

For every slide in the input directory the script:

1. opens the slide and reads a low-resolution thumbnail;
2. segments tissue against the bright background;
3. removes pen markings and smooths the mask morphologically;
4. walks the slide on a regular grid at the target magnification, keeping the
   tiles sufficiently covered by tissue that pass the blur and contrast checks;
5. optionally embeds the kept patches with UNI2-h;
6. appends coordinates, slide identifiers, patches and embeddings to a
   **case-level** HDF5 file, and writes a four-panel QC report.

Slides are grouped into cases by the first ``--case-id-length`` characters of the
filename (12 for a TCGA barcode such as ``TCGA-XX-XXXX``), so several slides of
one patient accumulate in a single HDF5 file.

Examples
--------
Patches only, no model (no GPU required)::

    python run_wsi_pipeline.py --wsi-dir slides/ --out-dir output/ --no-embeddings

Full pipeline with embeddings, discarding the patch pixels to save disk::

    python run_wsi_pipeline.py --wsi-dir slides/ --out-dir output/ --no-patches

Outputs
-------
``<out-dir>/h5/<case_id>.h5``   coords, slide_ids, patches, embeddings
``<out-dir>/qc/<slide>.png``    four-panel QC report
``<out-dir>/pipeline_summary.csv``  one row per slide, with the patch counts
"""

import argparse
import sys
import traceback
from pathlib import Path

import pandas as pd
from tqdm import tqdm

from wsi_processing import (
    extract_valid_patches,
    filter_tissue_artifacts,
    get_raw_tissue_mask,
    get_wsi_thumbnail,
    open_wsi,
    save_case_data_to_h5,
    save_qc_report,
)

SLIDE_EXTENSIONS = ("*.svs", "*.ndpi", "*.tif", "*.tiff", "*.mrxs", "*.scn")


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__.split("\n")[0],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    io_group = parser.add_argument_group("input / output")
    io_group.add_argument("--wsi-dir", type=Path, required=True, help="directory holding the slides")
    io_group.add_argument("--out-dir", type=Path, required=True, help="destination directory")
    io_group.add_argument("--pattern", default=None,
                          help="glob pattern of the slide files (default: every supported extension)")
    io_group.add_argument("--case-id-length", type=int, default=12,
                          help="characters of the filename that identify the case")
    io_group.add_argument("--overwrite", action="store_true",
                          help="reprocess slides whose case file already exists")

    patch_group = parser.add_argument_group("patch extraction")
    patch_group.add_argument("--patch-size", type=int, default=256, help="patch side in pixels")
    patch_group.add_argument("--magnification", type=float, default=20.0,
                             help="objective power at which patches are read")
    patch_group.add_argument("--overlap", type=float, default=0.0,
                             help="overlap fraction between adjacent tiles, in [0, 1)")
    patch_group.add_argument("--tissue-threshold", type=float, default=0.8,
                             help="minimum tissue ratio of a tile")
    patch_group.add_argument("--sharpness-threshold", type=float, default=0.0005,
                             help="minimum variance of the Laplacian (blur rejection)")
    patch_group.add_argument("--contrast-threshold", type=float, default=0.05,
                             help="low-contrast tolerance")
    patch_group.add_argument("--thumbnail-size", type=int, default=2048,
                             help="largest side of the thumbnail used for segmentation")
    patch_group.add_argument("--threshold-method", default="triangle", choices=["triangle", "otsu"],
                             help="thresholding method for tissue segmentation")
    patch_group.add_argument("--no-marker-filter", action="store_true",
                             help="skip pen-mark detection")

    model_group = parser.add_argument_group("feature extraction")
    model_group.add_argument("--no-embeddings", action="store_true",
                             help="extract patches only, without running UNI2-h")
    model_group.add_argument("--no-patches", action="store_true",
                             help="do not store the patch pixels (embeddings only)")
    model_group.add_argument("--weights-path", default=None,
                             help="local UNI2-h checkpoint; omit to use the Hugging Face Hub")
    model_group.add_argument("--device", default=None, help="'cuda' or 'cpu' (auto-detected if omitted)")
    model_group.add_argument("--batch-size", type=int, default=64, help="patches per forward pass")

    args = parser.parse_args(argv)
    if args.no_embeddings and args.no_patches:
        parser.error("--no-embeddings and --no-patches together would store nothing.")
    if not 0.0 <= args.overlap < 1.0:
        parser.error("--overlap must be in [0, 1).")
    return args


def find_slides(wsi_dir: Path, pattern):
    """Return the sorted list of slide files to process."""
    if pattern:
        return sorted(wsi_dir.glob(pattern))
    slides = []
    for ext in SLIDE_EXTENSIONS:
        slides.extend(wsi_dir.glob(ext))
    return sorted(set(slides))


def process_slide(slide_path: Path, args, extractor):
    """Run the whole pipeline on one slide and persist its outputs.

    Returns:
        A dictionary with the per-slide summary of the run.
    """
    tile_source = open_wsi(slide_path)

    thumb_rgb, rx, ry = get_wsi_thumbnail(tile_source, target_size=args.thumbnail_size)
    raw_mask = get_raw_tissue_mask(thumb_rgb, threshold_method=args.threshold_method)
    clean_mask = filter_tissue_artifacts(
        raw_mask, thumb_rgb, filter_markers=not args.no_marker_filter
    )

    coordinates, patches, slide_ids, quality_mask = extract_valid_patches(
        tile_source=tile_source,
        raw_tissue_mask=clean_mask,
        rx=rx,
        ry=ry,
        slide_id=slide_path.name,
        patch_width=args.patch_size,
        patch_height=args.patch_size,
        target_magnification=args.magnification,
        overlap=args.overlap,
        tissue_threshold=args.tissue_threshold,
        sharpness_threshold=args.sharpness_threshold,
        contrast_threshold=args.contrast_threshold,
    )

    save_qc_report(thumb_rgb, clean_mask, quality_mask, args.out_dir / "qc", slide_path.name)

    n_patches = len(coordinates)
    embeddings = None
    if n_patches and extractor is not None:
        embeddings = extractor.extract_features(patches, batch_size=args.batch_size)

    if n_patches:
        case_id = slide_path.stem[: args.case_id_length]
        save_case_data_to_h5(
            output_path=args.out_dir / "h5" / f"{case_id}.h5",
            coordinates=coordinates,
            slide_ids=slide_ids,
            patches=None if args.no_patches else patches,
            embeddings=embeddings,
        )

    return {
        "slide": slide_path.name,
        "case_id": slide_path.stem[: args.case_id_length],
        "tissue_pct_thumbnail": round(100 * float(clean_mask.mean()), 2),
        "accepted_pct_thumbnail": round(100 * float(quality_mask.mean()), 2),
        "n_patches": n_patches,
        "embed_dim": int(embeddings.shape[1]) if embeddings is not None else 0,
        "status": "ok" if n_patches else "no_valid_patches",
    }


def main(argv=None) -> int:
    args = parse_args(argv)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    slides = find_slides(args.wsi_dir, args.pattern)
    print(f"Slides found: {len(slides)}")
    if not slides:
        print(f"No slides matched in {args.wsi_dir}", file=sys.stderr)
        return 1

    if not args.overwrite:
        pending = [s for s in slides
                   if not (args.out_dir / "h5" / f"{s.stem[:args.case_id_length]}.h5").exists()]
        if len(pending) < len(slides):
            print(f"Skipping {len(slides) - len(pending)} slide(s) whose case file already exists "
                  f"(use --overwrite to reprocess).")
        slides = pending

    if not slides:
        print("Nothing to do: every slide already has a case file. Use --overwrite to reprocess.")
        return 0

    # The model is loaded once and reused across slides
    extractor = None
    if not args.no_embeddings:
        from wsi_processing import UNI2FeatureExtractor
        extractor = UNI2FeatureExtractor(weights_path=args.weights_path, device=args.device)

    records = []
    for slide_path in tqdm(slides, desc="Slides", unit="slide"):
        try:
            records.append(process_slide(slide_path, args, extractor))
        except Exception as exc:
            print(f"\nERROR on {slide_path.name}: {type(exc).__name__}: {exc}", file=sys.stderr)
            traceback.print_exc(file=sys.stderr)
            records.append({"slide": slide_path.name,
                            "case_id": slide_path.stem[: args.case_id_length],
                            "tissue_pct_thumbnail": 0.0, "accepted_pct_thumbnail": 0.0,
                            "n_patches": 0, "embed_dim": 0,
                            "status": f"error: {type(exc).__name__}"})

    summary = pd.DataFrame(records)
    summary_path = args.out_dir / "pipeline_summary.csv"
    summary.to_csv(summary_path, index=False)

    n_ok = int((summary["status"] == "ok").sum())
    print(f"\nProcessed: {n_ok} OK / {len(summary) - n_ok} failed or empty")
    print(f"Total patches: {int(summary['n_patches'].sum()):,}")
    print(f"Summary written to {summary_path}")
    return 0 if n_ok == len(summary) else 1


if __name__ == "__main__":
    sys.exit(main())
