# `wsi_processing` — WSI preprocessing and UNI2-h feature extraction

Independent pipeline used for **Approach E5 (UNI2-h + K-Medoids)**. Unlike the other approaches, patches are extracted from scratch directly at the target magnification, following the tiling methodology of the UNI2-h authors, rather than being cropped from the shared 1000 × 1000 px patches.

## Layout

```
wsi_processing/
├── __init__.py        public API (the model wrapper is imported lazily)
├── utils.py           slide opening, metadata, thumbnails, scale factors
├── segmentation.py    tissue segmentation against the bright background
├── artifacts.py       pen-mark removal, morphology, patch-level QC
├── storage.py         QC reports and case-level HDF5 persistence
└── models.py          UNI2-h feature extractor (torch + timm)
run_wsi_pipeline.py    command-line driver: slides in, HDF5 + QC reports out
```

## Installation

```bash
pip install large-image[openslide,tiff] numpy scikit-image h5py matplotlib pandas tqdm
pip install torch torchvision timm python-dotenv        # only for feature extraction
```

`large_image` needs the OpenSlide system library (`apt install openslide-tools` on Debian/Ubuntu). UNI2-h weights are gated on the Hugging Face Hub: accept the model licence and export `HF_TOKEN` (a `.env` file is read automatically), or pass `--weights-path` to load a local checkpoint.

## Usage

The command-line driver replaces the notebooks: one command processes a whole directory of slides.

```bash
# Patches only, no model (no GPU required)
python run_wsi_pipeline.py --wsi-dir slides/ --out-dir output/ --no-embeddings

# Full pipeline; store embeddings but discard patch pixels to save disk
python run_wsi_pipeline.py --wsi-dir slides/ --out-dir output/ --no-patches

# Reprocess everything, at 40x, with 256 px patches and 25 % overlap
python run_wsi_pipeline.py --wsi-dir slides/ --out-dir output/ \
    --magnification 40 --patch-size 256 --overlap 0.25 --overwrite
```

`python run_wsi_pipeline.py --help` lists every parameter (tissue, sharpness and contrast thresholds, thumbnail size, thresholding method, batch size, device).

### Outputs

| Path | Content |
|---|---|
| `<out-dir>/h5/<case_id>.h5` | `coords` (P, 2) int32, `slide_ids` (P,), `patches` (P, h, w, 3) uint8 gzip, `embeddings` (P, 1536) float32 |
| `<out-dir>/qc/<slide>.png` | Four panels: thumbnail, coarse mask, fine QC mask, accepted regions |
| `<out-dir>/pipeline_summary.csv` | One row per slide: tissue %, accepted %, patch count, status |

Rows are aligned across datasets: row *i* of `coords`, `slide_ids`, `patches` and `embeddings` refers to the same patch. Slides are grouped into cases by the first `--case-id-length` characters of the filename (12 for a TCGA barcode), so several slides of one patient accumulate in a single file.

The driver skips slides whose case file already exists, so an interrupted run can be resumed by re-issuing the same command; pass `--overwrite` to force reprocessing.

### Library use

```python
from wsi_processing import (open_wsi, get_wsi_thumbnail, get_raw_tissue_mask,
                            filter_tissue_artifacts, extract_valid_patches)

tile_source = open_wsi("slide.svs")
thumb, rx, ry = get_wsi_thumbnail(tile_source)
clean_mask = filter_tissue_artifacts(get_raw_tissue_mask(thumb), thumb)
coords, patches, slide_ids, qc_mask = extract_valid_patches(
    tile_source, clean_mask, rx, ry, slide_id="slide.svs")
```

Use `artifact_filtering` instead of `extract_valid_patches` when only the QC mask is needed: it applies the same filtering without keeping the patch pixels in memory.

## Pipeline stages

1. **Thumbnail and scale factors** (`get_wsi_thumbnail`) — the slide is read at low resolution; `rx`, `ry` map native coordinates to thumbnail space.
2. **Tissue segmentation** (`get_raw_tissue_mask`) — grayscale, Gaussian smoothing, inversion, then triangle or Otsu thresholding.
3. **Macro filtering** (`filter_tissue_artifacts`) — pen strokes (red, green, blue, black) are removed by RGB thresholds, then binary opening and erosion clean the mask.
4. **Micro filtering** (`extract_valid_patches`) — the slide is walked on a regular grid; tissue coverage is checked on the mask *before* any disk read, and tiles are then rejected if low-contrast or blurred (variance of the Laplacian).
5. **Feature extraction** (`UNI2FeatureExtractor`) — patches are resized to 224 × 224, normalised with the ImageNet statistics and embedded into 1536 dimensions.
6. **Persistence** (`save_case_data_to_h5`, `save_qc_report`).

## Changes made when reorganising

Behaviour is unchanged: on synthetic slides, thumbnails, tissue masks, cleaned masks, patch pixels, coordinates and QC masks are identical to those of the original code.

* `artifact_filtering` and `extract_valid_patches` were nearly identical copies of the same algorithm; the shared logic now lives in one private generator, and `artifact_filtering` is a thin wrapper over it.
* Imports are relative (`from .utils import ...`), so the package works regardless of where it sits on `sys.path`; `tqdm` and `os` moved to module level.
* `models.py` is imported lazily, so the package can be used without torch or timm installed.
* Fixed in `save_qc_report`: the fourth panel was indexed as `axes[4-1]`, and the docstring promised a `_qc.png` suffix the code never added (the code was kept, the docstring corrected).
* `save_case_data_to_h5` now refuses to append data whose dataset does not exist in the target file. Previously it created the dataset sized for the full file and wrote only the new rows, silently leaving the earlier rows as zeros and breaking the alignment between datasets.
* The mutable default argument `colors_to_filter=["red", ...]` became a tuple constant.
* `matplotlib` is set to the `Agg` backend, so QC reports render on a headless server.
* `run_wsi_pipeline.py` is new: it replaces the notebook-driven execution with one resumable command, and writes a per-slide summary.

## Known limitations

1. **Memory.** `extract_valid_patches` holds every accepted patch of a slide in RAM; a large slide at high magnification can need several GB. Use `--no-patches` (embeddings only) or `artifact_filtering` where possible.
2. **Incomplete edge tiles.** The grid is built with `range(0, native_w - w_native, stride_x)`, so the right and bottom edges of the slide are not tiled when the dimensions are not an exact multiple of the stride. This is the original behaviour and was left unchanged.
3. **Silent tile skips.** Read errors, wrong tile shapes and failed quality checks are all skipped without being counted separately, so the summary reports how many patches were kept but not why the rest were dropped.
4. **Magnification fallbacks.** When the slide metadata lacks the objective power and the pixel spacing, `get_wsi_metadata` assumes 20x and prints a warning. Since the patch size in native pixels is derived from it, a wrong assumption changes the physical area each patch covers. Check the warnings in the log.
