# m2nGAN

`m2nGAN` generates pseudo-nucleus fluorescence volumes from membrane fluorescence images for downstream cell-lineage analysis.

> **Important: do not skip the 3D NIfTI-to-2D TIFF conversion.**
> Model predictions are stored as 3D `.nii.gz` volumes, whereas the downstream lineage/ACEtree workflow expects a sequence of correctly named 2D ImageJ TIFF files. Run [`3dniigz_to_2dtif.py`](3dniigz_to_2dtif.py) after prediction and before lineage tracing. The actual conversion implementation is [`Utils/niigz3d_to_tiff2d.py`](Utils/niigz3d_to_tiff2d.py).

## Environment and reproducibility

Use **64-bit Python 3.7**. This is a legacy PyTorch/SimpleITK project; Python 3.8+—and especially the locally installed Python 3.13—must not be assumed to work with its historical dependency stack. A different Python interpreter can silently change image interpolation, TIFF metadata handling, or prevent binary packages from installing.

Create a clean environment and install the locked runtime dependencies:

```bash
conda create -n m2ngan python=3.7
conda activate m2ngan
python -m pip install --upgrade "pip<24.1"
python -m pip install -r requirements.txt
```

`requirements.txt` is part of the experiment configuration, not an optional convenience file. In particular, do not upgrade `numpy`, `scikit-image`, `tifffile`, `SimpleITK`, `torch`, or `torchvision` independently. `torch`/`torchvision` must be installed as a matched pair. For GPU execution, install the matching PyTorch 1.10.2 CUDA build from the PyTorch archive first, then install the remaining packages from `requirements.txt`.

After installation, verify the conversion-specific packages before processing a full experiment:

```bash
python -c "import sys, nibabel, numpy, skimage, tifffile; print(sys.version); print('nibabel', nibabel.__version__); print('numpy', numpy.__version__); print('scikit-image', skimage.__version__); print('tifffile', tifffile.__version__)"
```

## 3D `.nii.gz` to 2D TIFF: required export step

The entry point [`3dniigz_to_2dtif.py`](3dniigz_to_2dtif.py) batch-converts generated 3D nuclei volumes. It calls `seperate_3dniigz_to_2dtif()` in `Utils/niigz3d_to_tiff2d.py`, which:

1. loads each NIfTI volume with `nibabel`;
2. smooths it with `skimage.filters.gaussian(sigma=1.5)`;
3. resizes the XY plane with `skimage.transform.resize`, preserving the generated Z-depth;
4. samples and reverses the Z order to match the raw-image page convention;
5. writes one unsigned-8-bit ImageJ TIFF per page through `tifffile.imwrite`, including byte order, pixel resolution, size, and ImageJ metadata; and
6. names files as `<embryo>_L1-t<time>-p<page>.tif`, which is the downstream lookup convention.

Changing any of these conventions can produce TIFFs that open normally but are spatially misaligned, have inverted page order, or are not detected by lineage software. The versions of `scikit-image` and `tifffile` are therefore pinned: they control resampling behavior and serialized TIFF/ImageJ metadata, respectively.

### Configure and run

Edit the configuration block at the top of [`3dniigz_to_2dtif.py`](3dniigz_to_2dtif.py) for your dataset:

| Setting | Meaning |
| --- | --- |
| `embryo_names` | Embryo folder names to process. |
| `maxtimes` | Number of time points for each embryo. |
| `raw_xyzs` | Target raw-image shape in `(X, Y, Z)` order. These values define TIFF dimensions and page count. |
| `generativeNuc_dir_tem` | Root containing `<embryo>/RawNuc_m2nGAN_prediction/*.nii.gz`. |
| `saving_out_raw_root` | Root for the exported `<embryo>/tif/*.tif` files. |

Each input must follow this naming pattern:

```text
<prediction-root>/<embryo>/RawNuc_m2nGAN_prediction/<embryo>_<time>_rawNuc.nii.gz
```

Then run it from the repository root:

```bash
python 3dniigz_to_2dtif.py
```

The exporter is restart-safe: existing TIFFs are skipped and reported. Before starting lineage tracing, check one time point visually and confirm that (1) the TIFF dimensions equal the configured X/Y values, (2) there are exactly Z pages, (3) page numbering runs from `p01` to `pZZ`, and (4) the slice order aligns with the raw membrane stack.

## Fluorescence quantification

[`calculate_cellular_expression.py`](calculate_cellular_expression.py) measures the fluorescence signal of a protein of interest in every segmented cell at every time point. All paths are given on the command line.

### Inputs

| Input | Location and naming | Needed |
| --- | --- | --- |
| Cell segmentation | `<seg-root>/<embryo>/SegCell/<embryo>_<tp>_segCell.nii.gz` | always |
| Fluorescence image | `<expression-root>/<embryo>/Gene/<embryo>_<tp>_expression.nii.gz` | always |
| Lineage file | `<cd-root>/CD<embryo>.csv` | only for cell-name columns |
| Name dictionary | `name_dictionary.csv` (segmentation label to cell name) | only for cell-name columns |

- `<tp>` is the zero-padded time point (`001`, `002`, ...), numbered consecutively from `001`.
- The segmentation and the fluorescence image of a time point must have the same shape.
- The segmentations are the output of the cell segmentation pipeline ([CTransformer](https://github.com/chiellini/CTransformer)). If they are stored in a folder other than `SegCell`, pass `--seg-dir-name`.
- The lineage file and the name dictionary are optional. They are only used to label the output columns with cell names.

### Run

```bash
# Columns labelled with cell names
python calculate_cellular_expression.py \
    --seg-root /path/to/segmentation \
    --expression-root /path/to/expression \
    --cd-root /path/to/CDFiles \
    --name-dictionary /path/to/name_dictionary.csv \
    --output /path/to/ExpressionStat

# No lineage files: columns labelled with segmentation labels
python calculate_cellular_expression.py \
    --seg-root /path/to/segmentation \
    --expression-root /path/to/expression \
    --output /path/to/ExpressionStat
```

Optional arguments: `--embryo NAME` (repeatable; default is every embryo under `--seg-root`), `--workers N` for parallel processing, and `--assemble-only` to rebuild the tables from existing intermediate results. The script checks that all input files exist before the calculation starts.

### Outputs

Three tables per embryo are written to `--output`. Rows are time points.

| File | Content |
| --- | --- |
| `<embryo>_cellwise_expression.csv` | Summed fluorescence intensity within each cell. |
| `<embryo>_surfacewise_expression.csv` | Summed fluorescence intensity in the shell around each cell surface. |
| `<embryo>_contactwise_expression.csv` | Summed fluorescence intensity at each cell-cell contact (two header rows, `cell1` and `cell2`). |

Without lineage files the three tables end in `_by_label.csv`; cells are identified by segmentation label and contacts by `<label1>_<label2>`. Per-time-point intermediate results are kept in `<output>/<embryo>/`.

Expected outputs: the tables obtained for three embryos (`200710hmr1plc1p1`, `200710hmr1plc1p2`, `200710hmr1plc1p3`) are provided in [`Data/HMR1StyleExpressionStatHMR1/`](Data/HMR1StyleExpressionStatHMR1).

## Main project scripts

| Script | Purpose |
| --- | --- |
| `train.py` | Train the m2nGAN model. |
| `test.py` | Run patch-based inference on a 3D input. |
| `process_for_m2nGAN.py` | Prepare 3D training/testing volumes from segmentation data. |
| `3dniigz_to_2dtif.py` | **Required bridge** from predicted 3D NIfTI volumes to 2D TIFF lineage inputs. |
| `calculate_cellular_expression.py` | Quantify per-cell, cell-surface and cell-contact fluorescence (see [Fluorescence quantification](#fluorescence-quantification)). |
| `generate_lineage_tree_files.py` | Generate lineage-tree files. |

## Notes for future changes

- Keep the Python version and `requirements.txt` together with every experiment result.
- Test any dependency update on a known NIfTI input and compare TIFF shape, dtype, page order, metadata, and visual alignment before using it for a full lineage run.
- Preserve the TIFF filename template and `raw_xyzs` order unless all downstream lineage readers are updated at the same time.
