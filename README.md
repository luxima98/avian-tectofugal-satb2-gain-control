# Asymmetric Feedforward–Feedback Architecture and SATB2-Mediated Gain Control in the Avian Tectofugal Pathway

This repository contains analysis scripts and selected downstream results associated with the study **“Asymmetric Feedforward–Feedback Architecture and SATB2-Mediated Gain Control in the Avian Tectofugal Pathway.”**

- Each top-level folder corresponds to a specific experimental/analytical component of the study and contains the relevant scripts and generated outputs.
- For the bioinformatics and orthology components, this repository primarily stores **downstream artifacts** generated from publicly available raw data (see *Data provenance and statement*).

---

## Data provenance and statement (important)

The **bioinformatics analyses** and **orthology (cross-species mapping) analyses** in this repository are based on raw data originally released with:

- **Science (2024)**, DOI: `10.1126/science.adp5182`
- **Heidelberg Open Research Data (Heidelberg University)** dataset:
  - https://heidata.uni-eidelberg.de/dataset.xhtml?persistentId=doi:10.11588/DATA/BX6REK

The upstream raw data include (but are not limited to):

- **ISS (in situ sequencing) raw data**
- **Adult chicken single-cell/single-nucleus transcriptomic raw data**
- **Spatial in situ hybridization / spatial datasets**

This repository **does not redistribute the complete upstream raw datasets**. Instead, it contains selected **derived outputs** generated from analyses of those raw data, including:

- **NeuronChat** downstream `.rds` results
- Cell/spot **spatial coordinate tables** (e.g., `all_cells_coordinates.csv`)
- Intermediate tables and visualization outputs required for reproducing analysis steps

To reproduce the full pipeline from scratch, please download the raw data from the link above and update the absolute paths in the scripts to match your local environment.

---

## Repository structure

### `dendrite/`
Morphological processing and statistics for dendrites and spines.

- `process_dendrite_data.py`
  - Reads multiple `*.csv` files in the folder (morphometric exports)
  - Summarizes dendrite length/width, spine counts, mean/total spine length
  - Outputs:
    - `dendrite_analysis_results*.csv`
    - `dendrite_spine_details*.csv`
- `spine_classification.py`
  - Spine classification workflow (see *Dependencies*).
- `spine_density_analysis.py`
  - Spine density analysis.
- `dendrite_analysis_results_*.csv`
  - Example summary outputs.

### `rf_analysis/`
Receptive-field (RF) analysis summary, statistics, and visualization.

- `RF_MVL_E.xlsx`
  - RF-derived measurements/results for MVL and E nuclei.
- `visualize_mvl_e_analysis.py`
  - Reads the `MVL_分析结果` and `E_分析结果` sheets from `RF_MVL_E.xlsx`
  - Runs statistical tests (t-test / Welch / Mann–Whitney U) and produces summary figures
  - Outputs figures and a statistics CSV
- `more_visualizations.py`
  - Additional visualization utilities.
- `statistics_info.md`
  - Text summary of key statistical outcomes.

### `ISS/`
Spatial workflows for ISS data, including segmentation and mapping (Baysor + Tangram).

- `run_baysor_tangram_pipeline.py`
  - Converts spatial `h5ad` data into Baysor-compatible transcript tables
  - Calls **Julia + Baysor** for segmentation (script contains example local Julia/Baysor paths)
  - Generates Tangram input matrices, including:
    - `cell_coordinates.csv`
    - `cell_gene_expression.csv`
    - `spatial_cells.h5ad`
    - `cell_statistics.csv`
    - `segmentation_summary.csv`
- `run_tangram_official_flow.py`
  - Implements a Tangram mapping flow aligned with Tangram’s official notebooks (`map_cells_to_space`).
- `visualize_tangram_results.py`
  - Tangram result visualization.
- `run_cell_type_analysis.R`
  - Constructs NeuronChat objects from Seurat objects and runs groupwise analyses.
- `all_cells_coordinates.csv`
  - Spatial coordinate table (derived output).

### `orthology/`
Chicken–mouse orthology mapping and cross-species conversion utilities.

- `chicken_to_mouse_gene_conversion.R`
  - Builds chicken→mouse gene mappings (e.g., from OrthoFinder outputs such as `G.gallus__v__M.musculus.tsv`)
  - Optionally uses `biomaRt` to convert mouse Ensembl IDs to gene symbols
  - Produces `gene_conversion_mapping.csv` and converted matrices
- `complete_gene_conversion.R`
  - Extended conversion workflow: supports reuse of existing mappings, duplicate-resolution strategies (weighted average / best match), and `SCTransform`.
- `create_chicken_orth_mice_seurat.R`
  - Builds a cross-species converted Seurat object for downstream NeuronChat/CellChat analyses.
- `*_confidence_interactions.csv` / `neuronchat_best_hit_interactions.csv`
  - Stratified and best-hit interaction tables exported from NeuronChat post-processing.

### `neuronchat/`
NeuronChat inference, post-processing, and visualization.

- `complete_neuronchat_analysis.R`
  - Loads the cross-species Seurat object
  - Runs NeuronChat inference (`run_NeuronChat`)
  - Incorporates mapping confidence and expression support to weight and stratify interactions
  - Outputs `neuronchat_analysis_result.rds` and multiple `*.csv` tables
- `neuronchat_updated_visualization.R`
  - Visualization utilities for NeuronChat outputs.
- `E_MVL_integrated_module_analysis_results.R`
  - Integrated/module analysis script related to E/MVL.
- `neuronchat_analysis_result.rds`
  - Serialized NeuronChat result object (derived output).

---

## Software environment and dependencies

### Python packages (as used in this repository)

- `pandas`
- `numpy`
- `matplotlib`
- `seaborn`
- `scipy`
- `scanpy`
- `anndata`
- `tangram`
- `h5py`
- `Pillow` (`PIL`)

Notes:

- `ISS/run_baysor_tangram_pipeline.py` invokes **Julia + Baysor** via `subprocess`.

### R packages (as used in this repository)

- `Seurat`
- `SeuratObject` (layer-based APIs in Seurat v5)
- `dplyr`
- `NeuronChat`
- `CellChat`
- `biomaRt` (optional; used for Ensembl ID → gene symbol conversion)

### External tools / path dependencies

- **Julia** (example in scripts: `Julia-1.10.10`)
- **Baysor** (Julia package; scripts include example local paths to a Baysor project directory)

Several scripts contain author-specific absolute paths (e.g., `D:\r\data\...`). You will need to update these paths to match your local filesystem.

---

## Quick start (recommended)

### 1) Obtain the repository

Clone or download this repository to any local directory.

### 2) Obtain upstream raw data

- Download the upstream raw data from the heidata record listed above.
- Place the raw data in a local directory of your choice.
- Edit scripts to replace absolute paths (e.g., `base_dir`, `transcript_file`, `sc_file`) with your local paths.

### 3) Run example workflows

- **RF visualization**: run `rf_analysis/visualize_mvl_e_analysis.py`
- **Dendrite/spine summarization**: run `dendrite/process_dendrite_data.py` inside `dendrite/` (ensure the target `*.csv` inputs are present)
- **ISS → Baysor → Tangram**: see `ISS/run_baysor_tangram_pipeline.py` and `ISS/run_tangram_official_flow.py`
- **Orthology conversion + NeuronChat**:
  - run gene conversion scripts in `orthology/` to generate converted matrices/Seurat objects
  - then run `neuronchat/complete_neuronchat_analysis.R`

---

## Reproducibility notes

- **Paths**: many scripts are not parameterized and contain absolute paths; please edit them for your environment.
- **Runtime/compute**: Baysor segmentation, Tangram mapping, and NeuronChat inference can be computationally intensive.
- **Versioning**: Seurat v4 vs v5 differ in APIs (`GetAssayData`, `Layers()`, etc.). This repository includes partial compatibility handling, but a pinned environment is recommended.

---

## Citation

If you use the scripts and/or derived outputs in this repository, please cite the upstream data/source publication:

- Science (2024). DOI: `10.1126/science.adp5182`
- Heidelberg Open Research Data (Heidelberg University): https://heidata.uni-eidelberg.de/dataset.xhtml?persistentId=doi:10.11588/DATA/BX6REK
