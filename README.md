# nf-azimuth

A Nextflow workflow for single-cell RNA-seq cell type annotation using Azimuth reference datasets.

## Overview

This pipeline performs automated cell type annotation of single-cell RNA-seq data using the Azimuth algorithm. It can process data against reference datasets and optionally generate interactive visualizations using Vitessce.

## Quick Start

### Prerequisites

- [Nextflow](https://www.nextflow.io/) (≥21.04.0)
- [Docker](https://www.docker.com/) or [Singularity](https://sylabs.io/singularity/)
- [Reference Files](https://azimuth.hubmapconsortium.org/references/) (use Zenodo link to download)

### Basic Usage

```bash
nextflow run main.nf \
  --input /path/to/input.h5 \
  --reference /path/to/reference_folder \
  --output /path/to/output/
```

### With Interactive Visualization

```bash
nextflow run main.nf \
  --input /path/to/input.h5 \
  --reference /path/to/reference_folder \
  --output /path/to/output/ \
  --web_output /path/to/web/output/
```

## Parameters

### Required Parameters

| Parameter | Description | Example |
|-----------|-------------|---------|
| `--input` | Input single-cell data file (H5 Matrix or H5Seurat object) | `/data/sample.h5seurat` |
| `--reference` | Azimuth reference dataset file path | `/data/human_pbmc` |
| `--output` | Output directory for results | `/results/` |

### Optional Parameters

| Parameter | Description | Default | Example |
|-----------|-------------|---------|---------|
| `--web_output` | Directory for web visualization files | `false` | `/web/output/` |
| `--name` | Analysis name | `"Azimuth Analysis"` | `"My Analysis"` |
| `--description` | Analysis description | `"Mapping of gene expression data onto Azimuth reference"` | `"Custom description"` |

## Output Files

### Standard Output

- `*.h5seurat` - Annotated single-cell data in H5Seurat format
- `*.h5ad` - Annotated single-cell data in H5AD format (if conversion enabled)

### Web Visualization Output (when `--web_output` is specified)

- `vitessce.config.*.json` - Vitessce configuration file
- `chart.manifest.json` - Chart manifest for Cirro visualization
- `azimuth.zarr/` - Zarr-formatted data for web viewing
