# SpatialHub - Pipelines for Spatial Transcriptomics Analysis

> **⚠️ WARNING: SpatialHub is currently in pre-release and under active development. Although we have taken great care to ensure its quality and reliability, bugs and breaking changes may still occur. Please use it with caution and report any issues by opening a GitHub issue.**

## Installation

```bash
uv venv -p 3.12 spatialhub-py312
source spatialhub-py312/bin/activate
python --version  # Python 3.12.12
cd /path/to/spatialhub_dev/
uv pip install -e .
spatialhub --version
```

## Workflows

### CosMx makeZarr

This workflow generates one Zarr file per CosMx sample from raw data. It supports two branches depending on whether Ashlar is used.

- With Ashlar: cosmx_genBlankFOV, cosmx_completeGrid, cosmx_runAshlar, cosmx_makeZarr.
- Without Ashlar: cosmx_assembleFOVs, cosmx_makeZarr.

## Visium HD makeZarr

This workflow generates one Zarr file per Visium HD capture area from FASTQ files and images. It consists of two steps:

- Space Ranger count
- Zarr generation
