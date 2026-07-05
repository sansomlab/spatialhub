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

## Usage

SpatialHub provides the following command-line interface (CLI). A minimal example is shown below:

```bash
spatialhub $module_name full
```

The complete command-line interface is shown below:

```text
usage: spatialhub [-h] [--dry] [--cores CORES] [--jobs JOBS] [--lock] [-v] module task

An integrated platform for spatial transcriptomics analysis.

positional arguments:
  module         The module to run.
  task           The task to run, [config|full|<rulename>].

options:
  -h, --help     show this help message and exit
  --dry          Perform a dry run only.
  --cores CORES  Number of cores to use.
  --jobs JOBS    Number of parallel jobs to run.
  --lock         Lock output by chmod.
  -v, --version  show program's version number and exit
```

A special feature is provided via the `--lock` option. When enabled, SpatialHub applies `chmod -R a-w` to the output directory to prevent accidental modification or deletion of the results.

To unlock the output directory for modification or deletion, run:

```bash
chmod -R +w /folder/to/delete
```

## Workflows

### CosMx makeZarr

This workflow generates one Zarr file per CosMx sample from raw data. It supports two branches depending on whether Ashlar is used.

- With Ashlar: cosmx_genBlankFOV, cosmx_completeGrid, cosmx_runAshlar, cosmx_makeZarr.
- Without Ashlar: cosmx_assembleFOVs, cosmx_makeZarr.

### Visium HD makeZarr

This workflow generates one Zarr file per Visium HD capture area from FASTQ files and images. It consists of two steps:

- Space Ranger count,
- Zarr generation.
