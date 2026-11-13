# Pipelines

## cosmx_makeZarr
This workflow generates one Zarr file per CosMx sample from raw data. It supports two branches depending on whether Ashlar is used.

- With Ashlar: cosmx_genBlankFOV, cosmx_completeGrid, cosmx_runAshlar, cosmx_makeZarr.
- Without Ashlar: cosmx_assembleFOVs, cosmx_makeZarr.

## visiumhd_makeZarr
This workflow generates one Zarr file per Visium HD capture area from FASTQ files and images. It consists of two steps:

- Space Ranger count,
- Zarr generation.

## visiumhd_prepareBin2Cell
This workflow generates intermediate files up to the StarDist step of the Bin2Cell pipeline.

## runRCTD & runCell2Location
These two workflows deconvolve cell-type abundances in spatial transcriptomics data using a curated single-cell reference atlas, with either RCTD or Cell2Location.


## extractH5AD & aggregateH5AD

Both workflows generate an `AnnData` object stored in an `H5AD` file from a `SpatialData` object stored in a `Zarr` file.

They differ in whether a ready-to-use table is already present in the SpatialData object:

- `extractH5AD` is used when the table already exists, which is typically the case for Visium HD data.
- `aggregateH5AD` is used when the table is absent. It constructs an AnnData object by aggregating information from the points and shapes elements of the SpatialData object, which is typically the case for CosMx data.