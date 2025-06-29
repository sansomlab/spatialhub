# Miscellaneous tips

This documents points to more resources and lessons learned from working with `CosMx` data and developing `spatialhub`.


## Converting between python formats

### Generating `SpatialData` objects

- Main `spatialhub` script covering how to generate a `SpatialData` object: [`ashlar_zarr.py`]()
- Further examples of how to convert different output masks from different segmentation tools into a `GeoPandas` object compatible with `SpatialData` are covered in the [segmentation tutorial]()


### From `SpatialData` to `AnnData`

Once you have a `SpatialData` object including minimally a transcripts coordinates table and a segmentation mask, there are 2 main approaches to generate a cell (segment) x transcripts counts data table (further covered in [this example Jupyter notebook]()):

1. Using the [`SpatialData`]() library

2. Using the [`SOPA`]() library

> [!IMPORTANT]
> Make sure the index of your segmentation mask is set to what you want it to be, so that if you're importing some external metadata (e.g. from a cell type annotation tool), cells in the `AnnData` object have the correct index.


### From `AnnData` to `SpatialData`


### From `AnnData` to `spatial AnnData` (legacy format)
