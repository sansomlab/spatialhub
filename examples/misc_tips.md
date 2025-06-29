# Miscellaneous tips

This documents points to more resources and lessons learned from working with `CosMx` data and developing `spatialhub`.


## Converting between python formats

### Generating `SpatialData` objects

- Main `spatialhub` script covering how to generate a `SpatialData` object: [`ashlar_zarr.py`](../python/ashlar_zarr.py)
- Further examples of how to convert different output masks from different segmentation tools into a `GeoPandas` object compatible with `SpatialData` are covered in the [segmentation tutorial](./workflow_cosmx_mouse_gut_part2.md)

> [!NOTE]
> Overwriting of `.zarr` stacks is currently (as of June 2025) not fully supported. Whenever adding on to a `SpatialData` object, we recommend saving a copy, and then moving the new object to the old location.


### From `SpatialData` to `AnnData`

Once you have a `SpatialData` object including minimally a transcripts coordinates table and a segmentation mask, there are 2 main approaches to generate a cell (segment) x transcripts counts data table (further covered in [this example Jupyter notebook](https://github.com/sansomlab/spatialhub/blob/main/notebooks/transcripts_spatial_plots.ipynb)):

1. Using the [`SpatialData`](https://spatialdata.scverse.org/en/stable/) library (which requires the least pre-processing of the object)
2. Using the [`SOPA`](https://gustaveroussy.github.io/sopa/) library (which enables to calculate mean fluorescence intensities per channel and per cell (segment) if a microscopy image is available)

> [!TIP]
> Before aggregating the transcripts table, make sure the index of your segmentation mask is set to what you want it to be, so that if you're importing some external metadata (e.g. from a cell type annotation tool), cells in the `AnnData` object have the correct index.

> [!IMPORTANT]
> After aggregating the transcripts table, make sure the `adata.uns['spatialdata_attrs']` (a python dictionary) points to the correct keys - namely:
>    - `region` for the name of the segmentation mask
>    - `region_key` for the name of the variable in `adata.obs` where the `region` name is stored
>    - `instance_key` for the name of the variable in `adata.obs` where the cell (segment) unique ID (index) is stored
> If this parameters are incorrectly set, some plotting functionalities (e.g. colour-coding a cell by expression level of a transcript) will fail.


### From `AnnData` to `SpatialData`

If you've independently generated a `AnnData` and `SpatialData` object(s) for your project, as is likely to be the case if using `spatialhub`, you can easily link the two as long as the cell (segment) indexes in your `adata.obs_names` match those used in your `SpatialData`'s segmentation mask, and the `adata.uns['spatialdata_attrs']` dictionary is specified correctly (see section above).

**Example script to be added**


### From `AnnData` to `spatial AnnData` (legacy format)

Finally, you may want to use the legacy spatial `AnnData` (`squidpy`) format, which was the typical `python` structure for spatial transcriptomics datasets before the scverse introduced the `SpatialData` format. Some tools such as [`Squidpy`](https://squidpy.readthedocs.io/en/stable/) or [`SPIN`](https://github.com/wanglab-broad/spin) still expect this structure.

To turn a standard `AnnData` into a `spatial AnnData` object, simply store the cell (segment) centroid coordinates as a `numpy` array in `adata.obsm['spatial']` (you could use any key other than 'spatial' to store this array, but most tools will by default look for centroid coordinates in the 'spatial' slot), like this:

```
# Assuming df is a Pandas dataframe (e.g. adata.obs) listed coordinates for all cells
coords = df[['CenterX_sample_px', 'CenterY_sample_px']].copy()
coords = coords.to_numpy()
adata.obsm['spatial'] = coords
```

If you don't already have the cell (segment) centroid coordinates, you can generate them from any `SpatialData` object using the [`spatialdata.get_centroids()`](https://spatialdata.scverse.org/en/stable/api/operations.html#spatialdata.get_centroids) function (but **beware of coordinate system compatibilities between tools**: for example, we've noticed the coordinates output by `baysor` are not the same as those re-calculated by `spatialdata` using `get_centroids()`, due to different coordinate reference systems and spatial projections used - [see this GeoPandas guideline on the matter](https://geopandas.org/en/stable/docs/user_guide/projections.html)).

