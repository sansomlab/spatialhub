# CosMx mouse gut analysis workflow - Part 2

The following markdown provides a step-by-step workflow to reproduce the **re-segmentation of a CosMx study of the mouse gut in an IBD model** using the `spatialhub` pipeline. (See Part 3 for the next steps of the analysis, and Part 1 for the required preliminary steps.)

> [!IMPORTANT]
> At this stage (June 2025), this part is very much work in progress, as we haven't been able to generate a segmentation mask which is convincingly better than the default `AtoMx` one. This document thus provides an overview of what has been tested.

> [!TIP]
> Regardless of the segmenation tool used, we recommend first running a pilot `spatialhub` run on a subset of samples (this can be achieved by pointing to a subset `samples.tsv` table in the YAML file, instead of the `samples.tsv` table covering the full study dataset). **This piloting is essential since most segmentation tools rely on a lot of parameters, which will need to be optimized for each dataset.**  


* * *

## Study overview

The experiment covers . . .

### Data storage on BMRC

Raw and processsed data alongside key scripts to reproduce the analysis on the BMRC cluster in this directory: `/users/sansom/tme871/work/cosmx_mouse`

Here is an overview of where to find key files:

- Original raw data, as exported from `AtoMx SIP` through their user interface: `/users/sansom/tme871/work/cosmx_mouse/data/raw/atomx`

> [!NOTE]
> At this stage of the analysis (June 2025), we are confidently moving away from `AtoMx` default analyses tools, so these data files could be deleted (knowing they can always be downloaded again from our `AtoMx SIP`). I have kept them as an example of an `AtoMx` export, but the minimal relevant files are copied in the directories described in the bullet point below.

- **Raw data, re-organized for analysis with `spatialhub`**: `/users/sansom/tme871/work/cosmx_mouse/data/raw/${SLIDE_ID}`, with the following slide IDs available: NL4S3a, NL4S3b, NL5H1a, NL5H1b ('a' indicates a KIR run, 'b' indicates a CAMS run - see experimental design section above)
- **Metadata**, including `spatialhub` samples.tsv table (full and filtered for samples passing QC, see details below) and lists of CosMx probes (standard 'LBL-11176-05-Mouse-Universal-Cell-Characterization-Gene-List.txt' and custom add-on panel '1K_Oxford_Sansom_add-on_list.csv')

> [!WARNING]
> **In this study, some of the custom probes use special characters (e.g. *SiglecF (170)* instead of *Siglecf*) and/or non-conventional gene symbols (e.g. *Ly6G* instead of *Ly6g*; *Cd64* instead of *Fcgr1*). In addition, some probes in the standard panel are cross-reactive and match multiple genes. Whenever comparing a CosMx probe panel to reference scRNA-seq datasets/pathway libraries/etc., it is essential to ensure that gene names are properly encoded for the intended purpose.** This point is highlighted again where most relevant in the workflow description.

- **IBD_mouse_scRNA-seq_atlas**: Directory (and [GitHub repo](https://github.com/sansomlab/IBD_mouse_scRNA-seq_atlas)) collecting data and scripts used to generate a scRNA-seq atlas fit for the purpose of annotating this specific CosMx dataset.
- **spatialhub**: Main analysis folder, containing YAML, log out output files/directories from steps run with the `spatialhub` pipeline, as well as additional scripts for steps requiring more customization.
- Other directories contain pilot scripts (irrelevant to reproduce the final analysis), used to test different spatial transcriptomics tools before rolling them out to the dataset as a whole and/or adding them to the `spatialhub` suite of pipelines, if applicable

* * *


## Option 1: Using image-based tools

Image-based segmentation tools rely on the fluorescence microscopy image acquired on the `CosMx` device to define cell types. Note that any registered microscopy image could also be used (work in progress in the Hallou lab). The default `AtoMx` mask itself is derived from `cellpose`.


### SOPA-cellpose

For re-segmentation using image-based tools, **`spatialhub` currently builds on top of the [`SOPA`](https://gustaveroussy.github.io/sopa/) implementation of [`cellpose`](https://cellpose.readthedocs.io/en/latest/). We thus recommend reading the documentation for both tools for additional details.**

> [!NOTE]
> `SOPA` is well suited to enable fast re-segmentation, since it patches the image and parallelize the segmentation to run on each patch before (properly) stitching the segmentation mask pieces into one (unlike `AtoMx` tools, which just append pieces of the segmentation mask next to each other without resolving conflicts at the boundaries). At this stage, we have only implemented the `cellpose` tool from `SOPA` into `spatialhub`, but the `sopa_segment_img` pipeline could be expanded to cover more tools (please refer to `SOPA` Read The Docs for available image-based segmentation tools).

After running `spatialhub ashlar` with the `create_zarr` parameter set to True, you can readily run the `sopa_segment_img` pipeline within the same directory, since the input files will be fetched from the `zarr.dir` directory created by Ashlar, using the following commands:

```
spatialhub sopa_segment_img config
spatialhub sopa_segment_img make full -v5 -p20
```

The YAML file covers major parameters from the [`sopa.segmenation.cellpose` function](https://gustaveroussy.github.io/sopa/api/segmentation/#sopa.segmentation.cellpose). How those parameters impact the resulting segmentation mask are further described in the [`cellpose` Read The Docs](https://cellpose.readthedocs.io/en/latest/settings.html). Beyond running `cellpose`, `sopa.segmenation.cellpose` includes some helpful segmentation mask clean-up functionalities, such as the ability to remove cells under a certain area (again, see `SOPA` for details).

The output segmentation mask will be automatically added to the `zarr.dir` directory, using the `mask_key` specified in the YAML file.

> [!TIP]
> The best results with `cellpose` are likely to be achieved when running it based on a custom, re-trained model (see [`cellpose` tutorial here](https://cellpose.readthedocs.io/en/latest/gui.html#training-your-own-cellpose-model)). If using a pre-trained model, we've had most success with the `nuclei`, `cyto3` and `tissuenet_cp3` models (see [here for a full list of available `cellpose` pre-trained modes](https://cellpose.readthedocs.io/en/v3.1.1.1/models.html#full-built-in-models).

In the case of this specific study, the best results were obtained with the `cellpose nuclei` model, although we did not retain it as it misses a large area (cytoplasm) of manuy cells (notably, epithelial cells).

* * *


## Option 2: Using transcripts-based tools

Transcripts-based segmentation tools take advantage of transcripts information (on top of the microscopy image, if using an image-based mask as prior) to determine cell boundaries. Such tools notably include [`baysor`](https://github.com/kharchenkolab/Baysor) and [`proseg`](https://github.com/dcjones/proseg). 

> [!NOTE]
> 1. `spatialhub` currently includes a `sopa_segment_tx` pipeline, which builds on top of `SOPA` and its patching capability to perform segmentation using `baysor`. While this pipeline can work, it does not seem to interface correctly with `dask` in the backend, and thus requires a lot of memory, thereby defeating the purpose of using `SOPA-baysor` rather than `baysor` straightaway.
> 2. While we have successfully completed some `proseg` runs, we've consistently found `baysor` results of superior quality. At this stage, `spatialhub` thus only implements `baysor`. (See [here](https://github.com/sansomlab/spatial_tx/tree/main/M2_segmentation/proseg) for an example of how to run `proseg` on one sample from this CosMx mouse study outside of the `spatialhub` pipeline)

### `baysor`

After running `spatialhub ashlar`, you can readily run the `baysor` pipeline within the same directory, since the input files will be fetched from the `ashlar.dir` directory created by Ashlar, using the following commands:

```
spatialhub baysor config
spatialhub baysor make full -v5 -p20
```

In addition to a YAML file, which handles the `spatialhub` add-on functionalities (including, additional cell filtering criteria), you will also need a TOML file to specify the desired `baysor` parameters, such as [this TOML file](https://github.com/kharchenkolab/Baysor/blob/master/configs/example_config.toml). Please refer to the `baysor` documentation for details of what each of these parameters mean.

The output segmentation mask will be automatically added to the `zarr.dir` directory, using `baysor` as a key. (If a mask with key `baysor` already exists in the `SpatialData` object, this pre-existing mask will be renamed `baysor_old`.)

> [!TIP]
> If [running `baysor` with a prior](https://kharchenkolab.github.io/Baysor/dev/segmentation/#Using-a-prior-segmentation) (recommended), we've found that using a high-quality image-based segmentation mask with a high prior confidence score (>0.8) provides the best results. This differs from the default recommended setting, which assumes a low confidence score of 0.2 for the prior.

* * *

## Assessing segmenation masks

For a quick review and QC of a given segmentation mask (especially, as generated from running the `baysor` pipeline, see this [`jupyter` notebook](https://github.com/sansomlab/spatialhub/blob/main/notebooks/segmentation_review_baysor.ipynb).

* * *

## For developers: Misc lessons from piloting segmentation tools

Converting output segmentation masks to a compatible `SpatialData` shapes object can be tricky. If trying to solve this problem for other tools, see for example how it was solved for the `AtoMx` and `baysor` masks in the following scripts:

### a) `baysor_filter.py` lines 61-67 

Although the `GeoJSON` format should be straightforward to import as a `GeoPandas` dataframe, we've often run into issues (documented on GitHub for old versions of `baysor`), sovled by using the following filter for invalide geometries:

```
# Import JSON polygon file
print("Importing JSON segmentation file (removing invalid geometries)")
baysor_mask = gpd.read_file(os.path.join(path2seg, 'segmentation_polygons_2d.json'), on_invalid = 'warn')

# Filter out invalid geometries (there should be none in latest releases of baysor...)
baysor_mask = baysor_mask[baysor_mask['geometry'].is_valid]
print(baysor_mask)
```

### b) `ashlar_zarr.py` lines 156-227

This is an example of manually converting a segmentation mask provided as a flat data table to a `shapely` geometry object and eventually `SpatialData` object. Note that for samples/FOVs, we've found this script to fail due to invalid geometries... hence the `ashlar_zarr_DEBUG.py` script, which should **ONLY** be used when the default `ashlar_zarr.py` script failed, lest some relevant cells are lost from other samples/FOVs.

```
# import original segmentation mask for corresponding FOVs
print("Loading segmentation mask for relevant FOVs")

atomx_ls = []
for fov in fov_ls:
  # open file
  k = 5 - len(str(fov))
  file = os.path.join(path2seg, 'CellBoundaries_F' + '0'*k + str(fov) + '.csv')
  df_atomx = pd.read_csv(file)
  atomx_ls = atomx_ls + [df_atomx]

df_atomx = pd.concat(atomx_ls)
print(df_atomx)

# Convert coordinates from local to global in newly stitched image: 
print("Converting x,y shape points to new coordinates system")

# reformat as one list of paired coordinates per cell, as will be needed for later conversion to geopandas
df_atomx["atomx_index"] = "FOV" + df_atomx["fov"].astype(str) + "_C" + df_atomx["cellID"].astype(str)
cell_ids = list(set(df_atomx['atomx_index'].to_list()))

df_ls = [y for x, y in df_atomx.groupby("atomx_index")]
#df_ls[0]

coords_ls = []
for dfk in df_ls:
        
  # dfk = set of polygon coordinates for cell k
  dfk.index = range(0, len(dfk))  # re-index from zero to n-1

  # update coordinates from local to global, using Ashlar's FOV position file
  fov = dfk['fov'][0]
  f = (ashlar_df['FOV'] == fov)
  x = ashlar_df[f].iloc[0]['Position_X']; y = ashlar_df[f].iloc[0]['Position_Y']
  dfk['x_global'] = dfk['x_local'] + x; dfk['y_global'] = dfk['y_local'] + y

  # store as list for conversion to Shapely polygon
  coords_k = []
  for i in range(0, len(dfk)):
    coords_k = coords_k + [[dfk['x_global'][i], dfk['y_global'][i]]]
    coords_ls = coords_ls + [shapely.Polygon(coords_k)]

len(coords_ls)

# convert to geopandas format
print("Transforming to GeoPandas object for compatibility with SpatialData")
polygon_gdf = gpd.GeoDataFrame(geometry=coords_ls)
polygon_gdf["atomx_index"] = cell_ids
polygon_gdf

# convert to zarr compatible format
print("Adding shapes element to SpatialData")
dfs_atomx = models.ShapesModel.parse(polygon_gdf)
sdata['atomx'] = dfs_atomx
```

### c) other

`SOPA` also includes a tool to support conversion of different segmentation mask for compatibility with the `SpatialData` format - see for example [this script](https://github.com/sansomlab/spatial_tx/blob/main/M2_segmentation/cellpose/cellpose_conversion.py); however, they require a lot of resources!
