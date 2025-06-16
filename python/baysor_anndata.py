#!/usr/bin/env python3

#import sopa  # Note: spatialdata incompatible with latest dask, hence 'Future Warning'

import os
import sys
import shutil
import argparse
import spatialdata as sd
import pandas as pd
import anndata as ad
import scanpy as sc
import sopa

### PARAMETERS ###

print("Parsing arguments")
parser = argparse.ArgumentParser()

parser.add_argument("--sampleKey", default="None", type=str,
                    help="unique sample_id as in samples.tsv file, referring to the sample to execute the script on")
parser.add_argument("--fov2sample", default="samples.tsv", type=str,
                    help="path to the spatialhub samples TSV table")

args = parser.parse_args()


### SETUP ###

# Read in metadata file
path2meta = f"{args.fov2sample}"
df = pd.read_csv(path2meta, sep = "\t")

if "sample_name" not in df.columns:
    print("'sample_name' not found: setting it to 'sample_id'")
    df["sample_name"] = df["sample_id"]

# Subset to sample of interest
df = df[df["sample_id"] == args.sampleKey]
assert df.shape[0] == 1, "sample data frame can only have one single row"
df.index = [0]
print(df)
slideName = list(set(df["slide_id"]))[0]
sampleName = list(set(df["sample_name"]))[0]
print("Generating AnnData from baysor segmentation mask for:", sampleName, slideName)


# Define paths to source and destination directories for corresponding sample
path2seg = f"baysor.dir/{slideName}/{sampleName}_baysor/"
path2zarr = f"zarr.dir/{slideName}/{sampleName}.zarr"
#tempzarr = f"baysor.dir/{slideName}/{sampleName}_temp.zarr/"


### TASKS ###

# Import SpatialData object with Baysor segmentation mask
print("Sourcing SpatialData object with final Baysor segmentation mask")
sdata = sd.read_zarr(path2zarr)
print(sdata)

# Set baysor mask index to desired cell keys (stored as 'id' iif running the Baysor pipeline straightaway)
sdata['baysor'].index = sdata['baysor']['id']
sdata['baysor'].index.name = None
#sdata['baysor'].index

# Subsetting sdata object to one single transcripts table, and renaming it as expected by default in SOPA
sdata = sdata.subset(['image', 'baysor', 'tx_main'])
sdata['transcripts'] = sdata['tx_main']
del sdata['tx_main']
sdata

# Aggregate with SOPA, to enable aggregation by fluroescence intensity channel as well
print("Aggregating main transcripts table using Baysor segmentation mask")
sopa.aggregate(sdata, shapes_key = 'baysor', aggregate_channels = True)
print(sdata)


# Now, let's perform a little bit of cleaning up, to match AtoMx convention at best
print("Cleaning up AnnData object and its metadata")
adata = sdata['table']
df = adata.obs
print(df)

# Adding top-level metadata variables, to identify the source sample and data
df['baysor_index'] = sdata['baysor']['id']
df['segmentation_mask'] = 'baysor'
df['slide_id'] = slideName
df['sample_name'] = sampleName
df = df.drop(columns = ['slide', 'region', 'area'])
print(df)


# Retrieve cell coordinates and elongation from original Baysor metadata file
# (as SOPA aggregation triggers a WARNING that these metrics are likely to be inaccurate)
print("Retrieving cell coordinates and elongation from Baysor output")
baysor_stats = pd.read_csv(os.path.join(path2seg, 'segmentation_cell_stats.csv'))

# Filter for cells that passed QC
baysor_stats.index = baysor_stats['cell'].astype(str).astype('category')
baysor_stats = baysor_stats[baysor_stats.index.isin(df['baysor_index'])]
baysor_stats = baysor_stats.reindex(df['baysor_index'])
print(baysor_stats)

df.index = df['baysor_index']
if baysor_stats.index.equals(df.index):
    df['CenterX_sample_px'] = baysor_stats['x']
    df['CenterY_sample_px'] = baysor_stats['y']
    df['Area'] = baysor_stats['area']
    df['Elongation'] = baysor_stats['elongation']
    df['Density'] = baysor_stats['density']

df.index = df['cell_id']
print(df)


# Next, fetching relevant information from adata slots
print("Retrieving mean fluorescence intensities from .obsm slot")
if df.index.equals(adata.obs_names):
    # Adding mean fluorescence intensities to metadata
    intens = adata.obsm['intensities']
    intens = intens.add_prefix('Mean.')
    df = df.join(intens)

print(df)


# Update adata.obs with cleaned up df
if df.index.equals(adata.obs_names):
    adata.obs = df

del adata.obsm
del adata.uns


# Save AnnData object
print("Saving the following AnnData object:")
print(adata)
adata.write_h5ad(os.path.join(path2seg, 'anndata.h5ad'))

print("Done generating Baysor AnnData.")
