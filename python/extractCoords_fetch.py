#!/usr/bin/env python3

#import sopa  # Note: spatialdata incompatible with latest dask, hence 'Future Warning'

import os
import sys
import shutil
import argparse
import spatialdata as sd
from spatialdata import polygon_query
import pandas as pd

### PARAMETERS ###

print("Parsing arguments")
parser = argparse.ArgumentParser()

parser.add_argument("--sampleKey", default="None", type=str,
                    help="unique sample_id as in samples.tsv file, referring to the sample to execute the script on")
parser.add_argument("--fov2sample", default="samples.tsv", type=str,
                    help="path to the spatialhub samples TSV table")
parser.add_argument("--segMask", default='atomx', type=str,
                    help="Segmentation mask for which to retrieve centroid coordinates")
parser.add_argument("--segIndex", default='atomx_index', type=str,
                    help="Name of index key in segmentation GeoPandas dataframe storing unique segment ID (more explicit than 1:N)")

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
print("Extracting centroid coordinates for:", sampleName, slideName)


# Define paths to source and destination directories for corresponding sample
path2zarr = f"zarr.dir/{slideName}/{sampleName}.zarr"
outDir = f"coordinates.dir/{args.segMask}/{slideName}/"
os.makedirs(outDir, exist_ok=True)


### TASKS ###

# Import SpatialData object to fetch coordinates from
print("Importing SpatialData object")
sdata = sd.read_zarr(path2zarr)
print(sdata)

assert args.segMask in sdata.shapes.keys(), "SpatialData object does not include the requested segmentation mask"
assert args.segIndex in sdata[args.segMask].columns, "Segment (cell) index variable not found in GeoPandas data frame"

# Set atomx mask index to desired cell keys, to avoid indexing issues later
sdata[args.segMask].index = sdata[args.segMask][args.segIndex]
sdata[args.segMask].index.name = None
sdata[args.segMask].index

# Calculating centroids
print("Calculating cell centroid coordinates")
coords = sd.get_centroids(sdata[args.segMask])
coords = coords.compute()
coords = coords.rename(columns={"x": "CenterX_sample_px", "y": "CenterY_sample_px"})

# Adding top-level metadata variables, to identify the source sample and data
coords['segmentation_mask'] = args.segMask
coords['slide_id'] = slideName
coords['sample_name'] = sampleName
print(coords)

coords.to_csv(os.path.join(outDir, sampleName + '_sampleCoords.csv'), index = True)
print("Done retrieving centroid coordinates.")
