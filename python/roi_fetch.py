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
                    help="Segmentation mask for which to retrieve cells in the desired ROI")
parser.add_argument("--ROI", default='roi', type=str,
                    help="Name of the ROI in the SpatialData (.zarr) object")

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
print("Filtering baysor segmentation mask for:", sampleName, slideName)


# Define paths to source and destination directories for corresponding sample
path2zarr = f"zarr.dir/{slideName}/{sampleName}.zarr"
outDir = f"roi.dir/{args.segMask}/{slideName}/"
if not os.path.exists(outDir):
    os.mkdir(outDir)


### TASKS ###

# Import SpatialData object to fetch cells from
print("Importing SpatialData object")
sdata = sd.read_zarr(path2zarr)
print(sdata)

assert args.segMask in sdata.shapes.keys(), "SpatialData object does not include the requested segmentation mask"
assert args.ROI in sdata.shapes.keys(), "SpatialData object does not include the requested region of interest"

# subset to ROI
print("Subsetting SpatialData to its ROI")
roi_shape = sdata[args.ROI].geometry.iloc[0]
sdata_roi = polygon_query(
    sdata,
    polygon=roi_shape,
    target_coordinate_system="global",
)
print(sdata_roi)

# retrieve matching keys
print("Retrieving matching cell keys")
df = sdata_roi[args.segMask]
df = df.drop(columns = 'geometry')

# Adding top-level metadata variables, to identify the source sample and data
df['segmentation_mask'] = args.segMask
df['slide_id'] = slideName
df['sample_name'] = sampleName
print(df)

df.to_csv(os.path.join(outDir, sampleName + '_' + args.ROI + '_barcodes.csv'), index = False)
print("Done fetching cells in ROI.")
