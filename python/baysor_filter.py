#!/usr/bin/env python3

#import sopa  # Note: spatialdata incompatible with latest dask, hence 'Future Warning'

import os
import sys
import shutil
import argparse
import spatialdata as sd
from spatialdata import models
import pandas as pd
import shapely
import geopandas as gpd


### PARAMETERS ###

print("Parsing arguments")
parser = argparse.ArgumentParser()

parser.add_argument("--sampleKey", default="None", type=str,
                    help="unique sample_id as in samples.tsv file, referring to the sample to execute the script on")
parser.add_argument("--fov2sample", default="samples.tsv", type=str,
                    help="path to the spatialhub samples TSV table")
parser.add_argument("--minArea", default=500, type=int,
                    help="minimal area for a polygon to be considered as a valid cell")
parser.add_argument("--minCount", default=20, type=int,
                    help="minimal total transcript count for a polygon to be considered as a valid cell")

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
path2seg = f"baysor.dir/{slideName}/{sampleName}_baysor/"
path2zarr = f"zarr.dir/{slideName}/{sampleName}.zarr"
tempzarr = f"baysor.dir/{slideName}/{sampleName}_temp.zarr/"


### TASKS ###

# Import JSON polygon file
print("Importing JSON segmentation file (removing invalid geometries)")
baysor_mask = gpd.read_file(os.path.join(path2seg, 'segmentation_polygons_2d.json'), on_invalid = 'warn')

# Filter out invalid geometries (there should be none in latest releases of baysor...)
baysor_mask = baysor_mask[baysor_mask['geometry'].is_valid]
print(baysor_mask)


# Import cell statistics to filter those that are more likely to be artifacts
print("Importing polygon statistics file to filter out artifacts")
baysor_stats = pd.read_csv(os.path.join(path2seg, 'segmentation_cell_stats.csv'))
baysor_stats = baysor_stats[baysor_stats['cell'].isin(baysor_mask['id'])]
baysor_stats = baysor_stats[baysor_stats['area'] > args.minArea]
baysor_stats = baysor_stats[baysor_stats['n_transcripts'] >= args.minCount]
print(baysor_stats)

# Update baysor_mask accordingly
baysor_mask = baysor_mask[baysor_mask['id'].isin(baysor_stats['cell'])]
baysor_mask.to_file(os.path.join(path2seg, 'segmentation_polygons_filtered.geojson'), driver='GeoJSON')


# Import SpatialData object to add segmentation mask to
print("Appending baysor mask to SpatialData object")
sdata = sd.read_zarr(path2zarr)
baysor_mask = models.ShapesModel.parse(baysor_mask)
sdata['baysor'] = baysor_mask
print(sdata)


# Save mask to .zarr store
print("Saving updated SpatialData object")
sdata.write(file_path=tempzarr, overwrite=True)
    # unfortunately, overwriting does not work with spatialdata
    # problem documented here: https://github.com/scverse/spatialdata/discussions/520

# Move the relevant 'shapes' contents to 'master' sdata store
source = os.path.join(tempzarr, 'shapes', 'baysor')
destin = os.path.join(path2zarr, 'shapes')
shutil.move(source, destin)
shutil.rmtree(tempzarr)

print("Done cleaning up Baysor output.")