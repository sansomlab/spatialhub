#!/usr/bin/env python3

#import sopa  # Note: spatialdata incompatible with latest dask, hence 'Future Warning'

import os
import sys
import argparse
from spatialdata import SpatialData
from spatialdata import models
import skimage as ski
import pandas as pd

### PARAMETERS ###

print("Parsing arguments")
parser = argparse.ArgumentParser()

parser.add_argument("--stitchedDir", default="./ashlar.dir", type=str,
                    help="path to the project directory containing stitched data")
parser.add_argument("--sampleKey", default="None", type=str,
                    help="unique sample_id as in samples.tsv file, referring to the sample to execute the script on")
parser.add_argument("--fov2sample", default="samples.tsv", type=str,
                    help="path to the spatialhub samples TSV table")
parser.add_argument("--channels", default="PanCK,CD68,Membrane,CD45,DNA", type=str,
                    help="channel dictionary")
parser.add_argument("--rmProbes", default="Neg,Sys,Bac", type=str,
                    help="patterns of probes to exclude from main transcripts table")
parser.add_argument("--probeKey", default="target", type=str,
                    help="name of variable listing probes in transcripts file (e.g. target or gene)")
parser.add_argument("--x", default="x_global_px", type=str,
                    help="names of X spatial coordinate in transcript file")
parser.add_argument("--y", default="y_global_px", type=str,
                    help="names of Y spatial coordinate in transcript file")

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
print("Creating spatial data object for sample:", sampleName, slideName)


# Define paths to source and destination directories for corresponding sample
path2files = f"{args.stitchedDir}/{slideName}/Stitched2D/"
zarrDir = f"zarr.dir/{slideName}/"
os.makedirs(zarrDir, exist_ok=True)

outFile = os.path.join(zarrDir, sampleName + ".zarr")
if os.path.exists(outFile):
    print(outFile + " already exists - exiting .zarr creation script.")
    sys.exit()  # Terminate the script with status code 0


# Parse channel dictionary
x = str(args.channels)
x = x.split(",")
x = list(map(str, x))
channels_dict = x
print("Using the following channel dictionary:", channels_dict)

# Parse probes to exclude
x = str(args.rmProbes)
rm_probes_regex = x.replace(",", "|")
print("Using the following regex pattern for probes to exclude:", rm_probes_regex)


### TASKS ###

# Import using SOPA
#sdata = sopa.io.ome_tif(path2img, as_image=False)
# => here, we are losing the channel dictionary...

# Perform tissue-level segmentation, to remove outlier cells
#sopa.segmentation.tissue(sdata, mode = "staining", 
#                         level = 1,  # tissue-level segmentation can be performed on lower-quality image
#                         expand_radius_ratio = 0.05,  # expanding by 5% around defined shape
#                         drop_threshold = 0.01  # shapes covering >1% of the total area will be removed 
#                         )


# Import *stitched* image as array

print("Importing source image")
path2img = os.path.join(path2files, sampleName + "_stitched.ome.tiff")
img = ski.io.imread(path2img)
imgs = models.Image2DModel.parse(img, rgb = None, 
                                 c_coords = channels_dict)
    # this step converts the 'raw' image for compatibility with sdata format


# Import *stitched* transcripts file

print("Importing source transcripts file")
path2tx = os.path.join(path2files, sampleName + "_stitched_tx.csv")
tx = pd.read_csv(path2tx)

# split matrix between main panel and others (not to be considered in cell segmentation)
# => in case of bacterial probes on mouse panel, exclude these as well
idx0 = tx[args.probeKey].str.contains(rm_probes_regex)

tx0 = tx[idx0].copy()
tx0 = models.PointsModel.parse(tx0, feature_key = args.probeKey,
                               coordinates = {'x': args.x, 'y': args.y})

txs = tx[~idx0].copy()
txs = models.PointsModel.parse(txs, feature_key = args.probeKey,
                               coordinates = {'x': args.x, 'y': args.y})


# Create 'master' SpatialData object and write it to disk

print("Creating and writing SpatialData object")
sdata = SpatialData(images = {'image': imgs},
                    points = {'transcripts': txs,
                              'tx_other': tx0})
sdata.write(outFile, overwrite=True)
