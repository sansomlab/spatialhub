#!/usr/bin/env python3

import os
from pathlib import Path
from ashlar.fileseries import FileSeriesReader
from ashlar.reg import EdgeAligner, Mosaic, PyramidWriter
import numpy as np
import pandas as pd
import skimage.exposure, skimage.filters, skimage.io
import sys
import tifffile
import argparse
import ast
import glob


### PARAMETERS ###

print("Parsing arguments")
parser = argparse.ArgumentParser()

parser.add_argument("--sampleKey", default="None", type=str,
                    help="unique sample_id as in samples.tsv file, referring to the sample to execute the script on")
parser.add_argument("--fov2sample", default="samples.tsv", type=str,
                    help="path to the spatialhub samples TSV table")
parser.add_argument("--pxSize", default=0.120280945, type=str,
                    help="pixel to micrometer conversion factor")
parser.add_argument("--tileOvlp", default=0.01, type=str,
                    help="overlap between adjacent FOV tiles")
parser.add_argument("--keepChannels", default="0,2,4", type=str,
                    help="channels to keep in 3-channel saving mode")

args = parser.parse_args()

# CosMx microscope parameters
px_size = float(args.pxSize)
tile_ovlp = float(args.tileOvlp)

x = args.keepChannels
x = x.split(",")
x = list(map(int, x))
keepChannels = x
print(keepChannels)


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

# Extract FOV width and height
fov_width = list(set(df["fov_width"]))[0]
fov_height = list(set(df["fov_height"]))[0]

# Define path to data directories for corresponding sample
path2tiff = f"ashlar.dir/{slideName}/Morphology2D/{sampleName}"
assert os.path.exists(path2tiff), "Sample TIFF directory does not exist. Please run Ashlar setup first."

# Define path to source TIFF files
tif_files = glob.glob(os.path.join(path2tiff, "*.TIF"))

# Define and make directory to output stitched images
outDir = f"ashlar.dir/{slideName}/Stitched2D"
os.makedirs(outDir, exist_ok=True)


### TASKS ###

# Find pattern for series to stitch
fov_files = os.listdir(path2tiff)
assert len(fov_files) == fov_width * fov_height, "FOV grid incomplete. Please insert mock FOVs as placeholders to fill gaps."

f0 = Path(os.path.basename(fov_files[0])).stem  # extract basename without file extension
    # to get corresponding file naming convention (minus up to 5 digits of FOV number)
pattern = f0[:-5] + "{series}R.TIF"
print(pattern)


### Perform image stitching

print(f"Stitching TIFF tiles for image {path2tiff}")

reader = FileSeriesReader(
            path2tiff,
            pattern=pattern,
            width=fov_width,
            height=fov_height,
            overlap=tile_ovlp,
            pixel_size=px_size,
           )

aligner = EdgeAligner(reader, channel=4, verbose=True)
aligner.run()

# Save stitched image with all shannels
mosaic = Mosaic(aligner, aligner.mosaic_shape, channels=[0,1,2,3,4], verbose=True)
outFile = os.path.join(outDir, sampleName + "_stitched.ome.tiff")
writer = PyramidWriter([mosaic], outFile, verbose=True)
writer.run()
print()
print("Wrote full .ome.tif")

# Save stitched image with 3 shannels (for Cellpose GUI)
mosaic_light = Mosaic(aligner, aligner.mosaic_shape, channels=keepChannels, verbose=True)
outFile = os.path.join(outDir, sampleName + "_stitched_3-channel.ome.tiff")
writer = PyramidWriter([mosaic_light], outFile, verbose=True)
writer.run()
print()
print("Wrote 3-channel .ome.tif")


### Create table of corrected positions and shifts for each FOV

print(f"Writing corrected FOV positions for image {path2tiff}")

fields = [s for w, s in reader.metadata.all_series]
coords_df = pd.DataFrame(
        np.hstack([aligner.positions, aligner.shifts]),
        columns=["Position_Y", "Position_X", "Shift_Y", "Shift_X"],
        index=pd.Series(fields, name="Field"),
    ).sort_index(
        axis="columns"
    ).round(5)
        # the latter discards roundoff error in the low bits.

# Add back original FOV indexes to coordinates data frame
x = df["fov_sequence"][0]
x = x.split(",")  # convert string to list type
assert len(x) == coords_df.shape[0], "missing fields in Ashlar coordinates table"
coords_df["FOV"] = x

outFile = os.path.join(outDir, sampleName + "_ashlar_fov_positions.csv")
coords_df.to_csv(outFile)
print("Wrote corrected FOV positions file")
