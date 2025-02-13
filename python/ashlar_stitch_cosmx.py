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


### PARAMETERS ###

print("Parsing arguments")
parser = argparse.ArgumentParser()

parser.add_argument("--projDir", default="None", type=str,
                    help="path to the project directory")
parser.add_argument("--slideName", default="None", type=str,
                    help="`basename` of the CosMx slide used for directory and file naming")

args = parser.parse_args()


### SETUP ###

path2tiff = f"{args.projDir}/data/grouped/{args.slideName}/Morphology2D/"
path2meta = f"{args.projDir}/metadata/{args.slideName}_fov2sample.csv"
outDir = f"{args.projDir}/data/grouped/{args.slideName}/stitched.dir"

os.makedirs(outDir, exist_ok=True)

# Read in metadata file
df = pd.read_csv(path2meta)
if "sample_name" not in df.columns:
    print("'sample_name' not found: setting it to 'sample_id'")
    df["sample_name"] = df["sample_id"]

# Loop over different samples (in pipeline)


### PARAMETERS ###

# Input:
projDir = os.getenv("projDir")    # "/users/sansom/tme871/work/kir032_morrell_cosmx_cancer/shared/cosmx/"
slideName = os.getenv("slideName")  # "HNSCC_run1"
sampleName = os.getenv("sampleName")  # "HPV0_OXF046_1"
path2tiff = f"{projDir}/data/grouped/{slideName}/Morphology2D/{sampleName}/"

# Image format [have this information in the fov2sample file?]
fov_width = 3
fov_height = 3

# Microscope parameters
px_size = 0.120280945  # CosMx_default
tile_ovlp = 0.01

# Output:
outDir = f"{projDir}/data/grouped/{slideName}/ashlar.dir"
keepChannels = [0,2,4]  # Channel(s) to keep in 3-channel mode


### SETUP ###

os.makedirs(outDir, exist_ok=True)


### TASKS ###

# Find pattern for series to stitch
fov_files = os.listdir(path2tiff)
assert len(fov_files) == fov_width * fov_height, "FOV grid incomplete. Please insert mock FOVs as placeholders to fill gaps."

f0 = Path(os.path.basename(fov_files[0])).stem  # extract basename without file extension
    # to get corresponding file naming convention (minus up to 5 digits of FOV number)
pattern = f0[:-5] + "{series}R.TIF"


# Perform image stitching
print(f"Stitching TIF tiles for image {path2tiff}/{pattern}")

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

# Save stitched image with 3 shannels for Cellpose GUI
mosaic_light = Mosaic(aligner, aligner.mosaic_shape, channels=keepChannels, verbose=True)
outFile = os.path.join(outDir, sampleName + "_stitched_3-channel.ome.tiff")
writer = PyramidWriter([mosaic_light], outFile, verbose=True)
writer.run()
print()
print("Wrote 3-channel .ome.tif")


# Write out CSV with corrected position and shift distance for each FOV.                                                                                                                                   
fields = [s for w, s in reader.metadata.all_series]
df = pd.DataFrame(
            np.hstack([aligner.positions, aligner.shifts]),
                    columns=["Position_Y", "Position_X", "Shift_Y", "Shift_X"],
                    index=pd.Series(fields, name="Field"),
                ).sort_index(
                    axis="columns",
                    # Discard roundoff error in the low bits.
                ).round(5)
outFile = os.path.join(outDir, sampleName + "_ashlar_fov_positions.csv")
df.to_csv(outFile)
print("Wrote fov_positions.csv")
