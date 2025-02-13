#!/usr/bin/env python3

import os
import csv
import shutil
from pathlib import Path
import pandas as pd
import glob 
import re
import argparse

import numpy as np
from tifffile import imread, imwrite


### PARAMETERS ###

print("Parsing arguments")
parser = argparse.ArgumentParser()

parser.add_argument("--projDir", default="None", type=str,
                    help="path to the project directory")
parser.add_argument("--slideName", default="None", type=str,
                    help="`basename` of the CosMx slide used for directory and file naming")
parser.add_argument("--fov2sample", default="samples.tsv", type=str,
                    help="path to the spatialhub samples TSV table")

args = parser.parse_args()


### SETUP ###

path2tiff = f"{args.projDir}/data/raw/{args.slideName}/Morphology2D/"
path2meta = f"{args.fov2sample}"
outDir = f"{args.projDir}/data/grouped/{args.slideName}/"

os.makedirs(outDir, exist_ok=True)
os.makedirs(os.path.join(outDir, "Morphology2D"), exist_ok=True)

# Define path to each TIFF file
tif_files = glob.glob(os.path.join(path2tiff, "*.TIF"))

# Read in metadata file
df = pd.read_csv(path2meta, sep = "\t")
if "sample_name" not in df.columns:
    print("'sample_name' not found: setting it to 'sample_id'")
    df["sample_name"] = df["sample_id"]

# Subset to slide of interest
df = df[df["slide_id"] == args.slideName]


### TASKS ###

## Generate mock FOV file if missing

print("Testing if mock FOV tile exists.")

#assert "mock_fov.ome.tiff" in os.listdir(), "Missing mock FOV file in current working directory"
if "mock_fov.ome.tiff" not in os.listdir():
    print("Creating blank FOV tile.")
    # Read in any FOV to determine its dimensions
    exampleFOV = tif_files[0]
    fov = imread(exampleFOV)

    # Create matching array of zeros, with same dimensions and data type
    blank = np.zeros(fov.shape, fov.dtype)

    # Save as pyramidal .ome.tiff
    imwrite('mock_fov.ome.tiff', blank, bigtiff = True)
    print("Blank FOV tile created.")


## Read in FOV pattern from metadata, and add mock FOVs/rename files accordingly

print("Renaming source TIFF files for Ashlar.")
assert "fov_sequence" in df.columns, "Please specify the FOV sequence for stitching with Ashlar in the metadata (fov2sample) file"

f0 = Path(os.path.basename(tif_files[0])).stem  # extract basename without file extension
    # to get corresponding file naming convention (minus up to 5 digits of FOV number)

for sample in list(set(df["sample_name"])):
    print(sample)
    fov_pattern = list(set(df["fov_sequence"][df["sample_name"] == sample]))[0]
    fov_pattern = fov_pattern.split(",")
    
    # Set tif_files object to the list of TIFF files for the corresponding sample
    sample_directory = os.path.join(outDir, "Morphology2D", sample)
    tif_files = os.listdir(sample_directory)

    for i in range(0, len(fov_pattern)):
        print(i, fov_pattern[i])
        matching_files = [file for file in tif_files if re.search(fov_pattern[i], os.path.basename(file))]
        if len(matching_files) == 0:
            shutil.copy("mock_fov.ome.tiff", sample_directory)
            os.rename(os.path.join(sample_directory, "mock_fov.ome.tiff"),
                      os.path.join(sample_directory, f0[:-5] + f"{i:05}" + "R.TIF"))
        else:
            os.rename(os.path.join(sample_directory, f0[:-5] + f"{int(fov_pattern[i]):05}" + ".TIF"),
                      os.path.join(sample_directory, f0[:-5] + f"{i:05}" + "R.TIF"))
            # adding "R" after FOV number, to make it clear those are *R*enamed FOVs

print("FOV file renaming operation complete.")
