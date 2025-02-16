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
print(df)
slideName = list(set(df["slide_id"]))[0]
sampleName = list(set(df["sample_name"]))[0]

# Define path to data directories for corresponding sample
path2tiff = f"{args.projDir}/data/grouped/{slideName}/Morphology2D/{sampleName}"
assert os.path.exists(path2tiff), "Sample TIFF directory does not exist. Please run Ashlar setup first."

# Define path to source TIFF files
tif_files = glob.glob(os.path.join(path2tiff, "*.TIF"))


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

print("Working on", sampleName)
    
fov_pattern = list(set(df["fov_sequence"][df["sample_name"] == sampleName]))[0]
fov_pattern = fov_pattern.split(",")
    
# Set tif_files object to the list of TIFF files for the corresponding sample
tif_files = os.listdir(path2tiff)


# Check whether renamed TIFF files already exist in the directory
regex_pattern = f"R.TIF$" 
matching_files = [file for file in tif_files if re.search(regex_pattern, os.path.basename(file))]

if matching_files:
    print(f"WARNING: Renamed TIFF files are already present in the sample-specific directory. Skipping renaming operation.")

else:
    for i in range(0, len(fov_pattern)):
        print(i, fov_pattern[i])
        matching_files = [file for file in tif_files if re.search(fov_pattern[i], os.path.basename(file))]
        if len(matching_files) == 0:
            shutil.copy("mock_fov.ome.tiff", path2tiff)
            os.rename(os.path.join(path2tiff, "mock_fov.ome.tiff"),
                      os.path.join(path2tiff, f0[:-5] + f"{i:05}" + "R.TIF"))
        else:
            os.rename(os.path.join(path2tiff, f0[:-5] + f"{int(fov_pattern[i]):05}" + ".TIF"),
                      os.path.join(path2tiff, f0[:-5] + f"{i:05}" + "R.TIF"))
            # adding "R" after FOV number, to make it clear those are *R*enamed FOVs


print("FOV file renaming operation complete.")
