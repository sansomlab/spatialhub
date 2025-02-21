#!/usr/bin/env python3

import os
import csv
from pathlib import Path
import pandas as pd
import glob 
import re
import argparse
import ast

import numpy as np
from tifffile import imread, imwrite


### PARAMETERS ###

print("Parsing arguments")
parser = argparse.ArgumentParser()

parser.add_argument("--projDir", default="None", type=str,
                    help="path to the project directory containing raw data")
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

# Define paths to source and destination directories for corresponding sample
path2tiff = f"{args.projDir}/{slideName}/Morphology2D/"
outDir = f"ashlar.dir/{slideName}/Morphology2D/{sampleName}"
os.makedirs(outDir, exist_ok=True)

# Retrieve list of *all* raw TIFF files for corresponding *slide*
tif_files = glob.glob(os.path.join(path2tiff, "*.TIF"))
tif_basename = Path(os.path.basename(tif_files[0])).stem  # extract basename without file extension
    # to get corresponding file naming convention (minus up to 5 digits of FOV number)


### TASKS ###

## Generate mock FOV file if missing

print("Testing if mock FOV tile exists.")

if "mock_fov.ome.tiff" not in os.listdir("ashlar.dir"):
    print("Creating blank FOV tile.")
    # Read in any FOV to determine its dimensions
    exampleFOV = tif_files[0]
    fov = imread(exampleFOV)

    # Create matching array of zeros, with same dimensions and data type
    blank = np.zeros(fov.shape, fov.dtype)

    # Save as pyramidal .ome.tiff
    imwrite('ashlar.dir/mock_fov.ome.tiff', blank, bigtiff = True)
    print("Blank FOV tile created.")


## Symlink files into one folder per sample (group)

print("Creating symbolic links to source TIFF files.")

# Reformat 'fov_sequence' as iterable list
assert "fov_sequence" in df.columns, "Please specify the FOV sequence for stitching with Ashlar in the metadata (fov2sample) file"
x = df["fov_sequence"][0]
fov_pattern = x.split(',')

i = 0
for fov in fov_pattern:

    # Define destination symlink using iterator
    destination_directory = os.path.abspath(os.path.join(outDir, tif_basename[:-5] + f"{i:05}" + "R.TIF"))

    # Extract file name ending with "00{fov}.TIF" from source tiff directory, if applicable
    regex_pattern = f"00{fov}\\.TIF$" 
    matching_files = [file for file in tif_files if re.search(regex_pattern, os.path.basename(file))]

    if fov == "blank":
        os.symlink(os.path.abspath("ashlar.dir/mock_fov.ome.tiff"), destination_directory)
    
    elif matching_files:
        if len(matching_files) > 1:
           print(f"ERROR: Several matching TIFF files for FOV {fov}")
        else:
            source_file = os.path.abspath(matching_files[0])
            os.symlink(source_file, destination_directory)
    
    else:
        print(f"ERROR: File matching pattern '00{fov}.TIF' not found in {path2tiff}.")
    
    i += 1


print("Symlinks creation complete.")
