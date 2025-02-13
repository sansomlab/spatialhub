#!/usr/bin/env python3

import os
import csv
import shutil
from pathlib import Path
import pandas as pd
import glob 
import re
import argparse
import ast


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
#os.makedirs(os.path.join(outDir, "fov2sample.dir"), exist_ok=True)
os.makedirs(os.path.join(outDir, "Morphology2D"), exist_ok=True)

# Define path to each TIFF file
tif_files = glob.glob(os.path.join(path2tiff, "*.TIF"))


# Read in metadata file
df = pd.read_csv(path2meta, sep = "\t")
if "sample_name" not in df.columns:
    print("'sample_name' not found: setting it to 'sample_id'")
    df["sample_name"] = df["sample_id"]

# Reformat fov_sequence and "explode" dataframe
x = df["fov_sequence"].replace("blank", "", regex=True)
x = x.replace(",{2,}", ",", regex=True)
x = x.replace("^,", "", regex=True); x = x.replace(",$", "", regex=True)
x = x.replace("^", "[", regex=True); x = x.replace("$", "]", regex=True)
x = x.apply(ast.literal_eval)  # convert string to list type
fov = x.explode()

df_fov = pd.DataFrame({'fov': fov})
df_fov["key"] = df_fov.index
df["key"] = df.index
df = df_fov.merge(df, on='key')

# Subset to slide of interest
df = df[df["slide_id"] == args.slideName]


### TASKS ###

## Group by 'sample_name' and save each group as a separate CSV

#for sample_name, group in df.groupby("sample_name"):
#    # Get the list of FOVs for the current sample and create dataframe
#    fovs = group["fov"].tolist()
#    sample_df = pd.DataFrame({"sample_name": [sample_name], "fovs": [fovs]})
#    # Save dataframe to a CSV file
#    outFile = os.path.join(outDir, "fov2sample.dir", f"{sample_name}_fov2sample.csv")
#    sample_df.to_csv(outFile, index=False)


## Copy files into one folder per sample (group)

print("Copying source TIFF files.")

for index, row in df.iterrows():
    # Extract the file number and sample name
    tif_number = str(row['fov'])
    sample = str(row['sample_name'])

    # Build the destination path and create it
    destination_directory = os.path.join(outDir, "Morphology2D", sample)
    os.makedirs(destination_directory, exist_ok=True)

    # Extract file names ending with "00{tif_number}.TIF" in source tiff directory
    regex_pattern = f"00{tif_number}\\.TIF$" 
    matching_files = [file for file in tif_files if re.search(regex_pattern, os.path.basename(file))]
    if matching_files:
        if len(matching_files) > 1:
            print(f"ERROR: Several matching TIFF files for FOV {tif_number}")
        else:
            source_file = matching_files[0]
            shutil.copy(source_file, destination_directory)
    else:
        print(f"ERROR: File matching pattern '00{tif_number}.TIF' not found in {path2tiff}.")

print("FOV files copy operation complete.")
