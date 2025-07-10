#!/usr/bin/env python3

#import sopa  # Note: spatialdata incompatible with latest dask, hence 'Future Warning'

import os
import sys
import shutil
import argparse
import muspan as ms
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc


############################## PARAMETERS ##############################

print("Parsing arguments")
parser = argparse.ArgumentParser()

parser.add_argument("--sampleKey", default="None", type=str,
                    help="unique sample_id as in samples.tsv file, referring to the sample to execute the script on")
parser.add_argument("--fov2sample", default="samples.tsv", type=str,
                    help="path to the spatialhub samples TSV table")
parser.add_argument("--segMask", default='atomx', type=str,
                    help="Segmentation mask from which the cell types and centroids are derived")
parser.add_argument("--annotFile", default='None', type=str,
                    help="Path to the metadata file containing cell type annotations and spatial coordinates")
parser.add_argument("--celltypeKey", default='celltype_level1', type=str,
                    help="Name of the variable in metadata storing cell types to calculate spatial metrics for")
parser.add_argument("--minPopSize", default=30, type=int,
                    help="Minimum sample size of the cell type population for it to be considered in analyses")
parser.add_argument("--maxR", default=1000, type=int,
                    help="Maximum size of the radius to consider")
parser.add_argument("--step", default=10, type=int,
                    help="Annulus step")
parser.add_argument("--calcRipK", default=True, type=bool,
                    help="Whether to also calculate the Ripley's K statistic (correlated with croos-PCF)")

args = parser.parse_args()


############################## SETUP ##############################

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
sampleIndex = sampleName + '_' + slideName
print("Filtering segmentation mask for:", sampleName, slideName)


# Define paths to destination directories for corresponding sample
outPCF = f"muspan.dir/{args.segMask}/PCF/"
os.makedirs(outPCF, exist_ok=True)
outRipK = f"muspan.dir/{args.segMask}/RipK/"
os.makedirs(outRipK, exist_ok=True)


############################## TASKS ##############################

### Task 1: Generate MuSpan domain ###

# Import SpatialData metadata to calculate PCF from (could also load an AnnData and fetch its adata.obs slot)
print("Importing SpatialData metadata")
df0 = pd.read_csv(args.annotFile, index_col = 0)
print(df0)

if 'sample_index' not in df0.columns:
    assert 'sample_name' in df0.columns, "Metadata must include a unique 'sample_index' identifier column, or both 'sample_name' and 'slide_id'"
    assert 'slide_id' in df0.columns, "Metadata must include a unique 'sample_index' identifier column, or both 'sample_name' and 'slide_id'"
    df0['sample_index'] = df0['sample_name'] + '_' + df0['slide_id']

assert args.celltypeKey in df0.columns, "Cell type key not found in metadata"
assert 'CenterX_sample_px' in df0.columns, "Cell centroid coordinates for stitched samples not found in metadata"
assert 'CenterY_sample_px' in df0.columns, "Cell centroid coordinates for stitched samples not found in metadata"


# Filter metadata table for sample to create MuSpAn domain for

print("Creating MuSpAn domain for", sampleIndex)
dfx = df0[df0['sample_index'] == sampleIndex]
print(dfx)

mux = ms.domain(sampleIndex)

points = dfx[['CenterX_sample_px', 'CenterY_sample_px']].to_numpy(copy=True)
print(points)
mux.add_points(points, "cell_coordinates")

labels = dfx[args.celltypeKey].to_numpy(copy=True)
mux.add_labels("celltype", labels)

mux.estimate_boundary(method = 'convex hull')

print(mux)


### Task 2: Filter celltypes for population with at least N cells ###

all_celltypes = list(set(df0[args.celltypeKey]))
print("The global metadata file includes the following cell type annotations:")
print(all_celltypes)

pass_minPopSize = dfx[args.celltypeKey].value_counts() > args.minPopSize
pass_minPopSize = pass_minPopSize[pass_minPopSize == True]
celltypes_pass = [x for x in all_celltypes if x in pass_minPopSize.index.to_list()]
celltypes_pass = sorted(celltypes_pass)
print("The sample under study icnludes at least", args.minPopSize, "cells of each of the following populations:")
print(celltypes_pass)


### Task 3: Calculate cross-PCF statistic for each pairwise comparison

print("Calculating cross-PCF statistic for every possible pairwise comparison")
# https://docs.muspan.co.uk/latest/_collections/spatial_analysis_methods/Spatial%20stats%20-%201%20-%20pcf.html

# Initiate data frame to store results
radius = np.arange(0, args.maxR + 1, args.step)
pcf_df = pd.DataFrame(radius)

for i, celltype_i in enumerate(celltypes_pass):
    pop_A = ms.query.query(mux, ('label', 'celltype'), 'is', celltype_i)

    for j, celltype_j in enumerate(celltypes_pass):
        pop_B = ms.query.query(mux, ('label', 'celltype'), 'is', celltype_j)

        if j<=i:
            r, PCF_ij = ms.spatial_statistics.cross_pair_correlation_function(
                domain = mux,
                population_A = pop_A, population_B = pop_B,
                max_R = args.maxR,
                annulus_step = args.step,
                annulus_width = 10*args.step,  # use a large width compared to step for a smoother curve
                visualise_output = False
            )

            # Store PCF results in data frame
            pcf_df[celltype_i + '_vs_' + celltype_j] = pd.DataFrame(PCF_ij)

# Save final results
print("Saving PCF results")
pcf_df.index = radius
print(pcf_df)
pcf_df.to_csv(os.path.join(outPCF, mux.name + '_PCF.csv'))


### Task 4: Calculate Ripley's K statistic for each pairwise comparison

if args.calcRipK:

    print("Calculating Ripley's K statistic for every possible pairwise comparison")
    # https://docs.muspan.co.uk/latest/_collections/spatial_analysis_methods/Spatial%20stats%20-%205%20-%20RipleysK.html

    # Initiate data frame to store results
    radius = np.arange(0, args.maxR + 1, args.step)
    ripK_df = pd.DataFrame(np.pi*radius**2)

    for i, celltype_i in enumerate(celltypes_pass):
        pop_A = ms.query.query(mux, ('label', 'celltype'), 'is', celltype_i)

        for j, celltype_j in enumerate(celltypes_pass):
            pop_B = ms.query.query(mux, ('label', 'celltype'), 'is', celltype_j)

            if j<=i:

                r, ripK_ij = ms.spatial_statistics.cross_k_function(
                    domain = mux,
                    population_A = pop_A, population_B = pop_B,
                    max_R = args.maxR,
                    visualise_output = False
                )

                # Store PCF results in data frame
                ripK_df[celltype_i + '_vs_' + celltype_j] = pd.DataFrame(ripK_ij)

    # Save final results
    print("Saving Ripley's K results")
    ripK_df.index = radius
    print(ripK_df)
    ripK_df.to_csv(os.path.join(outRipK, mux.name + '_RipK.csv'))


print("Done calculating pairwise spatial statistics for cell populations of size", args.minPopSize, "cells or more in sample", sampleIndex)