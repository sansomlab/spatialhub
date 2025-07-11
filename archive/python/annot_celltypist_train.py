#!/usr/bin/env python3

import os
import argparse
import scanpy as sc
import pandas as pd
import celltypist
from celltypist import models


### PARAMETERS ###

print("Parsing arguments")
parser = argparse.ArgumentParser()

parser.add_argument("--atlasKey", default="None", type=str,
                    help="key of reference atlas to use for cell typing")
parser.add_argument("--atlasTSV", default="None", type=str,
                    help="path to the spatialhub TSV table describing reference dataset(s) to use for cell annotation")

args = parser.parse_args()


### SETUP ###

# Read in metadata file
path2meta = f"{args.atlasTSV}"
df0 = pd.read_csv(path2meta, sep = "\t")

# Subset to sample of interest
df = df0[df0["atlas_id"] == args.atlasKey]
assert df.shape[0] == 1, "atlas data frame can only have one single row"
df.index = [0]
atlas_key = df["atlas_id"][0]
print(df)

# Create out directory
outDir = os.path.join("annot.dir/celltypist", atlas_key)
print("Trained model will be saved in:", outDir)
if not os.path.exists(outDir):
    os.mkdir(outDir)


### TASKS ###

# Load reference AnnData (as subsetted for HVGs and/or probes in preliminary steps)
# WARNING: while Celltypist can work from a reference where only log1p data were shared
#          CosMx probe level aggregation should be handled carefully for such normalized data!
atlasPath = os.path.join("annot.dir/atlas.dir/", atlas_key + "_feature-subset.h5ad")
adata = sc.read_h5ad(atlasPath)

if 'log1p' in adata.layers.keys():
    print("using existing 'log1p' layer")
    adata.X = adata.layers['log1p']
else:
    print("calculating normalized counts")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

# Train model
model = celltypist.train(adata, labels = df['celltype_annot_key'], 
                         n_jobs = -1, feature_selection = True)

# Save!
model_name = atlas_key + '_' + df['celltype_annot_key'] + '_model'
print("Saving", model_name)
model.write(os.path.join(outDir, model_name + '.pkl'))
