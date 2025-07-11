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

#parser.add_argument("--celltypistModel", default="None", type=str,
#                    help="name of existing celltypist model to use")
parser.add_argument("--majorityVoting", default=False, type=bool,
                    help="name of existing celltypist model to use")
parser.add_argument("--ct_mode", default='best match', type=str,
                    help="see celltypist Read The Docs")
parser.add_argument("--ct_min_prop", default='best match', type=str,
                    help="see celltypist Read The Docs")
parser.add_argument("--atlasKey", default="None", type=str,
                    help="key of reference atlas to use for cell typing")
parser.add_argument("--atlasTSV", default="None", type=str,
                    help="path to the spatialhub TSV table describing reference dataset(s) to use for cell annotation")
parser.add_argument("--queryKey", type=str,
                    help="")

args = parser.parse_args()


### SETUP ###

# Read in metadata file
path2meta = f"{args.atlasTSV}"
df0 = pd.read_csv(path2meta, sep = "\t")

# Extract query dataset information
df_query = df0[df0["atlas_id"] == args.queryKey]
assert df_query.shape[0] == 1, "atlas data table can only list one single query (spatial) dataset"
df_query.index = [0]
query_key = df_query['atlas_id'][0]
print("Working with the following query (spatial) dataset to annotate: ")
print(df_query)

# Subset to atlas of interest
df = df0[df0["atlas_id"] == args.atlasKey]
assert df.shape[0] == 1, "atlas data frame can only have one single row"
df.index = [0]
atlas_key = df["atlas_id"][0]
print("Working with the following reference atlas for annotation: ")
print(df)

# Define path to pre-processed query dataset
queryPath = df_query['path'][0]

# Create out directories
model_name = atlas_key + '_' + df['celltype_annot_key'] + '_model'
outDir = os.path.join("annot.dir/celltypist", atlas_key)
print("Predictions output will be saved in:", outDir)
if not os.path.exists(outDir):
    os.mkdir(outDir)

# Celltypist options
models.models_path = outDir


### TASKS ###

sdata = sc.read_h5ad(queryPath)
print("Working from the following AnnData object:")
print(sdata)

if 'log1p' in sdata.layers:
    print("using existing 'log1p' layer")
    sdata.X = sdata.layers['log1p']
else:
    print("calculating normalized counts")
    sc.pp.normalize_total(sdata, target_sum=1e4)
    sc.pp.log1p(sdata)


### TASKS ###

# Load trained model

model = models.Model.load(model = os.path.join(outDir, model_name + '.pkl'))
print(model)


# Predict cell types

predictions = celltypist.annotate(sdata, 
                                  model = model_name, 
                                  mode = args.ct_mode,
                                  majority_voting = args.majorityVoting,
                                  min_prop = args.ct_min_prop)


# Save output

predictions.to_table(outDir, prefix = atlas_key + '_')
predictions.to_plots(outDir, prefix = atlas_key + '_')

if majorityVoting:
    sdata = predictions.to_adata(insert_conf_by = 'majority_voting', prefix = atlas_key + '_')
    df_scores = sdata.obs[[atlas_key + '_conf_score']]
    df_counts = predictions.summary_frequency(by = 'majority_voting')
else:
    sdata = predictions.to_adata(insert_conf_by = 'predicted_labels', prefix = atlas_key + '_')
    df_scores = sdata.obs[[atlas_key + '_conf_score']]
    df_counts = predictions.summary_frequency(by = 'predicted_labels')

scoresFile = os.path.join(outDir, atlas_key + '_confidence_scores.csv')
df_scores.to_csv(scoresFile)

countsFile = os.path.join(outDir, atlas_key + '_summary_counts.csv')
df_counts.to_csv(countsFile)
