#!/usr/bin/env python3

import os
import argparse
import scanpy as sc
import anndata as ad
import pandas as pd
from scipy.sparse import csr_matrix


############################## PARAMETERS ##############################

print("Parsing arguments")
parser = argparse.ArgumentParser()

parser.add_argument("--atlasKey", default="None", type=str,
                    help="key of reference atlas to use for cell typing")
parser.add_argument("--atlasTSV", default="None", type=str,
                    help="path to the spatialhub TSV table describing reference dataset(s) to use for cell annotation")
parser.add_argument("--featureSet", type=str,
                    help="type of features to subset the reference atlas to before training (one of: 'hvg', 'probes' or 'both')")
parser.add_argument("--probesMapping", type=str,
                    help="path to the mapping file matching probe names to Ensembl IDs")

args = parser.parse_args()




############################## SETUP ##############################

# Read in metadata file
path2meta = f"{args.atlasTSV}"
df = pd.read_csv(path2meta, sep = "\t")

# Subset to atlas of interest
df = df[df["atlas_id"] == args.atlasKey]
assert df.shape[0] == 1, "atlas data frame can only have one single row"
df.index = [0]
atlas_key = df["atlas_id"][0]
print(df)

# Define path to pre-processed reference dataset
# which has already been aggregated at the probe level where several genes match the same probe
atlasPath = os.path.join("annot.dir/atlas.dir/", atlas_key + ".h5ad")

# Read in probes mapping file
path2mapping = f"{args.probesMapping}"
df_map = pd.read_csv(path2mapping)
if 'species' in df_map.columns:
    df_map = df_map[df_map['species'] == df['species'][0]]

# Extract distint gene_ids (i.e. as they were before aggregating to probe names)
ens_version = df['ensembl_version'][0]
gene_key = 'gene_v' + f'{ens_version:.0f}'
panel_genes = set(df_map[gene_key].to_list())

# Define HVG batch key
hvg_batch = df['hvg_batch'][0]




############################## TASKS ##############################


print("Importing atlas " + atlasPath)
adata = sc.read_h5ad(atlasPath)

if 'counts' in adata.layers.keys():
  print("'counts' layer already stored in input reference dataset. Converting to CSR matrix.")
  adata.layers['counts'] = csr_matrix(adata.layers['counts'])
else:  # in this case, we'll assume counts are in X slot (but worth a manual check!)
  print("WARNING: no 'counts' layer stored in input reference dataset. Assuming adata.X slot is set to 'counts'")
  adata.layers['counts'] = csr_matrix(adata.X.copy())


# Perform feature selection

print("Defining feature subset to train model on.")

if args.featureSet in ['hvg', 'both']:

    # Retrieve HVGs from adata.var
    if 'highly_variable' in adata.var.columns:
        
        if type(adata.var['highly_variable'].iloc[0]) == str:
            adata.var['highly_variable'] = adata.var['highly_variable'] == 'True'
        
        hvg = set(adata.var.index[adata.var['highly_variable']].to_list())
        print("Number of HVGs in provided AnnData:", len(hvg))

        if len(hvg) > 2 * len(panel_genes):
            
            print("This exceeds 2x number of unique probes in panel (not recommended).",
                  "Re-calculating HVGs with " + hvg_batch + " as batch_key.")
            
            adata.var['original_hvg'] = adata.var['highly_variable']
            
            sc.pp.highly_variable_genes(adata,
                                        n_top_genes=min([2000, 2 * len(panel_genes)]),
                                        layer='counts',
                                        batch_key=hvg_batch,
                                        flavor="seurat_v3_paper",
                                        span=1,
                                        inplace=True)
            
            hvg = adata.var.index[adata.var['highly_variable']].to_list()
            print(len(hvg))

    else:
        print("No HVG found in provided AnnData. Re-calculating with " + hvg_batch + " as batch_key.")
        sc.pp.highly_variable_genes(adata,
                                        n_top_genes=min([2000, 2 * len(panel_genes)]),
                                        layer='counts',
                                        batch_key=hvg_batch,
                                        flavor="seurat_v3_paper",
                                        span=1,
                                        inplace=True)
        hvg = adata.var.index[adata.var['highly_variable']].to_list()
        print(len(hvg))
        
else:
    hvg = []
    
if args.featureSet in ['probes', 'both']:
    probes = df_map['probe_name'].to_list()
else:
    probes = []


features = set(hvg + probes)    
adata = adata[:, adata.var.index.isin(features)]

# Save resulting subset atlas
print("Saving the following subset", adata)
adata.write_h5ad(os.path.join("annot.dir/atlas.dir/", atlas_key + "_feature-subset.h5ad"))

