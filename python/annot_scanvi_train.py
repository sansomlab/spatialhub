#!/usr/bin/env python3

import os
import argparse
import scanpy as sc
import anndata as ad
import pandas as pd
import scvi
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
parser.add_argument("--scviNumWorkers", type=int,
                    help="number of workers for scVI/scANVI")
parser.add_argument("--scVI_pretrain", type=bool,
                    help="whether to first train model using scVI")

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

# Extract distint gene_ids (before aggregating to probe names)
ens_version = df['ensembl_version'][0]
gene_key = 'gene_v' + f'{ens_version:.0f}'
panel_genes = set(df_map[gene_key].to_list())


# scVI settings
scvi.settings.dl_num_workers = args.scviNumWorkers
#scvi.settings.num_threads = 7


# scVI/scANVI model covariates
hvg_batch = df['hvg_batch'][0]
scvi_batch = df['scvi_batch'][0]

x = df['categorical_covar'][0]
if x != 'none':
    scvi_categorical = x.split(",")
else:
    scvi_categorical = None

y = df['continuous_covar'][0]
if y != 'none':
    scvi_continuous = y.split(",")
else:
    scvi_continuous = None

# scVI/scANVI label to use for annotation
scvi_label = df['celltype_annot_key'][0]


# Create out directoruis
outDir = os.path.join("annot.dir/scanvi", atlas_key)
if not os.path.exists(outDir):
    os.mkdir(outDir)

# Define model name
x = df['categorical_covar'][0]
x = x.replace(',', '_')
y = df['continuous_covar'][0]
y = y.replace(',', '_')

model_name = 'BATCH_' + scvi_batch 
if x != 'none':
    model_name = model_name + '_CATEGORICAL_' + x
if y != 'none':
    model_name = model_name + '_CONTINUOUS_' + y

model_name = model_name + '_ANNOT_' + scvi_label
modelDir = os.path.join(outDir, model_name)
if not os.path.exists(modelDir):
    os.mkdir(modelDir)


# Set graphical options
plot_annot = [scvi_batch] + [scvi_label] + [df['lineage_key'][0]] + [df['celltype_other_key'][0]]

plot_covar = [scvi_batch] + [df['sample_key'][0]] + [df['donor_key'][0]]
if scvi_categorical is not None:
    plot_covar = plot_covar + scvi_categorical
if scvi_continuous is not None:
    plot_covar = plot_covar + scvi_continuous




############################## TASKS ##############################


# ---------- Task 1: Load and subset reference AnnData to HVGs and/or probes ---------- #

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
        
        if type(adata.var['highly_variable'][0]) == str:
            adata.var['highly_variable'] = adata.var['highly_variable'] == 'True'
        
        hvg = adata.var.index[adata.var['highly_variable']].to_list()

        if len(hvg) > 2 * len(panel_genes):
            
            print("Number of HVGs in provided AnnData exceeds 2x number of unique probes in panel (not recommended).",
                  "Re-calculating with " + hvg_batch + " as batch_key.")
            
            adata.var['original_hvg'] = adata.var['highly_variable']
            
            sc.pp.highly_variable_genes(adata,
                                        n_top_genes=min([2000, 2 * len(panel_genes)]),
                                        layer='counts',
                                        batch_key=hvg_batch,
                                        flavor="seurat_v3_paper",
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
print("Working with the following subset", adata)


# ---------- Task 2: Define scANVI model, and pre-train a scVI model if applicable ---------- #

# Make sure batch and categorical covariates are indeed categorical
catVarList = [scvi_batch] + scvi_categorical
for var in catVarList:
    adata.obs[var] = adata.obs[var].astype('category')
    print(adata.obs[var])

# scANVI supposedly performs better when starting from model pre-trained with scVI
# in addition, when defining scVI/scANVI model for label transfer, note the different parameters from typical SCVI recommendations - optimized for label transfer
# https://docs.scvi-tools.org/en/stable/tutorials/notebooks/scrna/scarches_scvi_tools.html#train-reference

if args.scVI_pretrain:

    print("Training scVI model")

    # setup AnnData for *scVI* modeling
    scvi.model.SCVI.setup_anndata(adata, 
                                  layer="counts", 
                                  batch_key=scvi_batch, 
                                  categorical_covariate_keys=scvi_categorical,
                                  continuous_covariate_keys=scvi_continuous)
    
    # Define *scVI* model
    scvi_ref = scvi.model.SCVI(
        adata,
        use_layer_norm="both",
        use_batch_norm="none",
        encode_covariates=True,
        dropout_rate=0.2,
        n_layers=2,
        n_latent=10,
        gene_likelihood='zinb'
    )
    
    # Train *scVI* reference
    scvi_ref.train(max_epochs = 50)
    scvi_ref.save(dir_path = modelDir, overwrite = True,
                  prefix = 'scVI_')
    
    # Define *scANVI* model from pre-trained *scVI* one
    scanvi_ref = scvi.model.SCANVI.from_scvi_model(scvi_ref,
                                                   labels_key = scvi_label,
                                                   unlabeled_category = 'unknown',
                                                   linear_classifier=True)  # may avoid over-fitting the training dataset
    
else:
    
    # setup AnnData directly for *scANVI* modeling
    scvi.model.SCANVI.setup_anndata(adata,
                                    layer="counts", 
                                    labels_key = scvi_label, 
                                    unlabeled_category = 'unknown',
                                    batch_key=scvi_batch, 
                                    categorical_covariate_keys=scvi_categorical,
                                    continuous_covariate_keys=scvi_continuous)
    
    # Define *scANVI* model from scratch
    scanvi_ref = scvi.model.SCANVI(
        adata,
        use_layer_norm="both",
        use_batch_norm="none",
        encode_covariates=True,
        dropout_rate=0.2,
        n_layers=2,
        n_latent=10,
        gene_likelihood='zinb',
        linear_classifier=True  # may avoid over-fitting the training dataset
    )


# ---------- Task 3: Train scanVI model as defined above (whether from scratch or from scVI) ---------- #

# For more parameter tweaking, see also: https://discourse.scverse.org/t/scvi-tools-label-transfer-accuracy/1503
print("Training scANVI model")
scanvi_ref.train(max_epochs=30, n_samples_per_label=100)

# Save model

#if args.scVI_pretrain:
#    scanvi_ref.save(dir_path = modelDir, overwrite = True,    
#                    prefix = 'scANVI_from_scVI_')
#else:
#    scanvi_ref.save(dir_path = modelDir, overwrite = True,    
#                    prefix = 'scANVI_')
#scanvi_ref = scvi.model.SCANVI.load(dir_path = modelDir, prefix = 'scANVI_', adata = adata)

scanvi_ref.save(dir_path = modelDir, overwrite = True)


# Explore latent space
print("Finding neighbours in latent space")
adata.obsm['X_scANVI'] = scanvi_ref.get_latent_representation()
sc.pp.neighbors(adata, use_rep='X_scANVI')
sc.tl.umap(adata)
adata.write_h5ad(os.path.join(modelDir, atlas_key + '_feature-subset.h5ad'))


# Generate UMAP figures to assess training performance
sc.settings.figdir = modelDir

sc.pl.umap(
    adata,
    color=plot_covar,
    frameon=False,
    ncols=2,
    save="_trained_reference_covar.pdf"
)

sc.pl.umap(
    adata,
    color=plot_annot,
    frameon=False,
    ncols=2,
    save="_trained_reference_annot.pdf"
)
