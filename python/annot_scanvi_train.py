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
atlasPath = os.path.join("annot.dir/atlas.dir/", atlas_key + "_feature-subset.h5ad")


# scVI settings
scvi.settings.dl_num_workers = args.scviNumWorkers
#scvi.settings.num_threads = 7


# scVI/scANVI model covariates
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


# Create out directories
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

plot_annot = list(filter(lambda x: x != 'none', plot_annot))
plot_covar = list(filter(lambda x: x != 'none', plot_covar))




############################## TASKS ##############################


# Load reference AnnData (subset to HVGs and/or probes)
print("Importing atlas " + atlasPath)
adata = sc.read_h5ad(atlasPath)
print("Training on the following reference dataset", adata)


# ---------- Task 1: Define scANVI model, and pre-train a scVI model if applicable ---------- #

# Make sure batch and categorical covariates are indeed categorical
catVarList = [scvi_batch] + scvi_categorical
for var in catVarList:
    adata.obs[var] = adata.obs[var].astype(str).astype('category')
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
    scvi_ref.train(max_epochs = 30)
    scvi_ref.save(dir_path = modelDir, overwrite = True,
                  prefix = 'scVI_')
    
    # Define *scANVI* model from pre-trained *scVI* one
    scanvi_ref = scvi.model.SCANVI.from_scvi_model(scvi_ref,
                                                   labels_key = scvi_label,
                                                   unlabeled_category = 'unknown',
                                                   linear_classifier = False)  # Setting to True may avoid over-fitting the training dataset
    
else:
    
    # setup AnnData directly for *scANVI* modeling
    scvi.model.SCANVI.setup_anndata(adata,
                                    layer="counts", 
                                    labels_key = scvi_label, 
                                    unlabeled_category = 'unknown',
                                    batch_key = scvi_batch, 
                                    categorical_covariate_keys = scvi_categorical,
                                    continuous_covariate_keys = scvi_continuous)
    
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
        linear_classifier=False  # Setting to True may avoid over-fitting the training dataset
    )


# ---------- Task 2: Train scanVI model as defined above (whether from scratch or from scVI) ---------- #

# For more parameter tweaking, see also: https://discourse.scverse.org/t/scvi-tools-label-transfer-accuracy/1503
print("Training scANVI model")
scanvi_ref.train(max_epochs=50, n_samples_per_label=100)

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
adata.write_h5ad(os.path.join(modelDir, atlas_key + '_ref-adata.h5ad'))


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
