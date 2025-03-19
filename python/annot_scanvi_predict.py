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

#parser.add_argument("--path2sdata", default="None", type=str,
#                    help="path to spatial AnnData file on which to predict labels")
parser.add_argument("--atlasKey", default="None", type=str,
                    help="key of reference atlas to use for cell typing")
parser.add_argument("--atlasTSV", default="None", type=str,
                    help="path to the spatialhub TSV table describing reference dataset(s) to use for cell annotation")
#parser.add_argument("--featureSet", type=str,
#                    help="type of features to subset the reference atlas to before training (one of: 'hvg', 'probes' or 'both')")
#parser.add_argument("--probesMapping", type=str,
#                    help="path to the mapping file matching probe names to Ensembl IDs")
#parser.add_argument("--batchKey_hvg", type=str,
#                    help="batch key to use when re-calculating HVGs, if applicable")
parser.add_argument("--scviNumWorkers", type=int,
                    help="number of workers for scVI/scANVI")
#parser.add_argument("--scVI_pretrain", type=bool,
#                    help="whether to first train model using scVI")

args = parser.parse_args()




############################## SETUP ##############################

# Read in metsdata file
path2meta = f"{args.atlasTSV}"
df = pd.read_csv(path2meta, sep = "\t")

# Extract query dataset information
df_query = df[df['type'] == "query"]
assert df_query.shape[0] == 1, "atlas data table can only list one single query (spatial) dataset"
df_query.index = [0]
print("Working with the following query (spatial) dataset to annotate: ")
print(df_query)

# Subset to atlas of interest
df = df[df["atlas_id"] == args.atlasKey]
assert df.shape[0] == 1, "atlas data frame can only have one single row"
df.index = [0]
atlas_key = df["atlas_id"][0]
print("Working with the following reference atlas for annotation: ")
print(df)

# Define path to pre-processed query dataset
queryPath = df_query['path'][0]

# Read in probes mapping file
#path2mapping = f"{args.probesMapping}"
#df_map = pd.read_csv(path2mapping)
#panel_probes = set(df_map['gene_v' + df['ensembl_version'][0].astype(str)].to_list())


# scVI settings
scvi.settings.dl_num_workers = args.scviNumWorkers
#scvi.settings.num_threads = 7


# scVI/scANVI model covariates
scvi_batch = df['scvi_batch'][0]

x = df['categorical_covar'][0]
if x != 'none':
    scvi_categorical = x.split(",")
    query_categorical = df_query['categorical_covar'][0].split(",")
    dict_categorical = dict(zip([df_query['scvi_batch'][0]] + query_categorical,
                                [scvi_batch] + scvi_categorical))
else:
    scvi_categorical = None
    dict_categorical = dict(zip([df_query['scvi_batch'][0]], [scvi_batch]))

y = df['continuous_covar'][0]
if y != 'none':
    scvi_continuous = y.split(",")
    query_continuous = df_query['continuous_covar'][0].split(",")
    dict_continuous = dict(zip(query_continuous, scvi_continuous))
else:
    scvi_continuous = None
    dict_continuous = None

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


# ---------- Task 1: Load and check spatial AnnData is ready for annotation ---------- #

print("Importing spatial dataset " + queryPath)
sdata = sc.read_h5ad(queryPath)

if 'counts' in sdata.layers.keys():
  print("'counts' layer already stored in input reference dataset. Converting to CSR array.")
  sdata.layers['counts'] = csr_matrix(sdata.layers['counts'])
else:  # in this case, we'll assume counts are in X slot (but worth a manual check!)
  print("WARNING: no 'counts' layer stored in input reference dataset. Assuming sdata.X slot is set to 'counts'")
  sdata.layers['counts'] = csr_matrix(sdata.X.copy())


# Make sure batch and categorical covariates exists, and are indeed categorical

print("Matching categorical variable names between reference and query")
sdata.obs = sdata.obs.rename(columns = dict_categorical) 
    # dictionary includes main scVI batch_key, and thus cannot be None

for key in dict_categorical.keys():
    # if no matching key in sdata for a covariate in sdata used to train model,
    # then create a mock categorical variable with fixed value (N/A)
    if key == 'none':
        cat_var = dict_categorical[key]
        sdata.obs[cat_var] = 'not applicable'

for value in dict_categorical.values():
    sdata.obs[value] = sdata.obs[value].astype('category')
    print(sdata.obs[value])

# continuous covariates are optional, therefore this dictionary may be empty
if dict_continuous is not None:
    print("Matching continuous variable names between reference and query")
    sdata.obs = sdata.obs.rename(columns = dict_continuous)
    for key in dict_continuous.keys():
        # if no matching key in sdata for a covariate in sdata used to train model,
        # then create a mock continuous variable with fixed value (0)
        if key == 'none':
            cont_var = dict_continuous[key]
            sdata.obs[cont_var] = 0

print("Working with the following spatial dataset", sdata)


# ---------- Task 2: Update model with query dataset and predict cell types ---------- #

# Prepare query using SCANVI model trained on reference dataset
scvi.model.SCANVI.prepare_query_anndata(sdata, modelDir)
scanvi_query = scvi.model.SCANVI.load_query_data(sdata, modelDir)

# Predict cell type in query dataset 
scanvi_query.train(
    max_epochs=100,
    plan_kwargs={"weight_decay": 0.0},
    check_val_every_n_epoch=10
)

sdata.obs['scANVI_predict_' + atlas_key] = scanvi_query.predict(soft=True)


# ---------- Task 3: Get marker genes ---------- #

# Obtain cell type-specific markers
#sc.tl.rank_genes_groups(sdata, groupby=SCANVI_PREDICTIONS_KEY, method="t-test")
#sc.get.rank_genes_groups_df(sdata, group=None).to_csv(SCANVI_MODEL_PATH + "/" + outFileBase + "_top-markers_t-test.csv")
#sc.tl.rank_genes_groups(sdata, groupby=SCANVI_PREDICTIONS_KEY, method="wilcoxon")
#sc.get.rank_genes_groups_df(sdata, group=None).to_csv(SCANVI_MODEL_PATH + "/" + outFileBase + "_top-markers_wilcox.csv")

# Generate dot plot of top markers
#sc.pl.rank_genes_groups_dotplot(
#    sdata, groupby=SCANVI_PREDICTIONS_KEY, standard_scale="var", n_genes=5,
#    save=outFileBase + ".pdf"
#)


# Clean up metsdata and write output
#sdata.obs.rename(columns={SPX_MOD_BATCH: "scanvi_batch"}, inplace = True)

#col_names = list(sdata.obs.columns)
#drop_cols = filter(lambda x: re.search(r'^X_scvi_', x), col_names)
#drop_cols = drop_cols + [SCANVI_LABELS_KEY]
#sdata.obs.drop(columns = drop_cols, inplace = True)

# Save
#sdata.write_h5ad(SCANVI_MODEL_PATH + "/" + outFileBase + ".h5ad")
sdata.obs.to_csv(os.path.join(modelDir, 'scANVI_predictions.csv'))
