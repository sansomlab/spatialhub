#!/usr/bin/env python3

import os
import re
import numpy as np
import pandas as pd
import anndata as ad
from scipy.sparse import csr_matrix


# ------------------------------ Functions => working from flat files ------------------------------ #

def adata_from_cosmx(path2flatFiles, 
                     slide_id,
                     metadata_vars = ['fov', 'CenterX_global_px', 'CenterY_global_px',  # essential columns for QC report
                                      'Area', 'AspectRatio', 'Width', 'Height',
                                      'Mean.PanCK', 'Max.PanCK', 'Mean.CD68', 'Max.CD68',
                                      'Mean.CD298.B2M', 'Max.CD298.B2M', 'Mean.CD45', 'Max.CD45',
                                      'Mean.DAPI', 'Max.DAPI', 'unassignedTranscripts']):
    '''
    Takes raw flat files from a CosMx experiment, as exported from AtoMx
    Returns a (standard) AnnData object matching the original AtoMx segmentation mask
    '''
    
    print("Generating AnnData from AtoMx-exported flat files with the following obs variables:", metadata_vars)
    flatFiles = os.listdir(path2flatFiles)

    ## Read in counts file

    countsFile = list(filter(lambda x: re.search('_exprMat_file.csv.gz', x), flatFiles))
    assert len(countsFile) == 1, "Provided path to flatFiles should contain one single expression matrix file"
    path2counts = os.path.join(path2flatFiles, countsFile[0])
    counts = pd.read_csv(path2counts)

    # Set index and remove non numeric columns
    counts["atomx_index"] = "FOV" + counts["fov"].astype(str) + "_C" + counts["cell_ID"].astype(str)
    counts.index = counts["atomx_index"]
    counts_mat = counts.drop(columns = ["fov", "cell_ID", "atomx_index"])

    # Convert to sparse matrix
    counts_sparse = csr_matrix(counts_mat.to_numpy(), dtype = 'int')

    # Create AnnData object and set indexes
    adata = ad.AnnData(counts_sparse)
    adata.obs_names = counts_mat.index
    adata.var_names = counts_mat.columns


    ## Read in metadata file

    metaFile = list(filter(lambda x: re.search('_metadata_file.csv.gz', x), flatFiles))
    assert len(countsFile) == 1, "Provided path to flatFiles should contain one single expression matrix file"
    path2meta = os.path.join(path2flatFiles, metaFile[0])
    meta = pd.read_csv(path2meta)

    # Set index and subset to columns of interest
    meta.index = "FOV" + meta["fov"].astype(str) + "_C" + meta["cell_ID"].astype(str)
    meta = meta[metadata_vars]
    meta.insert(loc = 0, column = 'slide_id', value = slide_id)
    adata.obs = meta

    return adata


## Example use, e.g. ahead of running probe QC pipeline

#import spatialhub.tasks.adata_utils as utils

#path2flatFiles = "/users/sansom/tme871/work/cosmx_mouse/data/raw/NL4S3a/flatFiles"
#adata = utils.adata_from_cosmx(path2flatFiles, slide_id = 'NL4S3a')
#path2flatFiles = "/users/sansom/tme871/work/cosmx_mouse/data/raw/NL4S3b/flatFiles"
#bdata = utils.adata_from_cosmx(path2flatFiles, slide_id = 'NL4S3b')

#adata_dict = {
#    'NL4S3a': adata, 'NL4S3b': bdata,
#}
#adata_combined = ad.concat(adata_dict, index_unique = "_")
#adata_combined.write('mouse_ibd_hh_atomx.h5ad', compression = 'gzip')




# ------------------------------ Useful code for updating (spatial) AnnData objects ------------------------------ #

# DRAFT: These need to be written up as proper functions


### 1) Adding per *sample* coordinates to an AnnData object collating several samples ###

# REQUIREMENT: First, run the 2 following pipeline tasks:
# spatialhub ashlar make full -v5 -p20
# spatialhub extractCoords make full -v5 -p20

# Import AnnData to which sample coordinates need to be appended
adata = sc.read_h5ad(path2adata)

# Retrieve per sample coordinates, as re-calculated with `ashlar`
# and extracted from each SpatialData object with `extractCoords`
seg_mask = 'atomx'
coords_files = glob.glob(os.path.join("coordinates.dir/", seg_mask, "*", "*_sampleCoords.csv"))
coords_list = []
for file in coords_files:
    df = pd.read_csv(file)
    coords_list = coords_list + [df]

coords_df = pd.concat(coords_list, axis=0, ignore_index=True)

# Current index is unique at the *sample* level => make it unique at the global level
coords_df.index = coords_df['Unnamed: 0'] + '_' + coords_df['slide_id']  # should work for all seg_mask, but need further testing

# Match indexes in coords_df and AnnData
coords_df = coords_df[coords_df.index.isin(adata.obs_names)]
coords_df = coords_df.reindex(adata.obs_names)
if adata.obs.index.equals(coords_df.index):
    adata.obs['CenterX_sample_px'] = coords_df['CenterX_sample_px']
    adata.obs['CenterY_sample_px'] = coords_df['CenterY_sample_px']

# return adata


### 2) Create composite slide coordinates ###

# Function to create a 'composite' spatial coordinates space,
# thus gathering all samples on one single 'slide' and removing empty space.

# REQUIREMENT: First, extract per-sample coordinates (as enabled by the function above) 
# RECOMMENDED: Run QC to remove exclude low-quality samples from the 'composite' slide
#   => this will remove as much blank space as possible from the new 'composite' coordinates system.
# NOTE: Different analyses may highlight different QC issues and require a different sample filtering
#   (e.g. some samples may be dropped after probe QC, after cell type annotation, after SPIN...)
#   => it may be worth re-creating a composite slide accordingly after each analysis.

# See corresponding R/_clean-up_utils.R function, which would be useful to translate as a python function here


### 3) Converting from spatial coordinates in adata.obsm to a 'spatial' obsm object  ###

# RECOMMENDED: Depending on needs for downstream analyses, use sample-level or composite slide coordinates (as may be generated with the functions above)
# rather than original coordinates (which only make sense for one microscopy slide, not for the whole project if it includes multiple slides).

# Import AnnData to which sample coordinates need to be appended
adata = sc.read_h5ad(path2adata)

x_coord_key = 'CenterX_composite_px'; y_coord_key = 'CenterY_composite_px' 
coords = adata.obs[[x_coord_key, y_coord_key]].copy()
coords = coords.to_numpy()

spatial_key = 'spatial'
    # NOTE: by default, anlyses tools relying on this .obsm slot (e.g. Squidpy or SPIN) expect the key to be named 'spatial'
    #   However, if you want to use this slot for visualization in cellxgene, you need to use 'X_spatial' as a key.
adata.obsm[spatial_key] = coords

# return adata


### 4) Filter AnnData based on a region of interest, e.g. manually drawn in Napari ###

# See corresponding R/_clean-up_utils.R function, which would be useful to translate as a python function here