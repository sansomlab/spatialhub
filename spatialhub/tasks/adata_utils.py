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




# ------------------------------ Functions => working from spatial data object ------------------------------ #