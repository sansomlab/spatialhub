'''
samples.py
======

Overview
--------

- read in samples to pandas table
- read in libraries to pandas table

'''

import yaml
import os
import shutil
import re
import copy
import pandas as pd
from pprint import pprint


# ------------------------------ utility functions -------------------------------- #

def check_cols(pd_frame, req_columns_list, table_name="default"):

    for col_name in req_columns_list:
        if col_name not in pd_frame.columns:
            raise ValueError("required column: '" + col_name + "' missing "
                             "in " + table_name + " table")


def check_values(pd_frame, col, allowed):
    '''utility function for sanity checking columns'''

    if not all([x in allowed for x in pd_frame[col].values]):
        raise ValueError("Only the following values are supported in column '"
                         + col + "': " + ", ".join(allowed))


# ------------------------------ classes -------------------------------- #

class samples():
    '''
    A class for defining the samples and libraries present in
    a spatial experiment.
    '''

    def __init__(self, sample_tsv = None):

        samples = pd.read_csv(sample_tsv, sep="\t")
        required_sample_cols = ["sample_id", "slide_id", "batch",
                                "fov_sequence", 
                                "fov_height", "fov_width"]
        
        if "sample_name" not in samples.columns:
            print("'sample_name' not found: setting it to 'sample_id'")
            samples["sample_name"] = samples["sample_id"]
        
        check_cols(samples, required_sample_cols, "samples.tsv")

        # Check for uniqueness of sample names per slide
        x = samples["sample_id"].astype(str) #+ "_" + samples["slide_id"].astype(str)
        if not x.is_unique:
            raise ValueError("Repeated sample_ids, please make sure these are unique (to identify same samples across different slides, use sample_name).")
        
        samples.index = x
        
        # < Add check for fov_sequence formatting (; and blank) >
        # < Use this when expanding pipeline to different spatial technology types >
        
        self.samples = samples.to_dict(orient="index")
        self.libs = samples.to_dict()

        
    def slide_ids(self):
        '''
        Returns a list of slide_ids for the project
        '''
        return(set(self.libs["slide_id"].values()))
        
    def sample_ids(self):
        '''
        Returns a list of sample_ids for the project
        '''
        return(set(self.libs["sample_id"].values()))
    
#    def slidexsample_name(self, sample_id):
#        '''
#        Extracts the slide sample names for a given sample_id
#        '''
#        self_subset = self[self['sample_id'] == 'sample_id']
#        sampleName = self_subset['sample_name']
#        slideName = self_subset['slide_id']
#        return({'sample': sampleName, 'slide': slideName})

#    def write_sample(self, outfile_path):
#        '''
#        Returns an input TSV file per sample_id, with matching variables to be used as input for the pipeline task
#        '''
#        group_df = pd.DataFrame({group: group_key}, index = [0])
#        group_df.to_csv(outfile_path, index = False, sep = "\t")


class atlases():
    '''
    A class for defining the reference atlases to be used for 
    cell type annotation.
    '''

    def __init__(self, atlas_tsv = None):

        refs = pd.read_csv(atlas_tsv, sep="\t")
        refs = refs[refs['type'] == "reference"]

        required_ref_cols = ["atlas_id", "path",
                             "celltype_annot_key",
                             "ensembl_version"]
        
        check_cols(refs, required_ref_cols, "atlas.tsv")

        # Check for uniqueness of sample names per slide
        x = refs["atlas_id"].astype(str)
        if not x.is_unique:
            raise ValueError("Repeated atlas_ids, please make sure these are unique.")
        
        refs.index = x
        
        check_values(refs, "ensembl_version", [87, 93, 98, 110])
        check_values(refs, "species", ["mouse", "human"])
        
        self.refs = refs.to_dict(orient="index")
        self.libs = refs.to_dict()
            
    def atlas_ids(self):
        '''
        Returns a list of atlas_ids for the project
        '''
        return(set(self.libs["atlas_id"].values()))
