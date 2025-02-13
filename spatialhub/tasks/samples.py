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

def check_cols(pd_frame, req_columns_list,
               table_name="default"):

    for col_name in req_columns_list:
        if col_name not in pd_frame.columns:
            raise ValueError("required column: '" + col_name + "' missing "
                             "in " + table_name + " table")


def check_values(pd_frame, col, allowed):
    '''utility function for sanity checking columns'''
    if not all([x in allowed for x in pd_frame[col].values]):
        raise ValueError("Only the following values are allowed in column '"
                         + col + "': " + ",".join(allowed))


# ------------------------------ classes -------------------------------- #

class samples():
    '''
    A class for defining the samples and libraries present in
    a spatial experiment.
    '''

    def __init__(self,
                 sample_tsv = None, 
                 library_tsv = None):

        samples = pd.read_csv(sample_tsv, sep="\t")
        required_sample_cols = ["sample_id", "slide_id", "batch",
                                "fov_sequence", 
                                "fov_height", "fov_width"]
        
        check_cols(samples, required_sample_cols,
                   "fov2samples.tsv")
        # < Add check for fov_sequence formatting (; and blank) >
                
        # < Use this when expanding pipeline to different spatial technology types >
        
        #libs = pd.read_csv(library_tsv, sep="\t")
        #required_library_cols = ["technology"]
        
        #check_cols(libs, required_library_cols,
        #           "libraries.tsv")
        
        #self.known_technology = ["cosmx",
        #                         "xenium"]
        
        #check_values(libs, "technology", 
        #             self.known_technology)
        
        # Check for uniqueness of sample names per slide
        x = samples["sample_id"].astype(str) + "_" + samples["slide_id"].astype(str)
        if not x.is_unique:
            raise ValueError("Repeated sample_ids within the same slide for at least one slide_id.")
        
        samples.index = x
        self.samples = samples.to_dict(orient="index")
        
        #libs.sort_values(["library_id",
        #            "technology",
        #            "sample"], inplace=True)
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
    
    def write_tsv(self, group, group_key, outfile_path):
        '''
        Returns a mock input TSV file per group (slide, sample, batch...) to iterate a task of the pipeline on
        '''
        group_df = pd.DataFrame({group: group_key}, index = [0])
        group_df.to_csv(outfile_path, index = False, sep = "\t")        
