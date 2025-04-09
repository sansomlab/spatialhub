#!/usr/bin/env python3

#import sopa  # Note: spatialdata incompatible with latest dask, hence 'Future Warning'

import os
import sys
import argparse
from spatialdata import SpatialData
from spatialdata import models
import skimage as ski
import pandas as pd
import shapely
import geopandas as gpd


### PARAMETERS ###

print("Parsing arguments")
parser = argparse.ArgumentParser()

parser.add_argument("--projDir", default="None", type=str,
                    help="path to the project directory containing raw data")
parser.add_argument("--stitchedDir", default="./ashlar.dir", type=str,
                    help="path to the project directory containing stitched data")
parser.add_argument("--sampleKey", default="None", type=str,
                    help="unique sample_id as in samples.tsv file, referring to the sample to execute the script on")
parser.add_argument("--fov2sample", default="samples.tsv", type=str,
                    help="path to the spatialhub samples TSV table")
parser.add_argument("--channels", default="PanCK,CD68,Membrane,CD45,DNA", type=str,
                    help="channel dictionary")
parser.add_argument("--rmProbes", default="Neg,Sys,Bac", type=str,
                    help="patterns of probes to exclude from main transcripts table")
parser.add_argument("--probeKey", default="target", type=str,
                    help="name of variable recording probes in transcripts file (e.g. target or gene)")
parser.add_argument("--x", default="x_global_px", type=str,
                    help="names of X spatial coordinate in transcript file")
parser.add_argument("--y", default="y_global_px", type=str,
                    help="names of Y spatial coordinate in transcript file")
parser.add_argument("--importAtomxMask", default=True, type=bool,
                    help="whether to add original AtoMx segmentation mask to sdata object")

args = parser.parse_args()


### SETUP ###

# Read in metadata file
path2meta = f"{args.fov2sample}"
df = pd.read_csv(path2meta, sep = "\t")

if "sample_name" not in df.columns:
    print("'sample_name' not found: setting it to 'sample_id'")
    df["sample_name"] = df["sample_id"]

# Subset to sample of interest
df = df[df["sample_id"] == args.sampleKey]
assert df.shape[0] == 1, "sample data frame can only have one single row"
df.index = [0]
print(df)
slideName = list(set(df["slide_id"]))[0]
sampleName = list(set(df["sample_name"]))[0]
print("Creating spatial data object for sample:", sampleName, slideName)


# Define paths to source and destination directories for corresponding sample
path2files = f"ashlar.dir/{slideName}/Stitched2D/"
path2seg = f"{args.projDir}/{slideName}/atomx_segmentation/"
if args.importAtomxMask:
    assert os.path.exists(path2seg), "Please store AtoMx segmentation mask in a 'atomx_segmentation' sub-directory for the corresponding slide"

zarrDir = f"zarr.dir/{slideName}/"
os.makedirs(zarrDir, exist_ok=True)

outFile = os.path.join(zarrDir, sampleName + ".zarr")
if os.path.exists(outFile):
    print(outFile + " already exists - exiting .zarr creation script.")
    sys.exit()  # Terminate the script with status code 0


# Parse channel dictionary
x = str(args.channels)
x = x.split(",")
x = list(map(str, x))
channels_dict = x
print("Using the following channel dictionary:", channels_dict)

# Parse probes to exclude
x = str(args.rmProbes)
rm_probes_regex = x.replace(",", "|")
print("Using the following regex pattern for probes to exclude:", rm_probes_regex)


### TASKS ###

# Import using SOPA
#sdata = sopa.io.ome_tif(path2img, as_image=False)
# => here, we are losing the channel dictionary...

# Perform tissue-level segmentation, to remove outlier cells
#sopa.segmentation.tissue(sdata, mode = "staining", 
#                         level = 1,  # tissue-level segmentation can be performed on lower-quality image
#                         expand_radius_ratio = 0.05,  # expanding by 5% around defined shape
#                         drop_threshold = 0.01  # shapes covering >1% of the total area will be removed 
#                         )


# ---------- Task 1: Import *stitched* image as array ---------- #

print("Importing source image")
path2img = os.path.join(path2files, sampleName + "_stitched.ome.tiff")
img = ski.io.imread(path2img)
imgs = models.Image2DModel.parse(img, rgb = None, 
                                 c_coords = channels_dict)
    # this step converts the 'raw' image for compatibility with sdata format


# ---------- Task 2: Import (and split) *stitched* transcripts file  ---------- #

print("Importing source transcripts file")
path2tx = os.path.join(path2files, sampleName + "_stitched_tx.csv")
tx = pd.read_csv(path2tx)

# split matrix between main panel, negative/control probes and others (not to be considered in cell segmentation)
# => for example, in case of bacterial probes on mouse panel, exclude these
idx0 = tx[args.probeKey].str.contains('Negative|System|Blank')
tx0 = tx[idx0].copy()
tx0 = models.PointsModel.parse(tx0, feature_key = args.probeKey,
                               coordinates = {'x': args.x, 'y': args.y})
tx = tx[~idx0].copy()  # updating tx to exclude negative/control probes
    

if rm_probes_regex != 'none':
    
    idx_other = tx[args.probeKey].str.contains(rm_probes_regex)
    tx_other = tx[idx_other].copy()
    tx_other = models.PointsModel.parse(tx_other, feature_key = args.probeKey,
                                        coordinates = {'x': args.x, 'y': args.y})
    
    tx = tx[~idx_other].copy()
    txs = models.PointsModel.parse(tx, feature_key = args.probeKey,
                                   coordinates = {'x': args.x, 'y': args.y})
    tx_dict = {'tx_main': txs, 'tx_other': tx_other, 'tx_ctrl': tx0}

else:
    txs = models.PointsModel.parse(tx, feature_key = args.probeKey,
                                   coordinates = {'x': args.x, 'y': args.y})
    tx_dict = {'tx_main': txs, 'tx_ctrl': tx0}


print("Creating SpatialData object")
sdata = SpatialData(images = {'image': imgs}, points = tx_dict)




# ---------- [OPTIONAL] Task 3: Import original segmentation mask  ---------- #

if args.importAtomxMask:

    print("Importing AtoMx default segmentation mask")

    # importing ashlar coordinates file for conversion of segmentation mask coordinates
    path2fov = os.path.join(path2files, sampleName + "_ashlar_fov_positions.csv")
    ashlar_df = pd.read_csv(path2fov)
    ashlar_df = ashlar_df.query("FOV != 'blank'")
    ashlar_df['FOV'] = ashlar_df['FOV'].astype('int')

    fov_ls = ashlar_df['FOV'].to_list()
    fov_ls = [int(x) for x in fov_ls if x != 'blank']
    print(sampleName, 'covers the following FOVs:', fov_ls)

    # import original segmentation mask for corresponding FOVs
    print("Loading segmentation mask for relevant FOVs")

    atomx_ls = []
    for fov in fov_ls:
        # open file
        k = 5 - len(str(fov))
        file = os.path.join(path2seg, 'CellBoundaries_F' + '0'*k + str(fov) + '.csv')
        df_atomx = pd.read_csv(file)
        atomx_ls = atomx_ls + [df_atomx]

    df_atomx = pd.concat(atomx_ls)
    print(df_atomx)


    # Convert coordinates from local to global in newly stitched image: 
    print("Converting x,y shape points to new coordinates system")

    # reformat as one list of paired coordinates per cell, as will be needed for later conversion to geopandas
    df_atomx["atomx_index"] = "FOV" + df_atomx["fov"].astype(str) + "_C" + df_atomx["cellID"].astype(str)
    cell_ids = list(set(df_atomx['atomx_index'].to_list()))

    df_ls = [y for x, y in df_atomx.groupby("atomx_index")]
    #df_ls[0]

    coords_ls = []
    for dfk in df_ls:
        
        # dfk = set of polygon coordinates for cell k
        dfk.index = range(0, len(dfk))  # re-index from zero to n-1

        # update coordinates from local to global, using Ashlar's FOV position file
        fov = dfk['fov'][0]
        f = (ashlar_df['FOV'] == fov)
        x = ashlar_df[f].iloc[0]['Position_X']; y = ashlar_df[f].iloc[0]['Position_Y']
        dfk['x_global'] = dfk['x_local'] + x; dfk['y_global'] = dfk['y_local'] + y

        # store as list for conversion to Shapely polygon
        coords_k = []
        for i in range(0, len(dfk)):
            coords_k = coords_k + [[dfk['x_global'][i], dfk['y_global'][i]]]
        coords_ls = coords_ls + [shapely.Polygon(coords_k)]

    len(coords_ls)
    #coords_ls

    # convert to geopandas format
    print("Transforming to GeoPandas object for compatibility with SpatialData")
    polygon_gdf = gpd.GeoDataFrame(geometry=coords_ls)
    polygon_gdf["atomx_index"] = cell_ids
    polygon_gdf

    # convert to zarr compatible format
    print("Adding shapes element to SpatialData")
    dfs_atomx = models.ShapesModel.parse(polygon_gdf)
    sdata['atomx'] = dfs_atomx




# ---------- Task 4: Write final SpatialData object to disk  ---------- #

print("Writing the following SpatialData object to disk")
print(sdata)
sdata.write(outFile, overwrite=True)
