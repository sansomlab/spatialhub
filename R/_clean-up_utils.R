#!/usr/bin/Rscript --vanilla

# Script to generate a composite microscopy slide for a cleaned-up dataset
# i.e. removing blank spaces left on the original microscopy slide after sample QC filtering
# and re-shuffling samples in an order that is more intuitive for visualization and/or MuSpAn analyses
# (e.g. one slide per contrast group, per donor, etc.)

# This script takes as an input a metadata table covering all cells and samples to be considered in the study
# and including `CenterX/Y_sample_px` *stitched* coordinates for each sample

#spx_file <- "./annotation/HNSCC_baysor/scANVI_conflict_resolution/HNSCC_baysor_metadata_scANVI_tumor_cellbender_annotations.csv"
#spx_meta <- read.csv(spx_file)
#spx_meta <- spx_meta |> dplyr::arrange(sample_index)

####### Function to apply a ROI filter to the dataset #######

# Filter dataset for cells within ROI, if applicable (and if not performed before)

filterROI <- function(spx_meta, sample_key = NA, 
                      path2roi, roi_key, cell_key) {

    # spx_meta: spatial metadata object
    # sample_key: name of column in spatial metadata object where the unique sample identified
    #   (taking into account the source slide) is stored
    #   Defaults to 'sample_name' + 'slide_id' if not specified
    # path2roi: path to a directory containing CSV tables of cell 'barcodes' (IDs) within the ROI of interest
    #   e.g. output from the spatialhub roi (`roi.dir`)
    # roi_key: key of the ROI to apply, as specified in the name of the ROI 'barcodes' file
    # cell_key: key used to identify a cell in the corresponding segmentation mask (e.g. 'baysor_index')

    print(paste0("Applying ROI filter."))
    
    roi_files <- list.files(path2roi, recursive = TRUE, full.names = TRUE)
    roi_files <- roi_files[grep(paste0(roi_key, "_barcodes.csv"), roi_files)]

    roi_ls <- list()
    for (file in roi_files) {
        roi_ls[[file]] <- read.csv(file)
    }
    roi_df <- do.call(rbind.data.frame, roi_ls)
    
    if (is.na(sample_key)) {
        roi_df$sample_index <- paste0(roi_df$sample_name, "_", roi_df$slide_id)
    } else {
        roi_df$sample_index <- roi_df[, sample_key]
    }
    roi_df <- roi_df |> dplyr::filter(sample_index %in% spx_meta$sample_index)
    
    # Update cell 'barcode' ID to match unique index appended by anndata when concatenating adatas from different samples
    spx_meta$segment_id <- paste0(spx_meta[, cell_key], "_", spx_meta$sample_index)
    roi_df$segment_id <- paste0(roi_df$id, "_", roi_df$sample_index)
    #roi_df <- plyr::join(roi_df, spx_meta[, c("X", "segment_id")], 
    #                     by = "segment_id", match = "first")
    
    spx_meta <- spx_meta |> 
        dplyr::filter(segment_id %in% roi_df$segment_id)
    
    spx_meta

}


####### Function to generate a composite slide on the full dataset #######

makeComposite <- function(spx_meta, sample_key = NA,
                          spacing = 500, n_col, composite_id = '0') {

    # spx_meta: spatial metadata object for which to create a composite slide
    #   =>  filter for samples to include on this composite slide only if creating several slides
    #       and make sure the metadata is sorted in the correct order
    # sample_key: name of column in spatial metadata object where the unique sample identified
    #   (taking into account the source slide) is stored
    #   Defaults to 'sample_name' + 'slide_id' if not specified
    # spacing: blank space to leave between samples (in same unit as cell coordinates)
    # n_col: number of columns in the composite slide
    # composite_id: identifier for the composite slide, if splitting the dataset into several

    if (is.na(sample_key)) {
        spx_meta$sample_index <- paste0(spx_meta$sample_name, "_", spx_meta$slide_id)
    } else {
        spx_meta$sample_index <- spx_meta[, sample_key]
    }
  
    # Split dataset by sample with its stitched coordinates
    spx_meta_ls <- split.data.frame(spx_meta, f = spx_meta$sample_index)
    #names(spx_meta_ls)

    # Define grid size to place samples on
    u <- c(max(spx_meta$CenterX_sample_px), max(spx_meta$CenterY_sample_px))
    # Let's place samples on a grid of U x U pixels per sample
    grid_size <- plyr::round_any(ceiling(max(u) + spacing), spacing*2, f = ceiling) 
    print("Per sample grid-size (parameter necessary for MuSpAn QCM analysis):")
    print(grid_size)

    for (i in 1:length(spx_meta_ls)) {
    
        df_tp <- spx_meta_ls[[i]]
        
        df_tp$CenterX_composite_px <- df_tp$CenterX_sample_px + grid_size * ((i-1) %% n_col)
        df_tp$CenterY_composite_px <- df_tp$CenterY_sample_px + grid_size * floor((i-1) / n_col)
        
        spx_meta_ls[[i]] <- df_tp
    
    }
    
    spx_meta <- do.call(rbind.data.frame, spx_meta_ls)
    spx_meta$composite_slide <- composite_id

    spx_meta

}

#spx_meta <- makeComposite(spx_meta, sample_key = 'sample_index', 
#                          spacing = 2500, n_col = 8)
write.csv(spx_meta[, c("X", "segment_id", "CenterX_composite_px", "CenterY_composite_px", "composite_slide")],
          file = "./coordinates/HNSCC_baysor_composite_full-dataset.csv",
          row.names = FALSE, quote = FALSE)