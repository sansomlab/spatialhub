#!/usr/bin/Rscript --vanilla

## < Add Description and usage info >

message("Converting original transcripts coordinates to stitched, per-sample coordinate system")
timestamp()


# ---------- Libraries ----------

stopifnot(
  require(optparse),
  require(tidyverse),
  require(data.table)
)


# ---------- Options ----------

option_list <- list(
    make_option(
        c("--projDir"),
        help="path to the project directory containing raw data"
    ),
    make_option(
        c("--sampleKey"),
        help="unique sample_id as in samples.tsv file, referring to the sample to execute the script on"
    ),
    make_option(
        c("--fov2sample"),
        help="path to the spatialhub samples TSV table"
    ),
    make_option(
        c("--pxSize"),
        default=0.120280945,
        help="pixel to micrometer conversion factor"
    )
)

opt <- parse_args(OptionParser(option_list=option_list))
set.seed(123456789)

cat("Running with options:\n")
print(opt)


# ---------- Useful functions ----------

### Read list of transcripts file

readSampleTx <- function(splitTxDir, slide, fovs,
                         simplify = FALSE) {  
    
    # if simplify = TRUE (default), removes cell_ids from AtoMx segmentation mask

    tx_ls <- list()
    for (f in fovs) {
    
    # import FOV raw file
    tx_file <- paste0(splitTxDir, "/", 
                      slide, "_tx_FOV", f, ".csv")
    
    if(file.exists(tx_file)) { 
      
      tx <- read.csv(tx_file) 
      
      # keep atomx cell index, to enable matching dataset if desired
      tx$atomx_index <- paste0("FOV", tx$fov, "_C", tx$cell_ID)
      #head(tx)
      
      # OPTIONAL: exclude negative/control probes
      tx <- tx |> #[, grep("Negative|SystemControl", colnames(tx), invert = TRUE)] |> 
        dplyr::select(-fov, -cell_ID, -cell) |>
        dplyr::rename(atomx_CellComp = CellComp)
      
      if (simplify == TRUE) {
        tx <- tx |> 
            dplyr::select(-atomx_index, -atomx_CellComp)
      }
      
      tx_ls[[paste0("FOV", f)]] <- tx
      
    } else {
      
      print(paste0("WARNING: No data for FOV", f))
      
    }
    
  }

  tx_ls

}

### Convert X/Y transcripts coordinates following image stitching with Ashlar

convertCoordsAshlar <- function(tx_ls, ashlar_df, 
                                px2um = 0.120280945)  # fixed value for CosMx
{

    # takes a list of transcripts files corresponding to one sample
    # returns a per-sample dataframe with corrected coordinates
  
    ashlar_df$FOV <- paste0("FOV", as.character(ashlar_df$FOV))
  
    for (f in names(tx_ls)) {
        
        print(paste0("Converting coords for FOV ", f))
        tx <- tx_ls[[f]]
        print(head(tx))
    
        # Update *global* coordinates (local coordinates within FOV should remain unchanged!)
        # Note: Here, we want to switch to a per-sample coordinate system, so we overwrite the original global coordinates
        tx$x_global_px <- ashlar_df$Position_X[ashlar_df$FOV == f] + tx$x_local_px
        tx$y_global_px <- ashlar_df$Position_Y[ashlar_df$FOV == f] + tx$y_local_px
        print(head(tx))
    
        # < OPTIONAL: Add px to um conversion >
    
        tx_ls[[f]] <- tx
    
    }
  
    print("Binding into one data frame")
    tx_df <- do.call(rbind.data.frame, tx_ls) 
    rownames(tx_df) <- 1:nrow(tx_df)

    tx_df
  
}


# ---------- Setup ----------

# Read in metadata file
path2meta <- opt$fov2sample
df <- read.delim(path2meta, sep = "\t")

if (!("sample_name" %in% names(df))) {
    print("'sample_name' not found: setting it to 'sample_id'")
    df$sample_name <- df$sample_id
}

# Subet to sample of interest
df <- df |> dplyr::filter(sample_id == opt$sampleKey)
print(df)
slideName <- as.character(unique(df$slide_id))
sampleName <- as.character(unique(df$sample_name))

# Extract list of fovs corresponding to sample
fovs <- unlist(strsplit(df$fov_sequence, ","))
fovs <- as.integer(sort(fovs[fovs != "blank"]))

# f"{args.projDir}/data/grouped/{slideName}/ashlar.dir"
ashlarDir <- paste("ashlar.dir", 
                   slideName, "Stitched2D", 
                   sep = "/")
stopifnot(dir.exists(ashlarDir))

# Define path to split AtoMx transcripts file
path2tx <- paste(opt$projDir,
                 slideName, "flatFiles", "split_txFile",
                 sep = "/")
#fileList <- list.files(path2tx, full.names = TRUE) 
#idx <- grep("_tx_file.csv.gz$", fileList)
#tx_file <- fileList[idx]


# ---------- Tasks ----------

# Import list of transcripts files for the necessary FOVs
tx_ls <- readSampleTx(path2tx, slide = slideName, fovs = fovs, 
                      simplify = FALSE)

# Import original Ashlar coordinates file
ashlarPath <- paste0(ashlarDir, "/", sampleName, "_ashlar_fov_positions.csv")
ashlar_df <- read.csv(ashlarPath)
#head(ashlar_df)
#ashlar_df |> dplyr::filter(Position_X == 0)

# Convert AtoMx coordinates to new stitched coordinates system
tx_df <- convertCoordsAshlar(tx_ls, ashlar_df, px2um = opt$pxSize)
head(tx_df)

# Save output
write.csv(tx_df, row.names = FALSE, quote = FALSE,
          file = paste0(ashlarDir, "/", sampleName, "_stitched_tx.csv"))


timestamp()
message("tx file stitching completed.")
