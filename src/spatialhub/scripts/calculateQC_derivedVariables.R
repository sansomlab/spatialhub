#!/usr/bin/Rscript --vanilla

# This script calculates derived variables from the output of calculateQC_fromCountsMatrix.R.

message("Computing derived probe detection QC metrics")
timestamp()


# ---------- Libraries ----------

stopifnot(
  require(optparse),
  require(tidyverse),
  require(data.table),
  require(plyr)
)


# ---------- Options ----------

option_list <- list(
    make_option(
        c("--run_dir"),
        required = TRUE,
        help="Directory to look up for per-sample QC metrics calculated from the counts matrix (where new outputs will also be written)"
    ),
    make_option(
        c("--sample_id"),
        required = TRUE,
        help="Sample ID to use to fetch QC metrics data tables"
    ),
    make_option(
        c("--runFOVqc"),
        default = FALSE,
        help="Whether to calculate summary metrics at the FOV-level"
    )
)

opt <- parse_args(OptionParser(option_list=option_list))
set.seed(123456789)

cat("Running with options:\n")
print(opt)




# --------- Setup ----------

cellQC_file <- paste0(opt$run_dir, "/", opt$sample_id, "_cellQCmetrics.csv")
stopifnot("file not found" = file.exists(cellQC_file))
df0 <- read.csv(cellQC_file, header = TRUE)

if (opt$runFOVqc) {

  # import probe classifier (to retrieve number of negative and study probes)
  probes_class_file <- paste0(opt$run_dir, "/", opt$sample_id, "_probeClassifier.csv")
  stopifnot("file not found" = file.exists(probes_class_file))
  probe_class <- read.csv(probes_class_file, header = TRUE)
  head(probe_class)

  # import metadata (to retrieve FOV)
  metadata_file <- paste0(opt$run_dir, "/", opt$sample_id, "_cellMetadata.csv")
  stopifnot("file not found" = file.exists(metadata_file))
  df_meta <- read.csv(metadata_file, header = TRUE)
  head(df_meta)
  
  if (!("FOV" %in% names(df_meta))) {
    
    # If no 'FOV' variable is present in the metadata
    if ("fov" %in% names(df_meta)) {
      # => first, check if it is a case-sensitivity issue
      col_idx <- which(names(df_meta) == "fov")
      names(df_meta)[col_idx] <- "FOV"
    } else if ("instance_id" %in% names(df_meta)) {
      # => next, check if 'instance_id' exists and follows the AtoMx pattern
      #    and extract FOV from it
      x <- df_meta$instance_id[1]
      if( grepl("^c_[0-9]+_[0-9]+_[0-9]+$", x) ) {
        df_meta$FOV <- tstrsplit(df_meta$instance_id, "_")[[3]]
      } else { 
        stop("cannot derive FOV-level QC metrics: 'FOV' variable not found in metadata")
      }
    } else {
      stop("cannot derive FOV-level QC metrics: 'FOV' variable not found in metadata")
    }
    
  }
  
  df_meta <- df_meta |> dplyr::select(instance_id, FOV)
  df0 <- plyr::join(df_meta, df0, by = "instance_id")
  
} 




# --------- FOV-level derived variables ----------

if (opt$runFOVqc) {
  
  print("Calculating FOV-level QC metrics.")
  # Aggregate cell-level metrics at the FOV level
  df <- df0 |> dplyr::select(-instance_id, -sample_id)
  
  # focusing on a subset of quantitative variables that make sense to aggregate at the FOV level
  # NOTE: could be interesting to have these FOV-level metrics for fluorescence data too...
  df <- df |> dplyr::select(FOV, nCount_neg, nCount_sys,
                            nCount_RNA, nFeature_RNA, complexity)
  df$nCell <- 1
  
  df_sum <- aggregate(. ~ FOV, df, sum)
  df_sum <- df_sum[, -which(names(df_sum) %in% c("nFeature_RNA", "complexity"))]
  names(df_sum) <- paste0(names(df_sum), "_FOVsum"); names(df_sum)[1] <- "FOV"
  
  df_avg <- aggregate(. ~ FOV, df, mean)
  df_avg <- df_avg[, -which(names(df_avg) == "nCell")]
  names(df_avg) <- paste0(names(df_avg), "_FOVavg"); names(df_avg)[1] <- "FOV"
  
  dfov <- plyr::join(df_sum, df_avg, by = "FOV")
  head(dfov)
  
  
  # Calculate derived variables as described in AtoMx user manual
  n_neg <- probe_class |> dplyr::filter(negative == TRUE) |> nrow()
  n_probes <- probe_class |> dplyr::filter(!negative & !positive) |> nrow()
  
  # Signal to noise ratio
  # defined as the FOV total number of detected probes normalized to the number of distinct probes in the panel, 
  # divided by the FOV total number negative probes normalized to the number of distinc negative probes.
  dfov$FOV_SNR <- (dfov$nCount_RNA_FOVsum / n_probes) / (dfov$nCount_neg_FOVsum / n_neg)
  
  # If it exists, retrieve existing FOV file and append FOV QC metrics to it
  fovQC_file <- paste0(opt$run_dir, "/", opt$sample_id, "_brukerQCresults.csv")
  if(file.exists(fovQC_file)) {
    dfov0 <- read.csv(fovQC_file, header = TRUE)
    dfov0 <- dfov0 |> dplyr::select(-gene, -failed_cycles)
    dfov0 <- dfov0[!duplicated(dfov0), ] |> dplyr::arrange(FOV)
    dfov <- plyr::join(dfov, dfov0, by = "FOV") |>
      dplyr::relocate(sample_id, .before = FOV)
  }

  # Save FOV-level QC metrics
  print("Saving FOV-level QC metrics:")
  print(head(dfov))
  write.csv(dfov, row.names = FALSE, quote = FALSE, 
            file = paste0(opt$run_dir, "/", opt$sample_id, "_fovQCmetrics.csv"))
 
}




# ---------- Sample-level derived variables ----------

print("Calculating sample-level QC metrics.")
# Aggregate cell-level metrics at the FOV level
df <- df0 |> dplyr::select(-instance_id, -FOV)

# focusing on quantitative variables that make sense to aggregate at the sample level
df <- df |> dplyr::select(sample_id, nCount_neg, nCount_sys,
                          nCount_RNA, nFeature_RNA, complexity)
df$nCell <- 1

df_sum <- aggregate(. ~ sample_id, df, sum)
df_sum <- df_sum[, -which(names(df_sum) %in% c("complexity", "nFeature_RNA"))]
names(df_sum) <- paste0(names(df_sum), "_sampleSum"); names(df_sum)[1] <- "sample_id"

df_avg <- aggregate(. ~ sample_id, df, mean)
df_avg <- df_avg[, -which(names(df_avg) == "nCell")]
names(df_avg) <- paste0(names(df_avg), "_sampleAvg"); names(df_avg)[1] <- "sample_id"

dfs <- plyr::join(df_sum, df_avg, by = "sample_id")


print("Saving sample-level QC metrics:")
head(dfs)
write.csv(dfs, row.names = FALSE, quote = FALSE, 
          file = paste0(opt$run_dir, "/", opt$sample_id, "_sampleQCmetrics.csv"))




# ----------

print("Done calculating QC metrics within the following R environment:")
sessionInfo()