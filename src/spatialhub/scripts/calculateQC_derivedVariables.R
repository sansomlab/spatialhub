#!/usr/bin/Rscript --vanilla

# This script calculates derived variables from the output of calculateQC_fromCountsMatrix.R.

message("Computing derived probe detection QC metrics")
timestamp()


# ---------- Libraries ----------

stopifnot(
  require(optparse),
  require(dplyr),
  require(data.table)
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
cat("Running with options:\n")
print(opt)




# --------- Setup ----------

instanceQC_file <- paste0(opt$run_dir, "/", opt$sample_id, "_instanceQCmetrics.csv")
stopifnot("Instance QC metrics file not found" = file.exists(instanceQC_file))
df0 <- data.table::fread(instanceQC_file)

if (anyDuplicated(df0$instance_id)) {
  stop("Duplicate instance_id values in instance QC metrics")
}

# Import optional FOV-level metadata if applicable
if (opt$runFOVqc) {

  probes_file <- file.path(
    opt$run_dir, paste0(opt$sample_id, "_probeQCmetrics.csv")
  )
  metadata_file <- file.path(
    opt$run_dir, paste0(opt$sample_id, "_instanceMetadata.csv")
  )
  stopifnot(
    "Probe QC metrics file not found" = file.exists(probes_file),
    "Instance metadata file not found" = file.exists(metadata_file)
  )
  probe_class <- data.table::fread(probes_file)
  df_meta <- data.table::fread(metadata_file)

  # Retrieve FOV identifiers
  if ("FOV" %in% names(df_meta)) {
    fov <- df_meta$FOV
  } else if ("fov" %in% names(df_meta)) {
    fov <- df_meta$fov
  } else {
    ids <- as.character(df_meta$instance_id)
    if (!all(grepl("^c_[0-9]+_[0-9]+_[0-9]+$", ids))) {
      stop("Cannot derive FOV from instance identifiers")
    }
    fov <- data.table::tstrsplit(ids, "_")[[3]]
  }
  df_meta$FOV <- fov
  if (anyDuplicated(df_meta$instance_id) || anyNA(df_meta$FOV)) {
    stop("Duplicate instance IDs or missing FOV assignments")
  }

  df0 <- dplyr::left_join(
    df0,
    df_meta |> dplyr::select(instance_id, FOV),
    by = "instance_id"
  )

  if (anyNA(df0$FOV)) {
    stop("Some instances could not be assigned to an FOV")
  }
}




# --------- FOV-level derived variables ----------

if (opt$runFOVqc) {
  print("Calculating FOV-level QC metrics.")
  
  dfov <- df0 |>
    dplyr::group_by(sample_id, FOV) |>
    dplyr::summarise(
      nInstance_FOVsum = dplyr::n(),
      nCount_neg_FOVsum = sum(nCount_neg),
      nCount_sys_FOVsum = sum(nCount_sys),
      nCount_RNA_FOVsum = sum(nCount_RNA),
      nCount_neg_FOVavg = mean(nCount_neg),
      nCount_sys_FOVavg = mean(nCount_sys),
      nCount_RNA_FOVavg = mean(nCount_RNA),
      nFeature_RNA_FOVavg = mean(nFeature_RNA),
      complexity_FOVavg = mean(complexity, na.rm = TRUE),
      .groups = "drop"
    )
  dfov$FOV <- as.character(dfov$FOV)
  
  # Calculate derived variables as described in AtoMx user manual
  n_neg <- probe_class |> dplyr::filter(negative == TRUE) |> nrow()
  n_probes <- probe_class |> dplyr::filter(!negative & !positive) |> nrow()
  
  if (n_neg == 0 || n_probes == 0) {
    stop("FOV SNR requires both negative and main probes to be present in the probe panel")
  }

  # Signal to noise ratio
  # defined as the FOV total number of detected probes normalized to the number of distinct probes in the panel, 
  # divided by the FOV total number negative probes normalized to the number of distinc negative probes.
  dfov$FOV_SNR <- ifelse(
    dfov$nCount_neg_FOVsum > 0,
    (dfov$nCount_RNA_FOVsum / n_probes) /
      (dfov$nCount_neg_FOVsum / n_neg),
    NA_real_
  )
  
# Append optional Bruker FOV QC results
  bruker_file <- file.path(opt$run_dir, paste0(opt$sample_id, "_brukerQCresults.csv"))
  if (file.exists(bruker_file)) {
    bruker <- data.table::fread(bruker_file)
    if (anyDuplicated(bruker$FOV)) {
      stop("Duplicate FOV identifiers in Bruker QC results")
    }
    bruker$FOV <- as.character(bruker$FOV)
    dfov <- dplyr::left_join(
      dfov,
      bruker |> dplyr::select(-sample_id),
      by = "FOV"
    )
  }

  # Save FOV-level QC metrics
  print("Saving FOV-level QC metrics:")
  print(head(dfov))
  data.table::fwrite(
    dfov,
    file.path(opt$run_dir, paste0(opt$sample_id, "_fovQCmetrics.csv"))
  )
 
}




# ---------- Sample-level derived variables ----------

print("Calculating sample-level QC metrics.")
dfs <- df0 |>
  dplyr::group_by(sample_id) |>
  dplyr::summarise(
    nInstance_sampleSum = dplyr::n(),
    nCount_neg_sampleSum = sum(nCount_neg),
    nCount_sys_sampleSum = sum(nCount_sys),
    nCount_RNA_sampleSum = sum(nCount_RNA),
    nCount_neg_sampleAvg = mean(nCount_neg),
    nCount_sys_sampleAvg = mean(nCount_sys),
    nCount_RNA_sampleAvg = mean(nCount_RNA),
    nFeature_RNA_sampleAvg = mean(nFeature_RNA),
    complexity_sampleAvg = mean(complexity, na.rm = TRUE),
    .groups = "drop"
  )


print("Saving sample-level QC metrics:")
head(dfs)
data.table::fwrite(
  dfs,
  file.path(opt$run_dir, paste0(opt$sample_id, "_sampleQCmetrics.csv"))
)




# ----------

print("Done calculating QC metrics within the following R environment:")
sessionInfo()