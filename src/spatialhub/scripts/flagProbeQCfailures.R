#!/usr/bin/Rscript --vanilla

# This script flags cells, samples and FOVs (if applicable) that failed QC.
# Threshold for pass/fail are defined by the user in a YAML file.
# All QC metrics were pre-computed in prior scripts.

message("Flagging QC failures")
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
        help="Directory to look up for per-sample QC metrics calculated"
    ),
    make_option(
        c("--metadata_file"),
        required = TRUE,
        help="Path to CSV file containing sample-level metadata for the study"
    ),
    make_option(
        c("--run_name"),
        required = TRUE,
        help="Name of the QC run (used to name outputs)"
    ),
    make_option(
        c("--outdir"),
        default = "current",
        help="Name of sub-directory where to save outputs from applying the specified QC thresholds"
    ),
    make_option(
        c("--qc_thresholds"),
        required = TRUE,
        help = "YAML file defining probe QC thresholds"
    ),
    make_option(
        c("--runFOVqc"),
        default = FALSE,
        help="Whether FOV QC was enabled for this QC run. If TRUE, FOV QC metrics will be imported and used to flag FOVs that failed QC."
    ),
    make_option(
        c("--runBrukerFOVqc"),
        default = FALSE,
        help="Whether Bruker-specific FOV QC was enabled for this QC run. If TRUE, additional Bruker-specific FOV QC metrics will be imported and used to flag FOVs that failed QC."
    )
)

opt <- parse_args(OptionParser(option_list=option_list))

cat("Running with options:\n")
print(opt)

qc <- yaml::read_yaml(opt$qc_thresholds)




# --------- Setup ----------

# Import QC tables for all samples and append at project level
message("Importing cell QC metrics tables.")
qc_tables <- list.files(path = opt$run_dir, pattern = "_cellQCmetrics\\.csv$",
                        full.names = TRUE)
if (length(qc_tables) == 0) {
  stop("No cell QC metrics tables found in the specified run directory.")
} else {
  qc_cell <- rbindlist(lapply(qc_tables, fread))
}

message("Importing cell metadata tables.")
meta <- list.files(path = opt$run_dir, pattern = "_cellMetadata\\.csv$",
                   full.names = TRUE)
if (length(meta) == 0) {
  stop("No cell QC metadata tables found in the specified run directory.")
} else {
  qc_meta <- rbindlist(lapply(meta, fread))
  #qc_meta <- qc_meta |> dplyr::select(sample_id, instance_id, FOV)
}
stopifnot(
  nrow(qc_meta) == nrow(dplyr::distinct(qc_meta, sample_id, instance_id))
)
qc_cell <- dplyr::left_join(qc_cell, qc_meta, by = c("sample_id", "instance_id"))


if(opt$runFOVqc) {
  message("FOV QC was enabled for this run. Importing FOV QC metrics tables.")
  qc_tables <- list.files(path = opt$run_dir, pattern = "_fovQCmetrics\\.csv$",
                          full.names = TRUE)
  if (length(qc_tables) == 0) {
    stop("No FOV QC metrics tables found in the specified run directory. Please run FOV QC to generate these.")
  } else {
    qc_fov <- rbindlist(lapply(qc_tables, fread))
  }
}

message("Importing sample QC metrics tables.")
qc_tables <- list.files(path = opt$run_dir, pattern = "_sampleQCmetrics.csv",
                        full.names = TRUE)
if (length(qc_tables) == 0) {
  stop("No sample QC metrics tables found in the specified run directory.")
} else {
  qc_sample <- rbindlist(lapply(qc_tables, fread))
  sample_vars <- names(qc_sample)
}

message("Importing user-defined sample-level metadata.")
qc_meta <- read.csv(opt$metadata_file, header = TRUE)
# Subset to included samples
qc_meta <- qc_meta[qc_meta$sample_id %in% qc_sample$sample_id, ]
stopifnot(
  nrow(qc_meta) == nrow(dplyr::distinct(qc_meta, sample_id))
)
qc_sample <- dplyr::left_join(qc_sample, qc_meta, by = c("sample_id"))


### Add cell-level QC flags

message("Adding cell-level QC flags.")
qc_cell$qcFlagCell_nCount <- ifelse(qc_cell$nCount_RNA >= qc$cell_nCount_cutoff, "Pass", "Fail")
qc_cell$qcFlagCell_nFeature <- ifelse(qc_cell$nFeature_RNA >= qc$cell_nFeature_cutoff, "Pass", "Fail")
qc_cell$qcFlagCell_percentNeg <- ifelse(qc_cell$percentNegCounts < qc$cell_percentNeg_cutoff, "Pass", "Fail")
qc_cell$qcFlagCell_complex <- ifelse(qc_cell$complexity > qc$cell_complexity_cutoff, "Pass", "Fail")

cell_pass_qc <- qc_cell$nCount_RNA >= qc$cell_nCount_cutoff &
  qc_cell$nFeature_RNA >= qc$cell_nFeature_cutoff &
  qc_cell$percentNegCounts < qc$cell_percentNeg_cutoff &
  qc_cell$complexity > qc$cell_complexity_cutoff
qc_cell$qcFlagCell_summary <- ifelse(cell_pass_qc, "Pass", "Fail")


### Add FOV-level QC flags

dfc <- qc_cell

if(opt$runFOVqc) {

  message("Adding FOV-level QC flags.")
  qc_fov$qcFlagFOV_nCount <- ifelse(qc_fov$nCount_RNA_FOVavg >= qc$fov_nCount_cutoff, "Pass", "Fail")
  qc_fov$qcFlagFOV_SNR <- ifelse(qc_fov$FOV_SNR > qc$fov_SNR_cutoff, "Pass", "Fail")

  if (opt$runBrukerFOVqc) {
  FOV_pass_qc <- qc_fov$nCount_RNA_FOVavg >= qc$fov_nCount_cutoff &
      qc_fov$FOV_SNR > qc$fov_SNR_cutoff &
      qc_fov$qcFlagGeneBias == "Pass"
  } else {
  FOV_pass_qc <- qc_fov$nCount_RNA_FOVavg >= qc$fov_nCount_cutoff &
      qc_fov$FOV_SNR > qc$fov_SNR_cutoff
  }
  qc_fov$qcFlagFOV_summary <- ifelse(FOV_pass_qc, "Pass", "Fail")

  message("Writing FOV QC metrics table to file.")
  write.csv(qc_fov, row.names = FALSE, quote = FALSE, 
            file = paste0(opt$outdir, "/", opt$run_name, "_fovQCmetrics.csv"))

  # Append FOV QC flags to cell-level QC table (necessary for sample-level QC)
  stopifnot( nrow(qc_fov) == nrow(dplyr::distinct(qc_fov, sample_id, FOV)) )
  dfc <- dplyr::left_join(qc_cell, qc_fov, by = c("sample_id", "FOV"))

  if (any(is.na(dfc$qcFlagFOV_summary) & !is.na(dfc$FOV))) {
    stop("Some cells could not be matched to an FOV QC record.")
  }

}


### Add sample-level QC flags

message("Adding sample-level QC flags.")
# Note: Sample-level QC requires integration of QC metrics across all (available) levels
stopifnot( nrow(qc_sample) == nrow(dplyr::distinct(qc_sample, sample_id)) )
dfa <- dplyr::left_join(dfc, qc_sample, by = c("sample_id"))

# Tally the cells that passed QC for each sample
dfa_ls <- split(dfa, f = dfa$sample_id)

for (sample in names(dfa_ls)) { 
  df <- dfa_ls[[sample]]
    
  # In summary, how many cells for that sample pass QC?
  # i.e. how many cells have enough counts, etc.
  sample_cell_pass_qc <- df$qcFlagCell_summary == "Pass"
  df$nCell_sampleSum_passQC <- sum(sample_cell_pass_qc)

  if (opt$runFOVqc) {

    # Tally the cells that pased *cell* AND *FOV* QC for each sample
    s <- sum(df$qcFlagFOV_summary == "Fail")
    df$samplePercentCellsInFailedFOV <- round(s / nrow(df), 4) * 100
     
    # In summary, how many cells for that sample pass all QC?
    # i.e. how many cells have enough counts, etc. AND fall within a high-quality FOV
    sample_cell_pass_qc <- df$qcFlagCell_summary == "Pass" &
      df$qcFlagFOV_summary == "Pass"
    df$nCellxnFOV_sampleSum_passQC <- sum(sample_cell_pass_qc)
  }

  dfa_ls[[sample]] <- df

}
dfa <- do.call(rbind.data.frame, dfa_ls)

if (opt$runFOVqc && qc$exclCellsInFailedFOVs) {
  # Define summary sample QC metric
  # i.e. sample with a minimum absolute count of cells and percent of cells passing QC
  # based on cell QC and FOV QC both
  dfa$samplePercentCellsPassQC <- round(dfa$nCellxnFOV_sampleSum_passQC/dfa$nCell_sampleSum, 4) * 100
  sample_pass_qc <- dfa$nCellxnFOV_sampleSum_passQC >= qc$sampleMinCells & dfa$samplePercentCellsPassQC > qc$sampleMinPercent
 } else {
  # Define summary sample QC metric based on cell QC only
  dfa$samplePercentCellsPassQC <- round(dfa$nCell_sampleSum_passQC/dfa$nCell_sampleSum, 4) * 100
  sample_pass_qc <- dfa$nCell_sampleSum_passQC >= qc$sampleMinCells & dfa$samplePercentCellsPassQC > qc$sampleMinPercent
}
dfa$qcFlagSample_summary <- ifelse(sample_pass_qc, "Pass", "Fail")

# Reducing data frame back again to one row per sample
keepVar <- c(sample_vars, 
             "nCell_sampleSum_passQC", 
             "nCellxnFOV_sampleSum_passQC", "samplePercentCellsInFailedFOV", 
             "samplePercentCellsPassQC", "qcFlagSample_summary")
keepVar <- names(dfa)[names(dfa) %in% keepVar]
keepVar <- keepVar[keepVar != "FOV"]
dfs <- dfa |> dplyr::select(all_of(keepVar))
head(dfs)
dfs <- dfs[!duplicated(dfs), ]
dfs <- dfs[!is.na(dfs$sample_id), ]
rownames(dfs) <- 1:nrow(dfs)

message("Writing sample-level QC metrics table to file.")
write.csv(dfs, row.names = FALSE, quote = FALSE, 
          file = paste0(opt$outdir, "/", opt$run_name, "_sampleQCmetrics.csv"))


### Writing final QC pass/fail decision rule

# REMINDER: We don't filter FOVs altogether to avoid patchy samples.
# Poor quality FOVs are excluded at the sample level.
dfa$qcFlagOverall_summary <- ifelse(
    dfa$qcFlagCell_summary == "Pass" &
    dfa$qcFlagSample_summary == "Pass",
    "Pass", "Fail"
)

message("Writing *enriched* cell QC metrics table to file.")
write.csv(dfa, row.names = FALSE, quote = FALSE, 
          file = paste0(opt$outdir, "/", opt$run_name, "_cellQCmetrics.csv"))




# ----------

print("Done flagging QC failures within the following R environment:")
sessionInfo()