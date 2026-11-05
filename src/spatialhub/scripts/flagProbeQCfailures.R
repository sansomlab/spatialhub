#!/usr/bin/Rscript --vanilla

# This script flags instances, samples and FOVs (if applicable) that failed QC.
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

if (opt$runBrukerFOVqc && !opt$runFOVqc) {
  stop("Bruker FOV QC requires FOV QC to be enabled")
}


# --------- Setup ----------

read_qc_tables <- function(pattern) {
  paths <- list.files(
    opt$run_dir, pattern = pattern, full.names = TRUE
  )
  if (!length(paths)) {
    stop("No QC files found matching: ", pattern)
  }
  data.table::rbindlist(lapply(paths, data.table::fread))
}

# Import QC tables for all samples and append at project level
message("Importing instance-level QC tables.")
qc_instance <- read_qc_tables("_instanceQCmetrics\\.csv$")
instance_meta <- read_qc_tables("_instanceMetadata\\.csv$")

if (anyDuplicated(qc_instance[, c("sample_id", "instance_id")]) ||
    anyDuplicated(instance_meta[, c("sample_id", "instance_id")])) {
  stop("Duplicate sample_id x instance_id combinations")
}
if (opt$runFOVqc && !"FOV" %in% names(instance_meta)) {
  stop(
    "FOV QC enabled but instance metadata has no FOV assignments"
  )
}

qc_instance <- dplyr::left_join(qc_instance, instance_meta, by = c("sample_id", "instance_id"))

# Import FOV QC metrics if applicable
if(opt$runFOVqc) {
  message("FOV QC was enabled for this run. Importing FOV QC metrics table.")
  qc_fov <- read_qc_tables("_fovQCmetrics\\.csv$")
  if (anyDuplicated(qc_fov[, c("sample_id", "FOV")])) {
    stop("Duplicate sample_id x FOV combinations")
  }
}

# Import spatialhub and potential user-defined, sample-level metadata table
message("Importing sample QC metrics table.")
qc_sample <- read_qc_tables("_sampleQCmetrics\\.csv$")

message("Importing user-defined sample-level metadata.")
sample_meta <- data.table::fread(opt$metadata_file)
sample_meta <- sample_meta |>
  dplyr::filter(sample_id %in% qc_sample$sample_id)

if (anyDuplicated(sample_meta$sample_id) ||
    anyDuplicated(qc_sample$sample_id)) {
  stop("Duplicate sample_id in sample-level metadata or QC metrics")
}


### Add instance-level QC flags

message("Adding instance-level QC flags.")
qc_instance <- qc_instance |>
  dplyr::mutate(
    qcFlagInstance_nCount = ifelse(
      !is.na(nCount_RNA) & nCount_RNA >= qc$instance_nCount_cutoff,
      "Pass", "Fail"
    ),
    qcFlagInstance_nFeature = ifelse(
      !is.na(nFeature_RNA) & nFeature_RNA >= qc$instance_nFeature_cutoff,
      "Pass", "Fail"
    ),
    qcFlagInstance_percentNeg = ifelse(
      !is.na(percentNegCounts) &
        percentNegCounts < qc$instance_percentNeg_cutoff,
      "Pass", "Fail"
    ),
    qcFlagInstance_complex = ifelse(
      !is.na(complexity) &
        complexity > qc$instance_complexity_cutoff,
      "Pass", "Fail"
    )
  ) |>
  dplyr::mutate(
    qcFlagInstance_summary = ifelse(
      qcFlagInstance_nCount == "Pass" &
        qcFlagInstance_nFeature == "Pass" &
        qcFlagInstance_percentNeg == "Pass" &
        qcFlagInstance_complex == "Pass",
      "Pass", "Fail"
    )
  )


### Add FOV-level QC flags

dfc <- qc_instance

if(opt$runFOVqc) {

  message("Adding FOV-level QC flags.")
  qc_fov <- qc_fov |>
    dplyr::mutate(
      qcFlagFOV_nCount = ifelse(
        !is.na(nCount_RNA_FOVavg) &
          nCount_RNA_FOVavg >= qc$fov_nCount_cutoff,
        "Pass", "Fail"
      ),
      qcFlagFOV_SNR = ifelse(
        !is.na(FOV_SNR) & FOV_SNR > qc$fov_SNR_cutoff,
        "Pass", "Fail"
      )
    )

  # Summarise FOV QC flags into a single summary flag
  fov_pass <- (
    qc_fov$qcFlagFOV_nCount == "Pass" &
      qc_fov$qcFlagFOV_SNR == "Pass"
  )
  if (opt$runBrukerFOVqc) {
    if (!"qcFlagGeneBias" %in% names(qc_fov)) {
      stop("Bruker FOV QC enabled but qcFlagGeneBias is missing")
    }
    fov_pass <- fov_pass &
      !is.na(qc_fov$qcFlagGeneBias) &
      qc_fov$qcFlagGeneBias == "Pass"
  }
  qc_fov$qcFlagFOV_summary <- ifelse(fov_pass, "Pass", "Fail")

  message("Writing FOV QC metrics table to file.")
  data.table::fwrite(
    qc_fov,
    file.path(opt$outdir, paste0(opt$run_name, "_fovQCmetrics.csv"))
  )

  # Append FOV QC flags to instance-level QC table (necessary for sample-level QC)
  qc_instance <- dplyr::left_join(
    qc_instance,
    qc_fov,
    by = c("sample_id", "FOV")
  )

  if (anyNA(qc_instance$qcFlagFOV_summary)) {
    stop("Some instances could not be matched to FOV QC results")
  }

}


### Add sample-level QC flags

message("Adding sample-level QC flags.")
sample_summary <- qc_instance |>
  dplyr::group_by(sample_id) |>
  dplyr::summarise(
    nInstance_sampleSum_passQC =
      sum(qcFlagInstance_summary == "Pass"),
    .groups = "drop"
  )

if (opt$runFOVqc) {
  fov_summary <- qc_instance |>
    dplyr::group_by(sample_id) |>
    dplyr::summarise(
      nInstancexnFOV_sampleSum_passQC =
        sum(qcFlagInstance_summary == "Pass" & qcFlagFOV_summary == "Pass"),
      samplePercentInstancesInFailedFOV =
        mean(qcFlagFOV_summary == "Fail") * 100,  # mean of logical vector gives proportion
      .groups = "drop"
  )

  sample_summary <- dplyr::left_join(
    sample_summary, fov_summary, by = "sample_id"
  )
}

qc_sample <- dplyr::left_join(
  qc_sample, sample_summary, by = "sample_id"
)


# Define rule for pass/fail for each sample

# Choose pass count according to the configured FOV policy
if (opt$runFOVqc && qc$exclInstancesInFailedFOVs) {
  passing <- qc_sample$nInstancexnFOV_sampleSum_passQC
} else {
  passing <- qc_sample$nInstance_sampleSum_passQC
}

# Derive percentage accordingly
qc_sample$samplePercentInstancesPassQC <- ifelse(
  qc_sample$nInstance_sampleSum > 0,
  100 * passing / qc_sample$nInstance_sampleSum,
  NA_real_
)

qc_sample$qcFlagSample_summary <- ifelse(
  !is.na(qc_sample$samplePercentInstancesPassQC) &
    passing >= qc$sampleMinInstances &
    qc_sample$samplePercentInstancesPassQC > qc$sampleMinPercent,
  "Pass", "Fail"
)

data.table::fwrite(
  qc_sample,
  file.path(opt$outdir, paste0(opt$run_name, "_sampleQCmetrics.csv"))
)


### Writing final QC pass/fail decision rule

# REMINDER: We don't filter FOVs altogether to avoid patchy samples.
# Poor quality FOVs are excluded at the sample level.
qc_instance <- dplyr::left_join(
  qc_instance,
  qc_sample,
  by = "sample_id"
)
qc_instance$qcFlagOverall_summary <- ifelse(
  qc_instance$qcFlagInstance_summary == "Pass" &
    qc_instance$qcFlagSample_summary == "Pass",
  "Pass", "Fail"
)

message("Writing *enriched* instance QC metrics table to file.")
data.table::fwrite(
  qc_instance,
  file.path(opt$outdir, paste0(opt$run_name, "_instanceQCmetrics.csv"))
)




# ----------

print("Done flagging QC failures within the following R environment:")
sessionInfo()