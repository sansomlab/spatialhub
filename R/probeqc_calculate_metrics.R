#!/usr/bin/Rscript --vanilla

# This script requires the following input from a CosMx experiment (QC variables may vary for other technologies):
#   - a *raw* gene x cell counts matrix
#   - FOV and sample-specific metadata
#   - *coordinates* and *FOV* metadata for each cell (for spatial plots in QC reports)
#   - panel-specific metadata, if using custom rather than universal panel

message("Calculating probe detection QC metrics for the provided cellxgene count matrix and associated metadata")
timestamp()


# ---------- Libraries ----------

stopifnot(
  require(optparse),
  require(tidyverse),
  require(data.table),
  require(zellkonverter),
  require(SingleCellExperiment)
)


# ---------- Options ----------

option_list <- list(
    make_option(
        c("--path2sce"),
        help="path to the counts matrix (AnnData.h5ad or SingleCellExperiment.rds) to use for calculation of cell-level QC metrics"
    ),
    make_option(
        c("--fov2sample"),
        help="path to the spatialhub samples TSV table describing samples covered in counts matrix"
    ),
    make_option(
        c("--background"),
        default="Negative",
        help="grep pattern for negative control (background) probe names"
    ),
    make_option(
        c("--poscontrol"),
        default="System",
        help="grep pattern for positive control (system) probe names"
    ),
    make_option(
        c("--customPanel"),
        default = "",
        help="path to description of custom panel used in the study, if applicable"
    ),
    make_option(
        c("--runFOVqc"),
        default = TRUE,
        help="whether to calculate FOV-level QC metrics"
    ),
    make_option(
        c("--FOV_nCount_cutoff"),
        default = 50,
        help="minimum FOV average nCount_RNA per cell for FOV to pass QC"
    ),
    make_option(
        c("--FOV_SNR_cutoff"),
        default = 3,
        help="minimum FOV signal to noise ratio for FOV to pass QC"
    ),
    make_option(
        c("--runBrukerFOVqc"),
        default = FALSE,
        help="whether to run Bruker's FOV QC tool"
    ),
    make_option(
        c("--brukerCodeFile"),
        default = "",
        help="path to Bruker FOV QC R script to source"
    ),
    make_option(
        c("--brukerDataFile"),
        default = "",
        help="path to Bruker RDS file containing panel-specific probes"
    ),
    make_option(
        c("--cell_nCount_cutoff"),
        default = 20,
        help=""
    ),
    make_option(
        c("--cell_nFeature_cutoff"),
        default = 5,
        help=""
    ),
    make_option(
        c("--cell_percentNeg_cutoff"),
        default = 5,
        help=""
    ),
    make_option(
        c("--cell_complexity_cutoff"),
        default = 1,
        help=""
    ),
    make_option(
        c("--sampleMinCells"),
        default = 2000,
        help="minimum number of cells from a sample that must pass QC for that sample to be retained"
    ),
    make_option(
        c("--sampleMinPercent"),
        default = 50,
        help="minimum percent of cells from a sample that must pass QC for that sample to be retained"
    ),
    make_option(
        c("--rmCellsInFailedFOVs"),
        default = FALSE,
        help="exclude cells located in a low-quality FOV, regardless of the individual quality metrics for this cell"
    )
)

opt <- parse_args(OptionParser(option_list=option_list))
set.seed(123456789)

cat("Running with options:\n")
print(opt)


# ---------- Setup ----------

# Read in (spatialhub) metadata file
path2meta <- opt$fov2sample
df <- read.delim(path2meta, sep = "\t")

if (!("sample_name" %in% names(df))) {
    print("'sample_name' not found: setting it to 'sample_id'")
    df$sample_name <- df$sample_id
}

# Convert to long format for merging with cell-level metadata
fov2sample <- df |>
  tidyr::separate_rows(fov_sequence, sep = ",") |>
  dplyr::rename(fov = fov_sequence) |>
  dplyr::filter(fov != "blank")
fov2sample$fov <- as.integer(as.character(fov2sample$fov))
fov2sample <- fov2sample |> dplyr::select(-fov_width, -fov_height)


### Read in `SingleCellExperiment` object on which to perform QC

stopifnot("file not found" = file.exists(opt$path2sce))
projName <- sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(opt$path2sce))
print(paste0("Working on project ", projName))

v <- unlist(strsplit(basename(opt$path2sce), "\\."))
file_extension <- tolower(v[length(v)])

if (file_extension == "h5ad") {
  sce <- zellkonverter::readH5AD(opt$path2sce, verbose = TRUE)
} else if (file_extension == "rds") {
  sce <- readRDS(opt$path2sce)
} else {
  stop("Please provide a `.h5ad` or `.rds` file as input.")
}


### Check required variables are saved in sce metadata, and add in sample-level ones if missing

stopifnot("'slide_id' not found in metadata" = "slide_id" %in% names(colData(sce)))

if (!("fov" %in% names(colData(sce)))) {  # e.g. if using a mask other than AtoMx default

  # Assert that there are not conflicts with the YAML options regarding FOV QC
  stopifnot("cannot run FOV QC: 'fov' variable not found in metadata" = opt$runFOVqc == FALSE)
  
  sample_key = ("sample_id" %in% names(colData(sce))) | ("sample_name" %in% names(colData(sce)))
  stopifnot("at least one of 'fov', 'sample_id' or 'sample_name' must appear in metadata" = sample_key)

  # Remove FOV information from sample metadata
  fov2sample <- fov2sample |> dplyr::select(-fov)
  fov2sample <- fov2sample[!duplicated(fov2sample), ]

}

# Append sample-level metadata

sce_meta <- as.data.frame(colData(sce))
sce_meta$cell_index <- colnames(sce)
sce_meta <- sce_meta |> dplyr::select(cell_index, everything())

matching_vars <- intersect(names(sce_meta), names(fov2sample))
spx_meta <- plyr::join(sce_meta, fov2sample, by = matching_vars, 
                       type = 'left', match = 'first')
rownames(spx_meta) <- spx_meta$cell_index

if (all(rownames(spx_meta) == rownames(colData(sce)))) {
  colData(sce) <- DataFrame(spx_meta)
}


### Extract counts matrix

if ("counts" %in% names(sce@assays@data)) {
  print("Assuming raw counts are stored in 'counts' slot")
  counts_mat <- sce@assays@data$counts
} else if ("X" %in% names(sce@assays@data)) {
  print("No 'counts' slot found. Assuming raw counts are stored in 'X' slot")
  counts_mat <- sce@assays@data$X
} else {
  stop("Raw counts not found. Please store as 'counts' or 'X' slot in sce@assays@data")
}

# Explicitely set matrix row and column names (essential when using 'reader = "R"' with zellkonverter)
rownames(counts_mat) <- rownames(sce)
colnames(counts_mat) <- colnames(sce)


### Handle optional arguments

if (opt$runFOVqc) {

  # Split sce object by slide and sample
  sce_ls <- list()

  colData(sce)$sample_key <- factor(paste0(as.character(colData(sce)$slide_id), "_",
                                           as.character(colData(sce)$sample_name)))

  for (sample in levels(colData(sce)$sample_key)) {
    
    keepIdx <- colData(sce)$cell_index[colData(sce)$sample_key == sample]
    sce_subset <- sce[, keepIdx]
    sce_ls[[sample]] <- sce_subset
    
  }
  
}

# OPTIONAL: Import custom probe panel, if applicable
if (nchar(opt$customPanel) > 0) {
  custom_probes <- read.csv(opt$customPanel)
  #custom_probes$group <- factor(as.character(custom_probes$species))  # To be refined
}



# ---------- Task 1: Calculate cell-level QC metrics ----------

# This set of metrics is independent of 'slide_id', and can thus be re-calculated on the dataset as a whole
print("Calculating cell-level QC metrics.")


## 1A: Re-calculate nCount/nFeature per cell from counts matrix

print(paste0("Total number of probes (including control and negative): ", dim(counts_mat)[1]))  
    # check this matches expectations for the corresponding panel
allc <- colSums(counts_mat) |> data.frame(); names(allc) <- "nCount_all"
allf <- colSums(counts_mat > 0) |> data.frame(); names(allf) <- "nFeature_all"


# negative (background) probes only
counts_neg <- counts_mat[grep(opt$background, rownames(counts_mat)), ]
print(paste0("Number of negative control probes: ", dim(counts_neg)[1]))
rownames(counts_neg)
negc <- colSums(counts_neg) |> data.frame(); names(negc) <- "nCount_neg"
n_neg <- nrow(counts_neg)  # storing for later

# system/positive (control) probes only
counts_sys <- counts_mat[grep(opt$poscontrol, rownames(counts_mat)), ]
print(paste0("Number of positive (system) control probes: ", dim(counts_sys)[1]))
rownames(counts_sys) |> head(n = 30)
sysc <- colSums(counts_sys) |> data.frame(); names(sysc) <- "nCount_sys"


# list of probes that do not belong to the universal panel
non_ucc_probes <- c(rownames(counts_neg), rownames(counts_sys))

if (nchar(opt$customPanel) > 0) {
  
  # custom probes (add-on panel) only
  counts_add <- counts_mat[rownames(counts_mat) %in% custom_probes$probe_name, ]
  print(paste0("Number of additional probes in custom panel: ", dim(counts_add)[1]))
  print(rownames(counts_add) |> head(n = 30))
  addc <- colSums(counts_add) |> data.frame(); names(addc) <- "nCount_add"
  addf <- colSums(counts_add > 0) |> data.frame(); names(addf) <- "nFeature_add"
  
  non_ucc_probes <- c(non_ucc_probes, rownames(counts_add))
  
}

# standard probes (universal panel) only
counts_ucc <- counts_mat[!(rownames(counts_mat) %in% non_ucc_probes), ]
default_probes <- rownames(counts_mat)
print(paste0("Number of probes in universal characterization panel: ", dim(counts_ucc)[1]))
uccc <- colSums(counts_ucc) |> data.frame(); names(uccc) <- "nCount_ucc"
uccf <- colSums(counts_ucc > 0) |> data.frame(); names(uccf) <- "nFeature_ucc"


# Bringing it all into one data frame
df <- cbind.data.frame(negc, sysc, uccc, uccf)
if (nchar(opt$customPanel) > 0) { df <- cbind.data.frame(df, addc, addf) }

# Derive total counts back again for custom panels
if (nchar(opt$customPanel) > 0) {
  df$nCount_RNA <- df$nCount_add + df$nCount_ucc
  df$nFeature_RNA <- df$nFeature_add + df$nFeature_ucc
  n_probes <- nrow(counts_ucc) + nrow(counts_add)
} else {
  df$nCount_RNA <- df$nCount_ucc
  df$nFeature_RNA <- df$nFeature_ucc
  df <- df |> dplyr::select(-nCount_ucc, -nFeature_ucc)
  n_probes <- nrow(counts_ucc)
}


## 1B: Calculate derived variables (as described in AtoMx user manual)

# Proportion of negative (does not consider system control probes)
df$percentNegCounts <- (df$nCount_neg / (df$nCount_RNA + df$nCount_neg)) * 100 

# "complexity" (transcripts per detected probe)
df$complexity <- df$nCount_RNA / df$nFeature_RNA

# Compute nCount for a random subset of probes (can help to spot FOVs with lower signal)
if (all(rownames(df) == colnames(counts_mat))) {
  df$nCount_rand20 <- colSums(counts_mat[sample(default_probes)[1:20], ])
} else {
  stop("rownames mismatch!")
}


## 1C: Add QC flags based on user-defined thresholds

df$qcFlagCell_nCount <- ifelse(df$nCount_RNA >= opt$cell_nCount_cutoff, "Pass", "Fail")
df$qcFlagCell_nFeature <- ifelse(df$nFeature_RNA >= opt$cell_nFeature_cutoff, "Pass", "Fail")
df$qcFlagCell_percentNeg <- ifelse(df$percentNegCounts < opt$cell_percentNeg_cutoff, "Pass", "Fail")
df$qcFlagCell_complex <- ifelse(df$complexity > opt$cell_complexity_cutoff, "Pass", "Fail")

cell_pass_qc <- df$nCount_RNA >= opt$cell_nCount_cutoff &
  df$nFeature_RNA >= opt$cell_nFeature_cutoff &
  df$percentNegCounts < opt$cell_percentNeg_cutoff &
  df$complexity > opt$cell_complexity_cutoff
df$qcFlagCell_summary <- ifelse(cell_pass_qc, "Pass", "Fail")


## 1D: Fetch back pre-existing cell-level metadata which is relevant for QC

df <- cbind("cell_index" = rownames(df), df)
dfc <- plyr::join(df, spx_meta, by = "cell_index")
if ("unassignedTranscripts" %in% names(spx_meta)) {
  dfc <- dfc |> dplyr::select(-unassignedTranscripts)  # This AtoMx variable is relevant at the FOV level
}
df0 <- dfc[, names(dfc) %in% c(names(df), names(fov2sample))]

# Re-order columns for more intuitive browsing
#names(dfc)[names(dfc) %in% names(fov2sample)]
#dfc <- dfc |> dplyr::select(TBC, everything())


## 1D: Save cell-level QC metrics

print("Saving cell-level QC metrics:")
head(dfc)
write.csv(dfc, row.names = FALSE, quote = FALSE, 
          file = paste0("probeqc.dir/", projName, "_cellQCmetrics.csv"))




# ---------- Task 2: Calculate FOV-level QC metrics ----------

# *WARNING*: An FOV is defined within the context of its 'slide_id' 
# => make sure the fov index takes the source slide into account

if (opt$runFOVqc) {

  print("Calculating FOV-level QC metrics.")

  ## 2A: Aggregate cell-level metrics at the FOV level

  df <- df0 |> dplyr::select(-cell_index)
  df$fov_index <- factor(paste0("FOV", df$fov, "_", df$slide_id))

  # focusing on quantitative variables that make sense to aggregate at the FOV level
  # NOTE: could be interesting to have these FOV-level metrics for fluorescence data too...
  df <- df |> dplyr::select(fov_index, 
                            nCount_RNA, nFeature_RNA, complexity,
                            nCount_neg, nCount_sys)
  df$nCell <- 1

  df_sum <- aggregate(. ~ fov_index, df, sum)
  df_sum <- df_sum[, -which(names(df_sum) %in% c("nFeature_RNA", "complexity"))]
  names(df_sum) <- paste0(names(df_sum), "_FOVsum"); names(df_sum)[1] <- "fov_index"

  df_avg <- aggregate(. ~ fov_index, df, mean)
  df_avg <- df_avg[, -which(names(df_avg) == "nCell")]
  names(df_avg) <- paste0(names(df_avg), "_FOVavg"); names(df_avg)[1] <- "fov_index"

  dfov <- plyr::join(df_sum, df_avg, by = "fov_index")


  ## 2B: Calculate derived variables (as described in AtoMx user manual)

  # Signal to noise ratio
  # defined as the FOV total number of detected probes normalized to the number of distinct probes in the panel, 
  # divided by the FOV total number negative probes normalized to the number of distinc negative probes.
  dfov$FOV_SNR <- (dfov$nCount_RNA_FOVsum / n_probes) / (dfov$nCount_neg_FOVsum / n_neg)


  # Fetch back pre-existing FOV-level metadata which is relevant for QC
  keepVar <- sort(which(names(spx_meta) %in% c("unassignedTranscripts", names(fov2sample))))
  spx_meta_fov <- spx_meta[, keepVar]

  # De-duplicating dataset at FOV-level, and adding fov_index
  spx_meta_fov <- spx_meta_fov[!duplicated(spx_meta_fov), ]
  fov_index <- paste0("FOV", spx_meta_fov$fov, "_", spx_meta_fov$slide_id)
  spx_meta_fov <- cbind("fov_index" = fov_index, spx_meta_fov)
  dfov <- plyr::join(dfov, spx_meta_fov, by = "fov_index")


  ## 2C: Run AtoMX FOV QC tool

  if (opt$runBrukerFOVqc) {
    print("Running Bruker FOV QC tool.")

    # Initiate relevant variables
    dfov$sample_key <- paste0(as.character(dfov$slide_id), "_", as.character(dfov$sample_name))
    #dfov$brukerFOVqc <- NA
    dfov$qcFlagInstr <- dfov$qcFlagGeneBias <- NA
    
    # Source functions and barcodes from Bruker Spatial Biology: https://github.com/Nanostring-Biostats/CosMx-Analysis-Scratch-Space/tree/Main/_code/FOV%20QC
    source(opt$brukerCode)
    all_panels <- readRDS(opt$brukerData)

    # Select barcodes corresponding to ucc panel used in the study
    ovlp <- c()
    for (i in 1:length(all_panels)) {
      v <- sum(default_probes %in% all_panels[[i]]$gene) / length(all_panels[[i]]$gene)
      ovlp <- c(ovlp, v)
    }
    stopifnot("Limited overlap between AtoMx panel and provided default probes. Check probe names were correctly imported." = max(ovlp) > 0.8)
    idx <- which(ovlp == max(ovlp))
    barcodes <- all_panels[[idx]]


    # Run FOV QC tool on sce object
    res_fov_qc <- list()
    for (sample in names(sce_ls)) {
      
      if (length(grep("_NA$", sample)) == 0) {

        print(paste0("Checking instrument failures for sample ", sample))
        res_fov_qc[[sample]] <- runFOVQC(counts = t(assay(sce_ls[[sample]])), 
                                         xy = data.frame(colData(sce_ls[[sample]])[, c("CenterX_global_px", "CenterY_global_px")]),
                                         fov = colData(sce_ls[[sample]])[, "fov"], 
                                         barcodemap = barcodes,  # REMINDER: this only considers probes in the custom panel (including Negative)
                                         max_prop_loss = 0.6, max_totalcounts_loss = 0.6)  # default 0.6, the higher, the more relaxed the QC
        
        # Append to QC metrics dataset for the corresponding slide x sample
        failedFOVs <- res_fov_qc[[sample]]$flaggedfovs; biasFOVs <- res_fov_qc[[sample]]$flaggedfovs_forbias
        dfov$qcFlagInstr[dfov$sample_key == sample] <- ifelse(dfov$fov[dfov$sample_key == sample] %in% failedFOVs, "Fail", "Pass")
        dfov$qcFlagGeneBias[dfov$sample_key == sample] <- ifelse(dfov$fov[dfov$sample_key == sample] %in% biasFOVs, "Fail", "Pass")
        #if ("sample_qc" %in% names(res_fov_qc[[sample]])) {
        #  dfov$brukerFOVqc <- res_fov_qc[[sample]]$sample_qc
        #}

      } else {
        unassigned_fovs <- unique(colData(sce_ls[[sample]])[, "fov"])
        paste(unassigned_fovs, collapse = ", ")
        print(paste0("No sample assigned to the following FOV(s): ", unassigned_fovs))
        print("Skipping instrument failure QC for these.")
      }

    }

    # Extract and save list of affected genes
    ls_bias <- list()
    for (sample in names(res_fov_qc)) {
      
      df <- res_fov_qc[[sample]]$flagged_fov_x_gene |> data.frame()
      
      if (nrow(df) > 0) {
        
        df$count <- 1
        df <- aggregate.data.frame(df, . ~ fov + gene, FUN = "sum")
        names(df) <- c("fov", "gene", "failed_cycles")
        df$sample_index <- sample
        
        df <- df |>
          dplyr::arrange(gene) |>
          dplyr::arrange(desc(failed_cycles)) |>
          dplyr::arrange(fov)
        
        ls_bias[[sample]] <- df
        
      }
      
    }
    df_bias <- do.call(rbind.data.frame, ls_bias)
    write.csv(df_bias, row.names = FALSE, quote = FALSE, 
              file = paste0("probeqc.dir/", projName, "_instrGeneBias.csv"))

    # Removing temporary sample_key
    dfov <- dfov |> dplyr::select(-sample_key)

  }


  ## 2D: Add QC flags based on user-defined thresholds

  dfov$qcFlagFOV_nCount <- ifelse(dfov$nCount_RNA_FOVavg >= opt$FOV_nCount_cutoff, "Pass", "Fail")
  dfov$qcFlagFOV_SNR <- ifelse(dfov$FOV_SNR > opt$FOV_SNR_cutoff, "Pass", "Fail")

  if (opt$runBrukerFOVqc) {
    FOV_pass_qc <- dfov$nCount_RNA_FOVavg >= opt$FOV_nCount_cutoff &
      dfov$FOV_SNR > opt$FOV_SNR_cutoff &
      dfov$qcFlagGeneBias == "Pass"
  } else {
    FOV_pass_qc <- dfov$nCount_RNA_FOVavg >= opt$FOV_nCount_cutoff &
      dfov$FOV_SNR > opt$FOV_SNR_cutoff
  }
  dfov$qcFlagFOV_summary <- ifelse(FOV_pass_qc, "Pass", "Fail")


  ## 2E: Save FOV-level QC metrics

  print("Saving FOV-level QC metrics:")
  print(head(dfov))
  write.csv(dfov, row.names = FALSE, quote = FALSE, 
            file = paste0("probeqc.dir/", projName, "_fovQCmetrics.csv"))

}




# ---------- Task 3: Calculate sample-level QC metrics ----------

# *WARNING*: A sample_*name* may be defined within the context of its 'slide_id' 
# => make sure the sample index takes the source slide into account

print("Calculating sample-level QC metrics.")


## 3A: Aggregate cell-level metrics at the sample level

df <- df0 |> dplyr::select(-cell_index)
df$sample_index <- factor(paste0(df$sample_name, "_", df$slide_id))

# focusing on quantitative variables that make sense to aggregate at the sample level
df <- df |> dplyr::select(sample_index,
                          nCount_RNA, nFeature_RNA, complexity,
                          nCount_neg, nCount_sys)
df$nCell <- 1

df_sum <- aggregate(. ~ sample_index, df, sum)
df_sum <- df_sum[, -which(names(df_sum) %in% c("complexity", "nFeature_RNA"))]
names(df_sum) <- paste0(names(df_sum), "_sampleSum"); names(df_sum)[1] <- "sample_index"

df_avg <- aggregate(. ~ sample_index, df, mean)
df_avg <- df_avg[, -which(names(df_avg) == "nCell")]
names(df_avg) <- paste0(names(df_avg), "_sampleAvg"); names(df_avg)[1] <- "sample_index"

dfs <- plyr::join(df_sum, df_avg, by = "sample_index")


## 3B: Additional sample-level variables, to check for group/batch effects

# Fetch back pre-existing FOV-level metadata which is relevant for QC
keepVar <- sort(which(names(spx_meta) %in% names(fov2sample)))
spx_meta_sample <- spx_meta[, keepVar]

# De-duplicating dataset at sample-level, and adding sample_index
spx_meta_sample <- spx_meta_sample[!duplicated(spx_meta_sample), ]
sample_index <- paste0(spx_meta_sample$sample_name, "_", spx_meta_sample$slide_id)
spx_meta_sample <- cbind("sample_index" = sample_index, spx_meta_sample)

dfs <- plyr::join(dfs, spx_meta_sample, by = "sample_index")
sample_qc_vars <- names(dfs)


## 3C: Add QC flags based on user-defined thresholds

# Note: This summary requires integration of QC metrics across all (available) levels
dfa <- plyr::join(dfc, dfs, by = intersect(names(dfc), names(dfs)))


# Tally the cells that passed *cell* QC for each sample
dfa_ls <- split(dfa, f = dfa$sample_index)

for (sample in names(dfa_ls)) { 
  df <- dfa_ls[[sample]]
    
  # In summary, how many cells for that sample pass QC?
  # i.e. how many cells have enough counts, etc.
  sample_cell_pass_qc <- df$qcFlagCell_summary == "Pass"
  df$nCell_sampleSum_passQC <- sum(sample_cell_pass_qc)

  dfa_ls[[sample]] <- df
}
dfa <- do.call(rbind.data.frame, dfa_ls)


# Tally the cells that pased *cell* AND *FOV* QC for each sample
if (opt$runFOVqc) {
  dfa <- plyr::join(dfa, dfov, by = intersect(names(dfa), names(dfov)))

  # Define sample-level FOV QC pass/fail
  dfa_ls <- split(dfa, f = dfa$sample_index)

  for (sample in names(dfa_ls)) {
    df <- dfa_ls[[sample]]
    
    s <- sum(df$qcFlagFOV_summary == "Fail")
    df$samplePercentCellsInFailedFOV <- round(s / nrow(df), 3) * 100
    #r <- length(unique(df$fov[df$qcFlagsFOV_summary == "Fail"]))
    #df$samplePercentFailedFOVs <- round(r / length(unique(df$fov)), 3) * 100
    
    # In summary, how many cells for that sample pass all QC?
    # i.e. how many cells have enough counts, etc. AND fall within a high-qaulity FOV
    sample_cell_pass_qc <- df$qcFlagCell_summary == "Pass" &
      df$qcFlagFOV_summary == "Pass"
    df$nCellxnFOV_sampleSum_passQC <- sum(sample_cell_pass_qc)

    dfa_ls[[sample]] <- df
  }
  dfa <- do.call(rbind.data.frame, dfa_ls)
}


if (opt$rmCellsInFailedFOVs) {

  # Define summary sample QC metric
  # i.e. sample with a minimum absolute count of cells and percent of cells passing QC
  # based on cell QC and FOV QC both
  dfa$samplePercentCellsPassQC <- round(dfa$nCellxnFOV_sampleSum_passQC/dfa$nCell_sampleSum, 3) * 100
  sample_pass_qc <- dfa$nCellxnFOV_sampleSum_passQC >= opt$sampleMinCells &
    dfa$samplePercentCellsPassQC > opt$sampleMinPercent
  dfa$qcFlagSample_summary <- ifelse(sample_pass_qc, "Pass", "Fail")

} else {

  # Define summary sample QC metric
  # i.e. sample with a minimum absolute count of cells and percent of cells passing QC
  # based on cell QC only

  dfa$samplePercentCellsPassQC <- round(dfa$nCell_sampleSum_passQC/dfa$nCell_sampleSum, 3) * 100
  sample_pass_qc <- dfa$nCell_sampleSum_passQC >= opt$sampleMinCells &
    dfa$samplePercentCellsPassQC > opt$sampleMinPercent
  dfa$qcFlagSample_summary <- ifelse(sample_pass_qc, "Pass", "Fail")

}


## 3D: Save sample-level QC metrics

# Reducing data frame back again to one row per sample
keepVar <- sort(which(names(dfa) %in% c(sample_qc_vars, 
                  "samplePercentCellsPassQC", "nCell_sampleSum_passQC", 
                  "samplePercentCellsInFailedFOV", "nCellxnFOV_sampleSum_passQC",
                  "qcFlagSample_summary")
                  )
                )
if ("fov" %in% keepVar) { keepVar <- keepVar[keepVar != "fov"] }
dfs <- dfa[, keepVar] #|> dplyr::select(-fov)
dfs <- dfs[!duplicated(dfs), ]
dfs <- dfs[!is.na(dfs$sample_id), ]
rownames(dfs) <- 1:nrow(dfs)

print("Saving sample-level QC metrics:")
head(dfs)
write.csv(dfs, row.names = FALSE, quote = FALSE, 
          file = paste0("probeqc.dir/", projName, "_sampleQCmetrics.csv"))

print("Saving exhaustive cell-level QC metrics (project metadata)")
rownames(dfa) <- 1:nrow(dfa); head(dfa)  # column order is a little random, but everything is there
write.csv(dfa, row.names = FALSE, quote = FALSE, 
          file = paste0("probeqc.dir/", projName, "_metadata_QC.csv"))




# ---------- Task 4: Find top expressed genes ----------

v <- rowSums(counts_mat) |> sort(decreasing = TRUE)
v10 <- v[-grep(opt$poscontrol, names(v))] |> head(n = 10)
print("Top most detected probes in this study (excluding positive controls): ")
v10

u <- rowSums(counts_mat) |> sort(decreasing = FALSE)
u10 <- u[-grep(opt$poscontrol, names(u))] |> head(n = 10)
print("Least detected probes in this study (excluding positive controls): ")
u10


print("Done calculating QC metrics within the following R environment:")
sessionInfo()
