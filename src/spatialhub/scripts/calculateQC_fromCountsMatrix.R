#!/usr/bin/Rscript --vanilla

# This script requires the following input from a spatial experiment
# (Note: This script was developed with CosMx data in mind - QC variables may vary slightly for other technologies):
#   - a *raw* gene x instance counts matrix (per sample - essential for FOV QC)
#   - regex pattern for negative and postivie control probes (default: "Negative" and "System")
#   - panel-specific metadata, if using custom rather than universal panel
#   - pointers to Bruker FOV QC tool, if applicable (CosMx technology only)

message("Calculating probe detection QC metrics")
timestamp()


# ---------- Libraries ----------

stopifnot(
  require(optparse),
  require(dplyr),
  require(tools),
  require(data.table),
  require(anndataR),
  require(SingleCellExperiment)
)


# ---------- Options ----------

option_list <- list(
    make_option(
        c("--path2sce"),
        help="path to raw H5AD counts matrix"
    ),
    make_option(
        c("--outdir"),
        required = TRUE,
        help = "Directory in which to write ProbeQC outputs"
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
    )
)

opt <- parse_args(OptionParser(option_list=option_list))
cat("Running with options:\n")
print(opt)


# --------- Setup ----------

### Read in `SingleCellExperiment` object on which to perform QC

stopifnot("file not found" = file.exists(opt$path2sce))
sample_id <- pathStem <- tools::file_path_sans_ext(basename(opt$path2sce))
print(paste0("Working on matrix ", pathStem))
sce <- anndataR::read_h5ad(opt$path2sce, as = "SingleCellExperiment")

# Use sample_id from metadata when available; otherwise stick to stem filename.
if ("sample_id" %in% names(colData(sce))) {
  s <- unique(as.character(colData(sce)$sample_id))
  if (length(s) != 1 || is.na(s)) {
    stop("Expected one non-missing sample_id per H5AD")
  }
  sample_id <- s
}


### Extract counts matrix

if ("counts" %in% assayNames(sce)) {
  print("Assuming raw counts are stored in 'counts' slot")
  counts_mat <- assay(sce, "counts")
} else if ("X" %in% assayNames(sce)) {
  print("No 'counts' slot found. Assuming raw counts are stored in 'X' slot")
  counts_mat <- assay(sce, "X")
} else {
  stop("Raw counts not found. Please store as 'counts' or 'X' slot in sce@assays@data")
}

# Ensure counts matrix and colData alignement are preserved
stopifnot(
  identical(colnames(counts_mat), colnames(sce)),
  !anyDuplicated(colnames(sce))
)
stopifnot(
  identical(rownames(counts_mat), rownames(sce)),
  !anyDuplicated(rownames(sce))
)


### OPTIONAL: Import custom probe panel, if applicable
if (nchar(opt$customPanel) > 0) {
  custom_probes <- read.csv(opt$customPanel)
}


### Initiate probe classifier

probe_class <- data.frame(probe_name = rownames(counts_mat), 
                          negative = grepl(opt$background, rownames(counts_mat)),
                          positive = grepl(opt$poscontrol, rownames(counts_mat)))
#tail(probe_class)




# ---------- Task 1: Calculate instance-level QC metrics ----------

# This set of metrics is independent of any metadata (sample, slide or FOV)
# and can thus be easily calculated on any chunk of the dataset
print("Calculating instance-level (segment or bin) QC metrics.")
print(paste0("Total number of probes (including control and negative): ", dim(counts_mat)[1]))  

# negative (background) probes only
counts_neg <- counts_mat[grep(opt$background, rownames(counts_mat)), , drop = FALSE]
print(paste0("Number of negative control probes: ", dim(counts_neg)[1]))
rownames(counts_neg)
negc <- data.frame(nCount_neg = colSums(counts_neg))

# system/positive (control) probes only
counts_sys <- counts_mat[grep(opt$poscontrol, rownames(counts_mat)), , drop = FALSE]
print(paste0("Number of positive (system) control probes: ", dim(counts_sys)[1]))
rownames(counts_sys) |> head(n = 30)
sysc <- data.frame(nCount_sys = colSums(counts_sys))


# define and derive counts for probes that do not belong to the universal panel
non_ucc_probes <- c(rownames(counts_neg), rownames(counts_sys))

if (nzchar(opt$customPanel) > 0) {
  
  # update probe classifier
  probe_class$custom <- (rownames(counts_mat) %in% custom_probes$probe_name)

  # custom probes (add-on panel) only
  counts_add <- counts_mat[rownames(counts_mat) %in% custom_probes$probe_name, , drop = FALSE]
  print(paste0("Number of additional probes in custom panel: ", dim(counts_add)[1]))
  print(rownames(counts_add) |> head(n = 30))
  addc <- data.frame(nCount_add = colSums(counts_add))
  addf <- data.frame(nFeature_add = colSums(counts_add > 0))
  
  non_ucc_probes <- c(non_ucc_probes, rownames(counts_add))
  
}
probe_class$ucc <- !as.logical(rowSums(probe_class[, 2:ncol(probe_class), drop = FALSE]))

# finally, derive counts for standard probes (universal panel) only
counts_ucc <- counts_mat[!(rownames(counts_mat) %in% non_ucc_probes), , drop = FALSE]
print(paste0("Number of probes in universal characterization panel: ", dim(counts_ucc)[1]))
uccc <- data.frame(nCount_ucc = colSums(counts_ucc))
uccf <- data.frame(nFeature_ucc = colSums(counts_ucc > 0))


# Bringing it all into one data frame
df <- cbind.data.frame('sample_id' = sample_id,
                       negc, sysc, uccc, uccf)
if (nchar(opt$customPanel) > 0) { df <- cbind.data.frame(df, addc, addf) }

# Derive total counts back again for custom panels
if (nzchar(opt$customPanel) > 0) {
  df$nCount_RNA <- df$nCount_add + df$nCount_ucc
  df$nFeature_RNA <- df$nFeature_add + df$nFeature_ucc
  n_probes <- nrow(counts_ucc) + nrow(counts_add)
} else {
  df$nCount_RNA <- df$nCount_ucc
  df$nFeature_RNA <- df$nFeature_ucc
  df <- df |> dplyr::select(-nCount_ucc, -nFeature_ucc)
  n_probes <- nrow(counts_ucc)
}

# Derived instance-level QC metrics
df$percentNegCounts <- ifelse(
  df$nCount_RNA + df$nCount_neg > 0,
  100 * df$nCount_neg / (df$nCount_RNA + df$nCount_neg),
  NA_real_
)
df$complexity <- ifelse(
  df$nFeature_RNA > 0,
  df$nCount_RNA / df$nFeature_RNA,
  NA_real_
)
# Important NOTE: the grid H5AD will contain occupied bins, while the complete QC metadata will also include empty imaged bins.
# Those additional empty bins will need their count metrics set to zero and their undefined ratios set to NA.


# Save instance-level QC metrics and metadata (avoids re-loading the entire dataset later)
print("Saving instance-level QC metrics:")
df <- df |>
  dplyr::mutate(instance_id = rownames(df),
                .after = sample_id)
head(df)
write.csv(df, row.names = FALSE, quote = FALSE, 
          file = paste0(opt$outdir, "/", pathStem, "_instanceQCmetrics.csv"))




# ---------- Task 2: Run AtoMX FOV QC tool ----------

# Unlike other FOV QC metrics (calculated elsewhere), which can be derived from instance-level QC metrics,
# Bruker's FOV QC tool is run using the counts matrix. To avoid reloading the entire dataset later,
# we will therefore run it here (if applicable).
 
if (opt$runBrukerFOVqc) {
  
  ### Pre-flight checks
  
  stopifnot(require(FNN))
  metadata <- as.data.frame(colData(sce))  
  
  # Retrieve FOV identifiers
  if ("FOV" %in% names(metadata)) {
    fov <- metadata$FOV
  } else if ("fov" %in% names(metadata)) {
    fov <- metadata$fov
  } else {
    ids <- if ("instance_id" %in% names(metadata)) {
      as.character(metadata$instance_id)
    } else { colnames(sce) }

    # check if 'instance_id' follows the AtoMx pattern and extract FOV from it
    if (!all(grepl("^c_[0-9]+_[0-9]+_[0-9]+$", ids))) {
      stop("Cannot derive FOV from observation identifiers")
    }
    fov <- data.table::tstrsplit(ids, "_")[[3]]
  }

  if (anyNA(fov)) {
    stop("Missing FOV identifiers in Bruker QC input")
  }

  metadata$FOV <- fov

  # Extract coordinates
  if (all(c("array_col", "array_row") %in% names(metadata))) {
    xy <- data.frame(
      x = metadata$array_col,
      y = metadata$array_row
    )
  } else if (all(c("x", "y") %in% names(metadata))) {
    xy <- metadata[, c("x", "y")]
  } else if (all(c("CenterX_global_px", "CenterY_global_px") %in% names(metadata))) {
    xy <- data.frame(
      x = metadata$CenterX_global_px,
      y = metadata$CenterY_global_px
    )
  } else {
    stop("Bruker FOV QC requires spatial coordinates")
  }
   
  # Source functions and barcodes from Bruker Spatial Biology:
  # https://github.com/Nanostring-Biostats/CosMx-Analysis-Scratch-Space/tree/Main/_code/FOV%20QC
  if (!file.exists(opt$brukerCodeFile)) {
    stop("Bruker code file not found: ", opt$brukerCodeFile)
  } else { source(opt$brukerCodeFile) }
  if (!file.exists(opt$brukerDataFile)) {
    stop("Bruker barcodes file not found: ", opt$brukerDataFile)
  } else { all_panels <- readRDS(opt$brukerDataFile) }
  
  # Select barcodes corresponding to universal (UCC) panel used in the study
  ovlp <- c()
  for (i in 1:length(all_panels)) {
    v <- sum(probe_class$probe_name[probe_class$ucc] %in% all_panels[[i]]$gene) / length(all_panels[[i]]$gene)
    ovlp <- c(ovlp, v)
  }
  stopifnot("Limited overlap between study panel and Bruker default probes. Check probe names were correctly imported." = max(ovlp) > 0.8)
  idx <- which(ovlp == max(ovlp))
  barcodes <- all_panels[[idx]]
  barcodes <- barcodes[barcodes$gene %in% probe_class$probe_name[probe_class$ucc], ]
  
  
  ### Ready to run!
  
  # Initialize a summary table where to save output
  dfov <- unique(data.frame(sample_id = sample_id, FOV = metadata$FOV))
  dfov$FOV <- as.numeric(as.character(dfov$FOV))
  dfov <- dfov |> dplyr::arrange(FOV)
  rownames(dfov) <- 1:nrow(dfov)
  
  print(paste0("Checking FOV instrument failures for sample ", sample_id))
  res_fov_qc <- runFOVQC(counts = t(counts_mat),
                         xy = xy, fov = metadata$FOV, 
                         barcodemap = barcodes,  # REMINDER: this only considers probes in the univeral panel (including Negative)
                         max_prop_loss = 0.6, max_totalcounts_loss = 0.6)  # default 0.6, the higher, the more relaxed the QC
  #summary(res_fov_qc)


  # Summarise flagged genes, retaining one row per FOV
  dfg <- as.data.frame(res_fov_qc$flagged_fov_x_gene)

  if (nrow(dfg) > 0) {
    dfg <- dfg |>
      dplyr::count(fov, gene, name = "cycles_gene_failed") |>
      dplyr::group_by(fov) |>
      dplyr::summarise(
        flagged_genes = paste(gene, collapse = ";"),
        nCycles_gene_failed = paste(cycles_gene_failed, collapse = ";"),
        total_flagged_genes = dplyr::n(),
        .groups = "drop"
      )

    dfg$FOV <- as.numeric(as.character(dfg$fov))
    dfg$fov <- NULL

    dfov <- dplyr::left_join(dfov, dfg, by = "FOV")
  } else {
    dfov$flagged_genes <- NA_character_
    dfov$nCycles_gene_failed <- NA_character_
    dfov$total_flagged_genes <- NA_integer_
  }

  # Evaluated FOVs with no flagged genes have zero failures
  dfov$total_flagged_genes[is.na(dfov$total_flagged_genes)] <- 0L

  # Ensure one row per FOV
  stopifnot(nrow(dfov) == nrow(dplyr::distinct(dfov, sample_id, FOV)))
  
  dfov$qcFlagInstr <- ifelse(dfov$FOV %in% res_fov_qc$flaggedfovs, "Fail", "Pass")
  dfov$qcFlagGeneBias <- ifelse(dfov$FOV %in% res_fov_qc$flaggedfovs_forbias, "Fail", "Pass")
    
  write.csv(dfov, row.names = FALSE, quote = FALSE, 
            file = paste0(opt$outdir, "/", pathStem, "_brukerQCresults.csv"))

  print(paste0("Bruker FOV QC completed for sample ", pathStem))
  
}




# ---------- Task 3: Probe QC metrics - Find top and least expressed genes ----------

if (all(rownames(counts_mat) == probe_class$probe_name)) {
  print("Adding summary statistics for each probe in the counts matrix to probe classifier table")
  probe_class <- probe_class |> 
    dplyr::mutate(nCells_detected = rowSums(counts_mat > 0)) |>
    dplyr::mutate(fractCells_detected = nCells_detected / ncol(counts_mat)) |>
    dplyr::mutate(nCount_total = rowSums(counts_mat)) |>
    dplyr::arrange(desc(nCount_total)) |>
    dplyr::mutate(nCount_rank = 1:nrow(probe_class)) |>
    dplyr::arrange(probe_name)
} else {
  stop("The probe classifier list and matrix row names do not match.")
}
write.csv(probe_class, row.names = FALSE, quote = FALSE, 
          file = paste0(opt$outdir, "/", pathStem, "_probeQCmetrics.csv"))


# Write instance-level metadata
metadata <- as.data.frame(colData(sce))
metadata$sample_id <- sample_id; metadata$instance_id <- colnames(sce)

# Retrieve FOV identifiers when available
if (!"FOV" %in% names(metadata)) {
  if ("fov" %in% names(metadata)) {
    metadata$FOV <- metadata$fov
  } else if (all(grepl(
    "^c_[0-9]+_[0-9]+_[0-9]+$",
    metadata$instance_id
  ))) {
    metadata$FOV <- data.table::tstrsplit(
      metadata$instance_id, "_"
    )[[3]]
  }
}

write.csv(metadata, row.names = FALSE, quote = FALSE, 
          file = paste0(opt$outdir, "/", pathStem, "_instanceMetadata.csv"))




# ----------

print("Done calculating QC metrics within the following R environment:")
sessionInfo()