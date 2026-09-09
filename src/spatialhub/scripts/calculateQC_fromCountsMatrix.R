#!/usr/bin/Rscript --vanilla

# This script requires the following input from a spatial experiment
# (Note: This script was developed with CosMx data in mind - QC variables may vary slightly for other technologies):
#   - a *raw* gene x cell counts matrix (per sample - essential for FOV QC)
#   - regex pattern for negative and postivie control probes (default: "Negative" and "System")
#   - panel-specific metadata, if using custom rather than universal panel
#   - pointers to Bruker FOV QC tool, if applicable (CosMx technology only)

message("Calculating probe detection QC metrics for the provided cellxgene count matrix and associated metadata")
timestamp()


# ---------- Libraries ----------

stopifnot(
  require(optparse),
  require(tidyverse),
  require(data.table),
  require(rhdf5),
  require(anndataR),
  require(SingleCellExperiment)
)


# ---------- Options ----------

option_list <- list(
    make_option(
        c("--path2sce"),
        help="path to the counts matrix (AnnData.h5ad or SingleCellExperiment.rds) to use for calculation of cell-level QC metrics"
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
set.seed(123456789)

cat("Running with options:\n")
print(opt)


# --------- Setup ----------

### Read in `SingleCellExperiment` object on which to perform QC

stopifnot("file not found" = file.exists(opt$path2sce))
pathStem <- sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(opt$path2sce))
print(paste0("Working on matrix ", pathStem))

v <- unlist(strsplit(basename(opt$path2sce), "\\."))
sample_id <- v[1]; file_extension <- tolower(v[length(v)])

if (file_extension == "h5ad") {
  sce <- anndataR::read_h5ad(opt$path2sce, as = "SingleCellExperiment")
} else if (file_extension == "rds") {
  sce <- readRDS(opt$path2sce)
} else {
  stop("Please provide a `.h5ad` or `.rds` file as input.")
}

# Check that relevant variables are present in sce metadata
stopifnot("'sample_id' not found in metadata" = 'sample_id' %in% names(colData(sce)))
stopifnot("'sample_id' not unique in metadata" = length(unique(colData(sce)$sample_id)) == 1)


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
#all(rownames(sce) == rownames(counts_mat)); all(colnames(sce) == colnames(counts_mat))


### OPTIONAL: Import custom probe panel, if applicable
if (nchar(opt$customPanel) > 0) {
  custom_probes <- read.csv(opt$customPanel)
  #custom_probes$group <- factor(as.character(custom_probes$species))  # To be refined
}


### Initiate probe classifier

probe_class <- data.frame(probe_name = rownames(counts_mat), 
                          negative = grepl(opt$background, rownames(counts_mat)),
                          positive = grepl(opt$poscontrol, rownames(counts_mat)))
#tail(probe_class)




# ---------- Task 1: Calculate cell-level QC metrics ----------

# This set of metrics is independent of any metadata (sample, slide or FOV)
# and can thus be easily calculated on any chunk of the dataset
print("Calculating cell-level QC metrics.")


# Re-calculate nCount/nFeature per cell from counts matrix
print(paste0("Total number of probes (including control and negative): ", dim(counts_mat)[1]))  

# negative (background) probes only
counts_neg <- counts_mat[grep(opt$background, rownames(counts_mat)), ]
print(paste0("Number of negative control probes: ", dim(counts_neg)[1]))
rownames(counts_neg)
negc <- colSums(counts_neg) |> data.frame(); names(negc) <- "nCount_neg"

# system/positive (control) probes only
counts_sys <- counts_mat[grep(opt$poscontrol, rownames(counts_mat)), ]
print(paste0("Number of positive (system) control probes: ", dim(counts_sys)[1]))
rownames(counts_sys) |> head(n = 30)
sysc <- colSums(counts_sys) |> data.frame(); names(sysc) <- "nCount_sys"


# define and derive counts for probes that do not belong to the universal panel
non_ucc_probes <- c(rownames(counts_neg), rownames(counts_sys))

if (nchar(opt$customPanel) > 0) {
  
  # update probe classifier
  probe_class$custom <- (rownames(counts_mat) %in% custom_probes$probe_name)

  # custom probes (add-on panel) only
  counts_add <- counts_mat[rownames(counts_mat) %in% custom_probes$probe_name, ]
  print(paste0("Number of additional probes in custom panel: ", dim(counts_add)[1]))
  print(rownames(counts_add) |> head(n = 30))
  addc <- colSums(counts_add) |> data.frame(); names(addc) <- "nCount_add"
  addf <- colSums(counts_add > 0) |> data.frame(); names(addf) <- "nFeature_add"
  
  non_ucc_probes <- c(non_ucc_probes, rownames(counts_add))
  
}
probe_class$ucc <- !as.logical(rowSums(probe_class[, 2:ncol(probe_class)]))

# finally, derive counts for standard probes (universal panel) only
counts_ucc <- counts_mat[!(rownames(counts_mat) %in% non_ucc_probes), ]
print(paste0("Number of probes in universal characterization panel: ", dim(counts_ucc)[1]))
uccc <- colSums(counts_ucc) |> data.frame(); names(uccc) <- "nCount_ucc"
uccf <- colSums(counts_ucc > 0) |> data.frame(); names(uccf) <- "nFeature_ucc"


# Bringing it all into one data frame
df <- cbind.data.frame('sample_id' = unique(colData(sce)$sample_id),
                       negc, sysc, uccc, uccf)
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

# Compute nCount for a random subset of probes (can help to spot FOVs with lower signal)
if (all(rownames(df) == colnames(counts_mat))) {
  df$nCount_rand20 <- colSums(counts_mat[sample(rownames(counts_mat))[1:20], ])
} else {
  stop("rownames mismatch!")
}

# NOTE: although the next 2 variables could be derived from previous summary metrics
# (without re-loading the counts matrix), we will calculate them in this script to streamline outputs
# (i.e. to gather all *cell* QC metrics in one CSV file)

# Proportion of negative (does not consider system control probes)
df$percentNegCounts <- (df$nCount_neg / (df$nCount_RNA + df$nCount_neg)) * 100 

# "complexity" (transcripts per detected probe)
df$complexity <- df$nCount_RNA / df$nFeature_RNA




# Save cell-level QC metrics and metadata (avoids re-loading the entire dataset later)
print("Saving cell-level QC metrics:")
df <- df |>
  dplyr::mutate(instance_id = rownames(df),
                .after = sample_id)
head(df)
write.csv(df, row.names = FALSE, quote = FALSE, 
          file = paste0(opt$outdir, "/", pathStem, "_cellQCmetrics.csv"))
write.csv(as.data.frame(colData(sce)), row.names = FALSE, quote = FALSE, 
          file = paste0(opt$outdir, "/", pathStem, "_cellMetadata.csv"))
write.csv(probe_class, row.names = FALSE, quote = FALSE, 
          file = paste0(opt$outdir, "/", pathStem, "_probeClassifier.csv"))




# ---------- Task 2: Run AtoMX FOV QC tool ----------

# Unlike other FOV QC metrics (calculated elsewhere), which can be derived from cell-level QC metrics,
# Bruker's FOV QC tool is run using the counts matrix. To avoid reloading the entire dataset later,
# we will therefore run it here (if applicable).
 
if (opt$runBrukerFOVqc) {
  
  ### Pre-flight checks
  
  stopifnot(require(FNN))
  
  if (!("x" %in% names(colData(sce))) | !("y" %in% names(colData(sce)))) {
    
    # If no x/y coordinate found in metadata
    if ("array_row" %in% names(colData(sce))) {
      # => check for synonyms for H5AD derived using SpatialHub
      col_idx <- which(names(colData(sce)) == "array_row")
      names(colData(sce))[col_idx] <- "x"
      col_idy <- which(names(colData(sce)) == "array_col")
      names(colData(sce))[col_idy] <- "y"
    } else if ("CenterX_global_px" %in% names(colData(sce))) {
      # => check for synonyms for H5AD derived directly from AtoMx flat files
      col_idx <- which(names(colData(sce)) == "CenterX_global_px")
      names(colData(sce))[col_idx] <- "x"
      col_idy <- which(names(colData(sce)) == "CenterY_global_px")
      names(colData(sce))[col_idy] <- "y"
    } else {
      stop("cannot run Bruker FOV QC: x/y coordinates not found in metadata")
    }
  }
  
  if (!("FOV" %in% names(colData(sce)))) {
    
    # If no 'FOV' variable is present in the metadata
    if ("fov" %in% names(colData(sce))) {
      # => first, check if it is a case-sensitivity issue
      col_idx <- which(names(colData(sce)) == "fov")
      names(colData(sce))[col_idx] <- "FOV"
    } else if ("instance_id" %in% names(colData(sce))) {
      # => next, check if 'instance_id' exists and follows the AtoMx pattern
      #    and extract FOV from it
      x <- colData(sce)$instance_id[1]
      if( grepl("^c_[0-9]+_[0-9]+_[0-9]+$", x) ) {
        colData(sce)$FOV <- tstrsplit(colData(sce)$instance_id, "_")[[3]]
      } else { 
        stop("cannot run Bruker FOV QC: 'FOV' variable not found in metadata")
        }
    } else {
      stop("cannot run Bruker FOV QC: 'FOV' variable not found in metadata")
    }
    
  }
  
  # Source functions and barcodes from Bruker Spatial Biology:
  # https://github.com/Nanostring-Biostats/CosMx-Analysis-Scratch-Space/tree/Main/_code/FOV%20QC
  if (!file.exists(opt$brukerCode)) {
    stop("Bruker code file not found: ", opt$brukerCode)
  } else { source(opt$brukerCode) }
  if (!file.exists(opt$brukerData)) {
    stop("Bruker barcodes file not found: ", opt$brukerData)
  } else { all_panels <- readRDS(opt$brukerData) }
  
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
  dfov <- as.data.frame(colData(sce)[, c("sample_id", "FOV")])
  dfov$FOV <- as.numeric(as.character(dfov$FOV))
  dfov <- dfov[!duplicated(dfov), ] |> dplyr::arrange(FOV)
  rownames(dfov) <- 1:nrow(dfov)
  
  print(paste0("Checking FOV instrument failures for sample ", pathStem))
  res_fov_qc <- runFOVQC(counts = t(assay(sce)), 
                         xy = data.frame(colData(sce)[, c("x", "y")]),
                         fov = colData(sce)[, "FOV"], 
                         barcodemap = barcodes,  # REMINDER: this only considers probes in the custom panel (including Negative)
                         max_prop_loss = 0.6, max_totalcounts_loss = 0.6)  # default 0.6, the higher, the more relaxed the QC
  #summary(res_fov_qc)

  # Extract and save list of affected genes
  dfg <- res_fov_qc$flagged_fov_x_gene |> data.frame()
  #dfg <- dfg[!duplicated(dfg), ]
  #head(dfg)
  
  if (nrow(dfg) > 0) {
    
    dfg$count <- 1
    dfg <- aggregate.data.frame(dfg, . ~ fov + gene, FUN = "sum")
    names(dfg) <- c("FOV", "gene", "failed_cycles")
    dfg$FOV <- as.numeric(as.character(dfg$FOV))

    dfg <- dfg |>
      dplyr::arrange(gene) |>
      dplyr::arrange(desc(failed_cycles)) |>
      dplyr::arrange(fov)

    dfov <- plyr::join(fov, dfg, by = "FOV", type = "full")

    # Compute fraction of failed genes per FOV
    #dfract <- data.frame(table(dfg$fov) / nrow(barcodes))
    #names(dfract) <- c("FOV", "fraction_genes_bias")
    #dfract$FOV <- as.numeric(as.character(dfract$FOV))
      ## WARNING: unclear how a FOV can fail for gene bias but have no gene flagged as biased?!
  
  }
  
  dfov$qcFlagInstr <- ifelse(dfov$FOV %in% res_fov_qc$flaggedfovs, "Fail", "Pass")
  dfov$qcFlagGeneBias <- ifelse(dfov$FOV %in% res_fov_qc$flaggedfovs_forbias, "Fail", "Pass")
    
  write.csv(dfov, row.names = FALSE, quote = FALSE, 
            file = paste0(opt$outdir, "/", pathStem, "_brukerQCresults.csv"))

  print(paste0("Bruker FOV QC completed for sample ", pathStem))
  
}




# ---------- Task 3: Find top and least expressed genes ----------

v <- rowSums(counts_mat) |> sort(decreasing = TRUE)
v10 <- v[-grep(opt$poscontrol, names(v))] |> head(n = 10)
print("Top most detected probes in this sample (excluding positive controls): ")
v10

u <- rowSums(counts_mat) |> sort(decreasing = FALSE)
u10 <- u[-grep(opt$poscontrol, names(u))] |> head(n = 10)
print("Least detected probes in this sample (excluding positive controls): ")
u10




# ----------

print("Done calculating QC metrics within the following R environment:")
sessionInfo()