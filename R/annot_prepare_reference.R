#!/usr/bin/Rscript --vanilla

## < Add Description and usage info >

message("Pre-processing reference dataset for cell type annotation")
timestamp()


# ---------- Libraries ----------

stopifnot(
  require(optparse),
  require(tidyverse),
  require(data.table),
  require(zellkonverter),
  require(SingleCellExperiment),
  require(Matrix)
)


# ---------- Options ----------

option_list <- list(
    make_option(
        c("--atlasKey"),
        help="key of reference dataset to pre-process"
    ),
    make_option(
        c("--atlasTSV"),
        help="path to the spatialhub TSV table describing reference dataset(s) to use for cell annotation"
    ),
    make_option(
        c("--reduceMeta"),
        help="whether to reduce the metadata to its minimal variables necessary for annotation"
    ),
    make_option(
        c("--path2mapping"),
        help="path to probe mapping metadata"
    ),
    make_option(
        c("--annothubDir"),
        help="path to directory of Ensembl annotation datasets"
    )
)

opt <- parse_args(OptionParser(option_list=option_list))
set.seed(123456789)

cat("Running with options:\n")
print(opt)

if (!dir.exists("annot.dir/atlas.dir")){
  dir.create("annot.dir/atlas.dir")
}


# ---------- Setup ----------

# Read in (spatialhub) metadata file
path2meta <- opt$atlasTSV
df <- read.delim(path2meta, sep = "\t")

# Subet to reference dataset of interest
df <- df |> dplyr::filter(atlas_id == opt$atlasKey)
print(df)

# Define metadata variables to keep accordingly
keepVar <- c(df$lineage_key, df$celltype_annot_key, df$celltype_other_key,
             df$sample_key, df$donor_key,
             df$hvg_batch, df$scvi_batch,
             df$categorical_covar, df$continuous_covar)
keepVar <- keepVar[keepVar != 'none']


### Read in `SingleCellExperiment` object to be processed

stopifnot("file not found" = file.exists(df$path))
print(paste0("Working on dataset ", df$atlas_id, " from ", df$source))

v <- unlist(strsplit(basename(df$path), "\\."))
file_extension <- tolower(v[length(v)])

if (file_extension == "h5ad") {

  sce <- zellkonverter::readH5AD(df$path, verbose = TRUE)

} else if (file_extension == "rds") {
  
  sce <- readRDS(df$path)
  
  # if Seurat object, convert to SingleCellExperiment
  if (length(grep("Seurat", class(sce))) > 0) {
    hvg <- Seurat::VariableFeatures(sce)
    sce <- SingleCellExperiment(assays = list(counts = sce@assays$RNA@counts),
                                colData = sce@meta.data)
    metadata(sce)$hvg <- hvg
  }

} else {
  stop("Please provide a `.h5ad` or `.rds` file as input.")
}

print("Working with the following `SingleCellExperiment` object:")
print(sce)




# ---------- Tasks ----------

### Clean up metadata (optional)

stopifnot("Variables specified in atlas table not all found in metadata of provided atlas object" = all(keepVar %in% names(colData(sce))))

coldf <- as.data.frame(colData(sce))

if (opt$reduceMeta) {
    coldf <- coldf |> dplyr::select(all_of(keepVar))
}
#coldf$barcode_id <- rownames(coldf)

if ( all(colnames(sce) == rownames(coldf)) ) {
  colData(sce) <- DataFrame(coldf)
} else {
    stop("Column mismatch!")
}


### Clean up gene annotations, mostly for matching with probe panel names

if (df$species == "mouse") {
    species <- "mm"
} else if (df$species == "human") {
    species <- "hs"
} else {
    stop("Species not supported!")
}
annotFile <- paste0(opt$annothubDir, '/', species, 
                    "_EnsDb_v", df$ensembl_version, ".tsv")

print(paste0("Using the following Ensembl database for gene mapping: ", annotFile))
ahub <- read.delim(annotFile, sep = "\t")
stopifnot("Key not found in provided Ensembl annotation dataset (must be one of gene, ensembl_id, entrez_id)." = df$gene_key %in% names(ahub))


# Using *rownames* to map reference dataset genes to Ensembl database
rowdf <- data.frame("original_name" = rownames(sce), "updated_name" = rownames(sce))
sce_genes <- rownames(sce)

idx_mismatch <- which(!(sce_genes %in% ahub[, df$gene_key]))
n_mismatch <- length(idx_mismatch)
print(paste0("Number of genes IDs not found in selected Ensembl database: ", n_mismatch))

if (df$gene_key == "gene") {

  # Handle the case of suspected 1:many mappings < gene_name.1 >
  # where multiple transcripts with distinct ensembl_ids map to the same gene
  for (gene in sce_genes[idx_mismatch]) {
    bgene <- tstrsplit(gene, "\\.")[[1]]  # retrieve gene "basename" 
    if (bgene %in% sce_genes) {
      rowdf$updated_name[rowdf$original_name == gene] <- bgene
    }
  }

}

# Now, let's match our rownames to the ensembl database
rowdf[, df$gene_key] <- rowdf$updated_name
rowdf <- plyr::join(rowdf, ahub, by = df$gene_key)

# Remove potential duplicates due to 1:many mappings in gene:ensembl_id, keeping only the first ensembl_id
rowdf <- rowdf |> dplyr::arrange(ensembl_id) |>
  dplyr::select(original_name, updated_name, gene, 
                ensembl_id, entrez_id, gene_biotype)
rowdf <- rowdf[!duplicated(rowdf$original_name), ]

# Set updated name to current *gene* (symbol) in provided Ensembl database
rowdf$updated_name <- rowdf$gene

n_mismatch <- sum(is.na(rowdf$updated_name))
print(paste0("Number of genes IDs that remain unmapped: ", n_mismatch))
print(rowdf$original_name[is.na(rowdf$updated_name)])

# Keep original name when no match was found
rowdf$updated_name[is.na(rowdf$updated_name)] <- rowdf$original_name[is.na(rowdf$updated_name)]
rowdf$gene[is.na(rowdf$gene)] <- rowdf$original_name[is.na(rowdf$gene)]


# Retrieve HVG information if known
if ("highly_variable" %in% names(rowData)) {
  rowdf$highly_variable <- rowData$highly_variable
  HVG_found <- TRUE
} else if ("hvg" %in% names(sce@metadata)) {
  rowdf$highly_variable <- 'False'
  rowdf$highly_variable[rowdf$original_name %in% sce@metadata$hvg] <- 'True'
  HVG_found <- TRUE
} else {
  print("No HVG found for this reference dataset.")
  HVG_found <- FALSE
}


### Retrieve CosMx probe mapping to specified Ensembl database

cosmx_df <- read.csv(opt$path2mapping)
cosmx_df <- cosmx_df[!duplicated(cosmx_df), ]
n_probes <- length(unique(cosmx_df$probe_name))

mapping_key <- paste0("gene_v", df$ensembl_version)
stopifnot("CosMx probe name to gene symbol mapping not supported for provided Ensembl database version." = mapping_key %in% names(cosmx_df))

cosmx_df$gene <- cosmx_df[, mapping_key]
cosmx_df <- cosmx_df |> dplyr::select(probe_name, gene)
rowdf <- plyr::join(rowdf, cosmx_df, by = "gene")

n_covered <- rowdf$probe_name[!is.na(rowdf$probe_name)] |> unique() |> length()
print(paste0("Genes in this reference dataset cover ", n_covered, 
             " distinct probes (out of ", n_probes, " in this panel)."))

# Define grouping variable, i.e. updated name with CosMx match
rowdf$group <- rowdf$updated_name
rowdf$group[!is.na(rowdf$probe_name)] <- rowdf$probe_name[!is.na(rowdf$probe_name)]

# Remove potential duplicates due to 1:many mappings in gene:ensembl_id, keeping only the first ensembl_id
rowdf <- rowdf |> dplyr::arrange(ensembl_id) |>
  dplyr::select(original_name, updated_name, gene, 
                ensembl_id, entrez_id, gene_biotype,
                probe_name, group)
rowdf <- rowdf[!duplicated(rowdf$original_name), ]


### Aggregate counts for gene_ids corresponding to the same probe name

print(paste0("Reducing counts matrix from ", length(unique(rowdf$original_name)), 
             " to ", length(unique(rowdf$group)), " distinct gene_ids"))

if ("counts" %in% names(sce@assays@data)) {
  counts <- sce@assays@data$counts
} else if ("X" %in% names(sce@assays@data)) {
  print("No counts layer found in provided reference atlas. Using 'X' slot instead")
  counts <- sce@assays@data$X
} else {
  stop("No counts layer found in provided reference atlas.")
}
v <- rowSums(counts)
stopifnot("Reference dataset contains normalized values. Please provide raw counts." = all(v - floor(v) == 0))

rownames(counts) <- rowdf$group
colnames(counts) <- colnames(sce)
print(counts[1:10, 1:10])

# Split original counts matrix between genes that need to be aggregated or not
dup_genes <- rowdf$group[duplicated(rowdf$group)]
counts0 <- counts[!(rownames(counts) %in% dup_genes), ]
counts_aggr <- counts[rownames(counts) %in% dup_genes, ]

counts_aggr <- t(sapply(by(as.matrix(counts_aggr), rownames(counts_aggr), colSums), identity))
counts_aggr <- as(counts_aggr, "sparseMatrix")  # convert back to sparse matrix

if (all(colnames(counts0) == colnames(counts_aggr))) {
  counts <- rbind(counts0, counts_aggr)
} else {
  stop("Mismatch in column (cell_id) names")
}
dim(counts)

# Update 'rowdf'
rowdf_aggr <- data.frame(gene = rownames(counts))
if (HVG_found) {
  rowdf_aggr$highly_variable <- 'False'
  
  hvg_group <- rowdf$group[rowdf$highly_variable == 'True']
  rowdf_aggr$highly_variable[rowdf_aggr$gene %in% hvg_group] <- 'True'
}
rowdf <- rowdf_aggr
rownames(rowdf) <- rowdf$gene


### Save updated H5AD

# Store it all back into one updated sce object
gene <- rowdf[match(rownames(counts), rownames(rowdf)), ]

if ( all(rownames(counts) == rownames(rowdf)) & all(colnames(counts) == rownames(colData(sce))) ) {
  sce <- SingleCellExperiment(assays = list(counts = counts),
                              colData = colData(sce),
                              rowData = DataFrame(gene))
} else {
    stop("Mismatch in row or column names!")
}

print("Saving processed reference dataset as .h5ad file for scANVI")
zellkonverter::writeH5AD(sce, compression = "gzip", verbose = TRUE,
                         file = paste0("annot.dir/atlas.dir/", df$atlas_id, ".h5ad"))


# Perform HVG selection, if applicable (?) - results vary greatly from scanpy approach...
#sce <- scuttle::logNormCounts(sce)
#var_sce <- scran::modelGeneVar(sce, block = colData(sce)$mouse)
#hvg_sce <- scran::getTopHVGs(var_sce, n = 2000)
