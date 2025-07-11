#!/usr/bin/Rscript --vanilla

## < Add Description and usage info >

message("Performing cell type annotation with SingleR")


# ---------- Libraries ----------

stopifnot(
  require(optparse),
  require(tidyverse),
  require(data.table),
  require(zellkonverter),
  require(SingleCellExperiment),
  require(SingleR),
  require(scuttle)
)

# ---------- Options ----------

option_list <- list(
    make_option(
        c("--atlasKey"),
        help="key of reference atlas to use for annotation"
    ),
    make_option(
        c("--atlasTSV"),
        help="path to the spatialhub atlas TSV table describing reference datasets to use for annotation"
    ),
    make_option(
        c("--numWorkers"), default = 12,
        help="number of parallel processes to use"
    )
)

opt <- parse_args(OptionParser(option_list=option_list))
set.seed(123456789)

cat("Running with options:\n")
print(opt)


# ---------- Setup ----------

multicoreParam <- MulticoreParam(workers = opt$numWorkers)

# Read in (spatialhub) metadata file
path2meta <- opt$atlasTSV
df <- read.delim(path2meta, sep = "\t")

if (!("sample_name" %in% names(df))) {
    print("'sample_name' not found: setting it to 'sample_id'")
    df$sample_name <- df$sample_id
}

# Set paths to files
df_query <- df |> dplyr::filter(type == 'query')
stopifnot("spatialhub annot pipeline can only handle one query dataset at a time" = nrow(df_query) == 1)
path2spx <- df$path[df$type == 'query']

df <- df |> dplyr::filter(atlas_id == opt$atlasKey)
stopifnot("atlas_key is not unique" = nrow(df) == 1)
path2ref <- df$path

# Define label to use for annotations
annot_label <- df$celltype_annot_key

# Create SingleR out directory if it does not exist
outDir <- paste0("annot.dir/singleR/", opt$atlasKey)
if (!dir.exists(outDir)){
    dir.create(outDir)
}


# ---------- Tasks ----------

### Read in query (spatial) dataset

print(paste0("Importing (spatial) query dataset from: ", path2spx))
sce <- zellkonverter::readH5AD(path2spx, verbose = TRUE)

if (!('counts' %in% names(assays(sce)))) {
    print("No 'counts' assay found. Assuming raw counts are in 'X' slot.")
    idx <- which(names(assays(sce)) == 'X')
    names(assays(sce))[idx] <- 'counts'
}
names(sce)

sce <- scuttle::logNormCounts(sce)
    ## in theory, query dataset does not need to be log-transformed, as only the ranks will be used by SingleR()
    ## in practice, raw counts triggered an error message...


### Read in reference dataset

print(paste0("Importing reference atlas from: ", path2ref))
ref <- zellkonverter::readH5AD(path2ref, verbose = TRUE)
ref <- scuttle::logNormCounts(ref)
    ## reference dataset must be normalized and log-transformed


### Predict cell types in query

matching_genes <- intersect(rownames(sce), rownames(ref))
print(paste0("Found the following number of matching gene_ids between query and reference: ", length(matching_genes)))

print(paste0("Predicting cell type using the following annotation key: ", annot_label))
preds <- SingleR(test = sce, ref = ref, 
                 labels = annot_label, 
                 de.method = 'wilcox',
                 BPPARAM = mutlicoreParam)


### Save results

print("Saving predicted labels")

labels <- data.frame(cell_index = rownames(preds),
                     query = rep(basename(path2spx), nrow(preds)),
                     labels = preds$labels,
                     pruned_labels = preds$pruned.labels)

write.table(labels,
            gzfile(file.path(outDir, paste(basename(path2spx), "labels.tsv.gz", sep="_"))),
            quote=FALSE, col.names=TRUE, row.names=FALSE,
            sep="\t")


print("Saving prediction confidence scores")

scores <- data.frame(cell_index = rownames(preds),
                     query = rep(basename(path2spx), nrow(preds)),
                     delta_next = preds$delta.next)
scores <- cbind(scores, data.frame(preds$scores))

write.table(scores,
            gzfile(file.path(outDir, paste(basename(path2spx), "scores.tsv.gz", sep="_"))),
            quote=FALSE, col.names=TRUE, row.names=FALSE,
            sep="\t")
