#!/usr/bin/Rscript --vanilla

library(AnnotationHub)
library(tidyverse)


hub = AnnotationHub()
#yes  # if directory does not already exist

# Mouse - v93, 98 and 110 (i.e. matching major Cellranger reference releases)
query(hub, c("EnsDb", "Mus musculus", "110"))
edb = hub[["AH113713"]]
edb
AnnotationDbi::keytypes(edb)

# Extract desired information
keys <- keys(edb, "GENEID")  # GENEID = Ensembl gene ID
columns <- c("SYMBOL", "ENTREZID", "GENEBIOTYPE", "DESCRIPTION")
ahub <-
  ensembldb::select(edb, keys, columns, keytype = "GENEID") |>
  as_tibble()
ahub <- ahub[!duplicated(ahub), ] |> dplyr::arrange(GENEID)
names(ahub) <- c("ensembl_id", "gene", "entrez_id", "gene_biotype", "gene_description")

write.table(ahub, row.names = FALSE, quote = FALSE, sep = "\t",
            file = 'data/AnnotHub_databases/mm_EnsDb_v110.csv')


# Human - v93, 98 and 110 (i.e. matching major Cellranger reference releases)
query(hub, c("EnsDb", "Homo sapiens", "110"))
edb = hub[["AH113665"]]
edb
AnnotationDbi::keytypes(edb)

# Extract desired information
keys <- keys(edb, "GENEID")  # GENEID = Ensembl gene ID
columns <- c("SYMBOL", "ENTREZID", "GENEBIOTYPE", "DESCRIPTION")
ahub <-
  ensembldb::select(edb, keys, columns, keytype = "GENEID") |>
  as_tibble()
ahub <- ahub[!duplicated(ahub), ] |> dplyr::arrange(GENEID)
names(ahub) <- c("ensembl_id", "gene", "entrez_id", "gene_biotype", "gene_description")

#data/AnnotHub_databases
write.table(ahub, row.names = FALSE, quote = FALSE, sep = "\t",
            file = 'reference-data/hs_EnsDb_v110.csv')