#!/usr/bin/Rscript --vanilla

# Script to support the mapping of CosMx probe panels to key Ensembl genome databases
# (versions 93, 98 and 110, to match major Cellranger reference genome releases)


############################## Mouse UCC 1k panel ##############################

# Import CosMx probe panel description (negative control removed)
df <- read.delim("data/LBL-11176-05-Mouse-Universal-Cell-Characterization-Gene-List.txt")
names(df) <- c("probe_name", "human_ortholog", "gene_symbol", "description", "gene_alias", "panel_type")
df[df$gene_symbol == "", ]  # check there are no more negative control probes
head(df)

# Convert to long format => one gene_symbol per row
df_long <- df |> 
  tidyr::separate_rows(gene_symbol, sep = ",")
df


### v93

ahub <- read.delim("data/AnnotHub_databases/mm_EnsDb_v93.csv")
df_long$gene_symbol[which(!(df_long$gene_symbol %in% ahub$gene))]

df_long$gene_v93 <- df_long$gene_symbol

"1110008P14Rik" %in% ahub$gene
df_long$gene_v93[df_long$gene_symbol == "Bbln"] <- "1110008P14Rik"

ahub$gene[grep("Ccl21", ahub$gene)] |> sort() |> unique()
df_long$gene_v93[df_long$gene_symbol == "Ccl21d"] <- NA

"Cd244" %in% ahub$gene
df_long$gene_v93[df_long$gene_symbol == "Cd244a"] <- "Cd244"

ahub$gene[grep("^H2a", ahub$gene)] |> sort() |> unique()
df_long$gene_v93[df_long$gene_symbol == "H2az1"] <- "H2afz"

ahub$gene[grep("^Hbb", ahub$gene)] |> sort() |> unique()
df_long$gene_v93[df_long$gene_symbol == "Hbb"] <- NA

ahub$gene[grep("^Ifn", ahub$gene)] |> sort() |> unique()
df_long$gene_v93[df_long$gene_symbol == "Ifna"] <- "Ifna1"

"Il1f9" %in% ahub$gene
df_long$gene_v93[df_long$gene_symbol == "Il36g"] <- "Il1f9"

"9530003J23Rik" %in% ahub$gene
df_long$gene_v93[df_long$gene_symbol == "Lyz3"] <- "9530003J23Rik"

ahub$gene[grep("^Plac9", ahub$gene)] |> sort() |> unique()
df_long$gene_v93[df_long$gene_symbol == "Plac9"] <- "Plac9a,Plac9b"

"Ddx58" %in% ahub$gene
df_long$gene_v93[df_long$gene_symbol == "Rigi"] <- "Ddx58"

# Convert to long format => one gene_symbol per row
df_long <- df_long |> 
  tidyr::separate_rows(gene_v93, sep = ",")
df_long


### v98

ahub <- read.delim("data/AnnotHub_databases/mm_EnsDb_v98.csv")
df_long$gene_symbol[which(!(df_long$gene_symbol %in% ahub$gene))]

df_long$gene_v98 <- df_long$gene_symbol

"1110008P14Rik" %in% ahub$gene
df_long$gene_v98[df_long$gene_symbol == "Bbln"] <- "1110008P14Rik"

ahub$gene[grep("Ccl21", ahub$gene)] |> sort() |> unique()
df_long$gene_v98[df_long$gene_symbol == "Ccl21d"] <- NA

ahub$gene[grep("^H2a", ahub$gene)] |> sort() |> unique()
df_long$gene_v98[df_long$gene_symbol == "H2az1"] <- "H2afz"

ahub$gene[grep("^Hbb", ahub$gene)] |> sort() |> unique()
df_long$gene_v98[df_long$gene_symbol == "Hbb"] <- NA

ahub$gene[grep("^Ifn", ahub$gene)] |> sort() |> unique()
df_long$gene_v98[df_long$gene_symbol == "Ifna"] <- "Ifna1"

"Il1f9" %in% ahub$gene
df_long$gene_v98[df_long$gene_symbol == "Il36g"] <- "Il1f9"

"9530003J23Rik" %in% ahub$gene
df_long$gene_v98[df_long$gene_symbol == "Lyz3"] <- "9530003J23Rik"

ahub$gene[grep("^Plac9", ahub$gene)] |> sort() |> unique()
df_long$gene_v98[df_long$gene_symbol == "Plac9" & df_long$gene_v93 == "Plac9a"] <- "Plac9a"
df_long$gene_v98[df_long$gene_symbol == "Plac9" & df_long$gene_v93 == "Plac9b"] <- "Plac9b"

"Ddx58" %in% ahub$gene
df_long$gene_v98[df_long$gene_symbol == "Rigi"] <- "Ddx58"


### v110

ahub <- read.delim("data/AnnotHub_databases/mm_EnsDb_v110.csv")
df_long$gene_symbol[which(!(df_long$gene_symbol %in% ahub$gene))]

df_long$gene_v110 <- df_long$gene_symbol

ahub$gene[grep("^Hbb", ahub$gene)] |> sort() |> unique()
df_long$gene_v110[df_long$gene_symbol == "Hbb"] <- NA

ahub$gene[grep("^Ifn", ahub$gene)] |> sort() |> unique()
df_long$gene_v110[df_long$gene_symbol == "Ifna"] <- "Ifna1"

"9530003J23Rik" %in% ahub$gene
df_long$gene_v110[df_long$gene_symbol == "Lyz3"] <- "9530003J23Rik"


# Save mapping keys

df <- df_long |> dplyr::select(probe_name, gene_v93, gene_v98, gene_v110)
write.csv(df, row.names = FALSE, quote = FALSE, 
          file = "data/cosmx_mouse_ucc_1k_ensembl_mapping.csv")
