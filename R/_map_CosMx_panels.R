#!/usr/bin/Rscript --vanilla

# Script to support the mapping of CosMx probe panels to key Ensembl genome databases
# (versions 93, 98 and 110, to match major Cellranger reference genome releases)

setwd("~/devel/spatialhub/")


############################## Mouse UCC 1k panel ##############################

# Import CosMx probe panel description (negative control removed)
df <- read.delim("data/LBL-11176-05-Mouse-Universal-Cell-Characterization-Gene-List.txt")
names(df) <- c("probe_name", "human_ortholog", "gene_symbol", "description", "gene_alias", "panel_type")
df[df$gene_symbol == "", ]  # check there are no more negative control probes
head(df)

# Convert to long format => one gene_symbol per row
df_long <- df |> 
  tidyr::separate_rows(gene_symbol, sep = ",")


### v87

ahub <- read.delim("data/AnnotHub_databases/mm_EnsDb_v87.tsv")
df_long$gene_symbol[which(!(df_long$gene_symbol %in% ahub$gene))]

df_long$gene_v87 <- df_long$gene_symbol

"1110008P14Rik" %in% ahub$gene
df_long$gene_v87[df_long$gene_symbol == "Bbln"] <- "1110008P14Rik"

ahub$gene[grep("Ccl21", ahub$gene)] |> sort() |> unique()
df_long$gene_v87[df_long$gene_symbol == "Ccl21d"] <- NA

"Cd244" %in% ahub$gene
df_long$gene_v87[df_long$gene_symbol == "Cd244a"] <- "Cd244"

ahub$gene[grep("^H2a", ahub$gene)] |> sort() |> unique()
df_long$gene_v87[df_long$gene_symbol == "H2az1"] <- "H2afz"

ahub$gene[grep("^Hbb", ahub$gene)] |> sort() |> unique()
df_long$gene_v87[df_long$gene_symbol == "Hbb"] <- "Hbb-bh0,Hbb-bh1,Hbb-bh2,Hbb-bh3,Hbb-bs,Hbb-bt"

ahub$gene[grep("^Ifn", ahub$gene)] |> sort() |> unique()
df_long$gene_v87[df_long$gene_symbol == "Ifna"] <- "Ifna1"

"Il1f9" %in% ahub$gene
df_long$gene_v87[df_long$gene_symbol == "Il36g"] <- "Il1f9"

"9530003J23Rik" %in% ahub$gene
df_long$gene_v87[df_long$gene_symbol == "Lyz3"] <- "9530003J23Rik"

ahub$gene[grep("^Plac9", ahub$gene)] |> sort() |> unique()
df_long$gene_v87[df_long$gene_symbol == "Plac9"] <- "Plac9a,Plac9b"

"Ddx58" %in% ahub$gene
df_long$gene_v87[df_long$gene_symbol == "Rigi"] <- "Ddx58"

# Convert to long format => one gene_symbol per row
df_long <- df_long |> 
  tidyr::separate_rows(gene_v87, sep = ",")
df_long


### v93

ahub <- read.delim("data/AnnotHub_databases/mm_EnsDb_v93.tsv")
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
df_long$gene_v93[df_long$gene_symbol == "Hbb" & df_long$gene_v87 == "Hbb-bh0"] <- "Hbb-bh0"
df_long$gene_v93[df_long$gene_symbol == "Hbb" & df_long$gene_v87 == "Hbb-bh1"] <- "Hbb-bh1"
df_long$gene_v93[df_long$gene_symbol == "Hbb" & df_long$gene_v87 == "Hbb-bh2"] <- "Hbb-bh2"
df_long$gene_v93[df_long$gene_symbol == "Hbb" & df_long$gene_v87 == "Hbb-bh3"] <- "Hbb-bh3"
df_long$gene_v93[df_long$gene_symbol == "Hbb" & df_long$gene_v87 == "Hbb-bs"] <- "Hbb-bs"
df_long$gene_v93[df_long$gene_symbol == "Hbb" & df_long$gene_v87 == "Hbb-bt"] <- "Hbb-bt"

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


### v98

ahub <- read.delim("data/AnnotHub_databases/mm_EnsDb_v98.tsv")
df_long$gene_symbol[which(!(df_long$gene_symbol %in% ahub$gene))]

df_long$gene_v98 <- df_long$gene_symbol

"1110008P14Rik" %in% ahub$gene
df_long$gene_v98[df_long$gene_symbol == "Bbln"] <- "1110008P14Rik"

ahub$gene[grep("Ccl21", ahub$gene)] |> sort() |> unique()
df_long$gene_v98[df_long$gene_symbol == "Ccl21d"] <- NA

ahub$gene[grep("^H2a", ahub$gene)] |> sort() |> unique()
df_long$gene_v98[df_long$gene_symbol == "H2az1"] <- "H2afz"

ahub$gene[grep("^Hbb", ahub$gene)] |> sort() |> unique()
df_long$gene_v98[df_long$gene_symbol == "Hbb" & df_long$gene_v93 == "Hbb-bh0"] <- "Hbb-bh0"
df_long$gene_v98[df_long$gene_symbol == "Hbb" & df_long$gene_v93 == "Hbb-bh1"] <- "Hbb-bh1"
df_long$gene_v98[df_long$gene_symbol == "Hbb" & df_long$gene_v93 == "Hbb-bh2"] <- "Hbb-bh2"
df_long$gene_v98[df_long$gene_symbol == "Hbb" & df_long$gene_v93 == "Hbb-bh3"] <- "Hbb-bh3"
df_long$gene_v98[df_long$gene_symbol == "Hbb" & df_long$gene_v93 == "Hbb-bs"] <- "Hbb-bs"
df_long$gene_v98[df_long$gene_symbol == "Hbb" & df_long$gene_v93 == "Hbb-bt"] <- "Hbb-bt"

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

ahub <- read.delim("data/AnnotHub_databases/mm_EnsDb_v110.tsv")
df_long$gene_symbol[which(!(df_long$gene_symbol %in% ahub$gene))]

df_long$gene_v110 <- df_long$gene_symbol

ahub$gene[grep("^Hbb", ahub$gene)] |> sort() |> unique()
df_long$gene_v110[df_long$gene_symbol == "Hbb" & df_long$gene_v98 == "Hbb-bh0"] <- "Hbb-bh0"
df_long$gene_v110[df_long$gene_symbol == "Hbb" & df_long$gene_v98 == "Hbb-bh1"] <- "Hbb-bh1"
df_long$gene_v110[df_long$gene_symbol == "Hbb" & df_long$gene_v98 == "Hbb-bh2"] <- "Hbb-bh2"
df_long$gene_v110[df_long$gene_symbol == "Hbb" & df_long$gene_v98 == "Hbb-bh3"] <- "Hbb-bh3"
df_long$gene_v110[df_long$gene_symbol == "Hbb" & df_long$gene_v98 == "Hbb-bs"] <- "Hbb-bs"
df_long$gene_v110[df_long$gene_symbol == "Hbb" & df_long$gene_v98 == "Hbb-bt"] <- "Hbb-bt"

ahub$gene[grep("^Ifn", ahub$gene)] |> sort() |> unique()
df_long$gene_v110[df_long$gene_symbol == "Ifna"] <- "Ifna1"

"9530003J23Rik" %in% ahub$gene
df_long$gene_v110[df_long$gene_symbol == "Lyz3"] <- "9530003J23Rik"


# Save mapping keys

df <- df_long |> dplyr::select(probe_name, gene_v87, 
                               gene_v93, gene_v98, gene_v110,
                               panel_type)
df <- df[!duplicated(df), ]
write.csv(df, row.names = FALSE, quote = FALSE, 
          file = "data/cosmx_mouse_ucc_1k_ensembl_mapping.csv")




############################## Human UCC 1k panel ##############################

# Import CosMx probe panel description (negative control removed)
df <- read.delim("data/LBL-11178-03-Human-Universal-Cell-Characterization-Panel-Gene-Target-List.txt")
names(df) <- c("probe_name", "gene_symbol", "description", "gene_alias", "panel_type")
df <- df[, c("probe_name", "gene_symbol", "description", "gene_alias", "panel_type")]
df[df$gene_symbol == "", ]  # check there are no more negative control probes
head(df)

# Convert to long format => one gene_symbol per row
df_long <- df |> 
  tidyr::separate_rows(gene_symbol, sep = " | ") |>
  dplyr::filter(gene_symbol != "|")


### v87

ahub <- read.delim("data/AnnotHub_databases/hs_EnsDb_v87.tsv")
df_long$gene_symbol[which(!(df_long$gene_symbol %in% ahub$gene))]

df_long$gene_v87 <- df_long$gene_symbol

ahub$gene[grep("^ATP5B|^ATPSB", ahub$gene)] |> sort() |> unique()
df_long$gene_v87[df_long$gene_symbol == "ATP5F1B"] <- "ATP5B"

ahub$gene[grep("^ATP5E|^ATPSE", ahub$gene)] |> sort() |> unique()
df_long$gene_v87[df_long$gene_symbol == "ATP5F1E"] <- "ATP5E"
df_long$gene_v87[df_long$gene_symbol == "SLMO2-ATP5E"] <- "ATP5E"

"C9orf16" %in% ahub$gene
df_long$gene_v87[df_long$gene_symbol == "BBLN"] <- "C9orf16"

ahub$gene[grep("^SCYA4L|^CCL4L", ahub$gene)] |> sort() |> unique()
df_long$gene_v87[df_long$gene_symbol == "CCL4L1"] <- NA  # probe capturing CCL4L2 as well

ahub$gene[grep("^CD8B", ahub$gene)] |> sort() |> unique()
df_long$gene_v87[df_long$gene_symbol == "CD8B2"] <- "CD8BP"

ahub$gene[grep("CUZD1", ahub$gene)] |> sort() |> unique()
df_long$gene_v87[df_long$gene_symbol == "FAM24B-CUZD1"] <- "CUZD1"

"DDX58" %in% ahub$gene
df_long$gene_v87[df_long$gene_symbol == "RIGI"] <- "DDX58"

ahub$gene[grep("^FCGR3", ahub$gene)] |> sort() |> unique()
df_long$gene_v87[df_long$gene_symbol == "LOC124905743"] <- "FCGR3B"

ahub$gene[grep("^FYB|^SLAP|^ADAP", ahub$gene)] |> sort() |> unique()
df_long$gene_v87[df_long$gene_symbol == "FYB1"] <- "FYB"

ahub$gene[grep("^H2AZ|H2AFZ", ahub$gene)] |> sort() |> unique()
df_long$gene_v87[df_long$gene_symbol == "H2AZ1"] <- "H2AFZ"

ahub$gene[grep("^H4FG|HIST1H4C", ahub$gene)] |> sort() |> unique()
df_long$gene_v87[df_long$gene_symbol == "H4C3"] <- "HIST1H4C"

ahub$gene[grep("^HMGN2$|^TRMT5|C15orf21", ahub$gene)] |> sort() |> unique()
df_long$gene_v87[df_long$gene_symbol == "LOC107986383"] <- "HMGN2"

"ICOSLG" %in% ahub$gene
df_long$gene_v87[df_long$gene_symbol == "LOC102723996"] <- "ICOSLG"

ahub$gene[grep("^LGALS9", ahub$gene)] |> sort() |> unique()
df_long$gene_v87[df_long$gene_symbol == "LOC124907928"] <- "LGALS9C"

ahub$gene[grep("HSALNG000469|NONHSAG001906.2", ahub$gene)] |> sort() |> unique()
df_long$gene_v87[df_long$gene_symbol == "LINC01781"] <- NA

"AC079767.4" %in% ahub$gene
df_long$gene_v87[df_long$gene_symbol == "LINC01857"] <- "AC079767.4"

ahub$gene[grep("HSALNG0089242|NONHSAG010486.2", ahub$gene)] |> sort() |> unique()
df_long$gene_v87[df_long$gene_symbol == "LINC02446"] <- NA

ahub$gene[grep("^NACA", ahub$gene)] |> sort() |> unique()
df_long$gene_v87[df_long$gene_symbol == "NACA4P"] <- "NACAP1"

ahub$gene[grep("^PLAC8", ahub$gene)] |> sort() |> unique()
df_long$gene_v87[df_long$gene_symbol == "LOC107987169"] <- "PLAC8"
df_long$gene_v87[df_long$gene_symbol == "LOC124905457"] <- "PLAC8"

ahub$gene[grep("^PPIA", ahub$gene)] |> sort() |> unique()
df_long$gene_v87[df_long$gene_symbol == "PPIAL4H"] <- "PPIAL4C"
df_long$gene_v87[df_long$gene_symbol == "PPIAP40"] <- "PPIA"
df_long$gene_v87[df_long$gene_symbol == "PPIAP46"] <- "PPIA"
df_long$gene_v87[df_long$gene_symbol == "PPIAP59"] <- "PPIA"
df_long$gene_v87[df_long$gene_symbol == "PPIAP60"] <- "PPIA"

ahub$gene[grep("^SCG5", ahub$gene)] |> sort() |> unique()
df_long$gene_v87[df_long$gene_symbol == "ARHGAP11A-SCG5"] <- "SCG5"

ahub$gene[grep("^VEGF", ahub$gene)] |> sort() |> unique()
df_long$gene_v87[df_long$gene_symbol == "PIR-FIGF"] <- "VEGFD"

ahub$gene[grep("^C10orf54", ahub$gene)] |> sort() |> unique()
df_long$gene_v87[df_long$gene_symbol == "VSIR"] <- "C10orf54"

ahub$gene[grep("^CSDAP1", ahub$gene)] |> sort() |> unique()
df_long$gene_v87[df_long$gene_symbol == "YBX3P1"] <- "CSDAP1"

df_long$gene_v87[which(!(df_long$gene_v87 %in% ahub$gene))]


### v91

ahub <- read.delim("data/AnnotHub_databases/hs_EnsDb_v91.tsv")
df_long$gene_symbol[which(!(df_long$gene_symbol %in% ahub$gene))]

df_long$gene_v91 <- df_long$gene_symbol

ahub$gene[grep("^ATP5B|^ATPSB", ahub$gene)] |> sort() |> unique()
df_long$gene_v91[df_long$gene_symbol == "ATP5F1B"] <- "ATP5B"

ahub$gene[grep("^ATP5E|^ATPSE", ahub$gene)] |> sort() |> unique()
df_long$gene_v91[df_long$gene_symbol == "ATP5F1E"] <- "ATP5E"
df_long$gene_v91[df_long$gene_symbol == "SLMO2-ATP5E"] <- "ATP5E"

"C9orf16" %in% ahub$gene
df_long$gene_v91[df_long$gene_symbol == "BBLN"] <- "C9orf16"

ahub$gene[grep("^SCYA4L|^CCL4L", ahub$gene)] |> sort() |> unique()
df_long$gene_v91[df_long$gene_symbol == "CCL4L1"] <- NA  # probe capturing CCL4L2 as well

ahub$gene[grep("CUZD1", ahub$gene)] |> sort() |> unique()
df_long$gene_v91[df_long$gene_symbol == "FAM24B-CUZD1"] <- "CUZD1"

"DDX58" %in% ahub$gene
df_long$gene_v91[df_long$gene_symbol == "RIGI"] <- "DDX58"

ahub$gene[grep("^FCGR3", ahub$gene)] |> sort() |> unique()
df_long$gene_v91[df_long$gene_symbol == "LOC124905743"] <- "FCGR3B"

ahub$gene[grep("^H2AZ|H2AFZ", ahub$gene)] |> sort() |> unique()
df_long$gene_v91[df_long$gene_symbol == "H2AZ1"] <- "H2AFZ"

ahub$gene[grep("^H4FG|HIST1H4C", ahub$gene)] |> sort() |> unique()
df_long$gene_v91[df_long$gene_symbol == "H4C3"] <- "HIST1H4C"

ahub$gene[grep("^HMGN2$|^TRMT5|C15orf21", ahub$gene)] |> sort() |> unique()
df_long$gene_v91[df_long$gene_symbol == "LOC107986383"] <- "HMGN2"

"ICOSLG" %in% ahub$gene
df_long$gene_v91[df_long$gene_symbol == "LOC102723996"] <- "ICOSLG"

ahub$gene[grep("^LGALS9", ahub$gene)] |> sort() |> unique()
df_long$gene_v91[df_long$gene_symbol == "LOC124907928"] <- "LGALS9C"

ahub$gene[grep("^NACA", ahub$gene)] |> sort() |> unique()
df_long$gene_v91[df_long$gene_symbol == "NACA4P"] <- "NACAP1"

ahub$gene[grep("^PLAC8", ahub$gene)] |> sort() |> unique()
df_long$gene_v91[df_long$gene_symbol == "LOC107987169"] <- "PLAC8"
df_long$gene_v91[df_long$gene_symbol == "LOC124905457"] <- "PLAC8"

ahub$gene[grep("^PPIA", ahub$gene)] |> sort() |> unique()
df_long$gene_v91[df_long$gene_symbol == "PPIAL4H"] <- "PPIAL4C"

ahub$gene[grep("^SCG5", ahub$gene)] |> sort() |> unique()
df_long$gene_v91[df_long$gene_symbol == "ARHGAP11A-SCG5"] <- "SCG5"

ahub$gene[grep("^VEGF", ahub$gene)] |> sort() |> unique()
df_long$gene_v91[df_long$gene_symbol == "PIR-FIGF"] <- "VEGFD"

ahub$gene[grep("^YBX3", ahub$gene)] |> sort() |> unique()
df_long$gene_v91[df_long$gene_symbol == "YBX3P1"] <- "YBX3"

df_long$gene_v91[which(!(df_long$gene_v91 %in% ahub$gene))]


### v93

ahub <- read.delim("data/AnnotHub_databases/hs_EnsDb_v93.tsv")
df_long$gene_symbol[which(!(df_long$gene_symbol %in% ahub$gene))]

df_long$gene_v93 <- df_long$gene_symbol

ahub$gene[grep("^ATP5F1B|^ATP5F1E|^ATP5E|^ATPSE", ahub$gene)] |> sort() |> unique()
df_long$gene_v93[df_long$gene_symbol == "SLMO2-ATP5E"] <- "ATP5F1E"

"C9orf16" %in% ahub$gene
df_long$gene_v93[df_long$gene_symbol == "BBLN"] <- "C9orf16"

ahub$gene[grep("^SCYA4L|^CCL4L", ahub$gene)] |> sort() |> unique()
df_long$gene_v93[df_long$gene_symbol == "CCL4L1"] <- NA  # probe capturing CCL4L2 as well

ahub$gene[grep("CUZD1", ahub$gene)] |> sort() |> unique()
df_long$gene_v93[df_long$gene_symbol == "FAM24B-CUZD1"] <- "CUZD1"

"DDX58" %in% ahub$gene
df_long$gene_v93[df_long$gene_symbol == "RIGI"] <- "DDX58"

ahub$gene[grep("^FCGR3", ahub$gene)] |> sort() |> unique()
df_long$gene_v93[df_long$gene_symbol == "LOC124905743"] <- "FCGR3B"

ahub$gene[grep("^H2AZ|H2AFZ", ahub$gene)] |> sort() |> unique()
df_long$gene_v93[df_long$gene_symbol == "H2AZ1"] <- "H2AFZ"

ahub$gene[grep("^H4FG|HIST1H4C", ahub$gene)] |> sort() |> unique()
df_long$gene_v93[df_long$gene_symbol == "H4C3"] <- "HIST1H4C"

ahub$gene[grep("^HMGN2$|^TRMT5|C15orf21", ahub$gene)] |> sort() |> unique()
df_long$gene_v93[df_long$gene_symbol == "LOC107986383"] <- "HMGN2"

"ICOSLG" %in% ahub$gene
df_long$gene_v93[df_long$gene_symbol == "LOC102723996"] <- "ICOSLG"

ahub$gene[grep("^LGALS9", ahub$gene)] |> sort() |> unique()
df_long$gene_v93[df_long$gene_symbol == "LOC124907928"] <- "LGALS9C"

ahub$gene[grep("^NACA", ahub$gene)] |> sort() |> unique()
df_long$gene_v93[df_long$gene_symbol == "NACA4P"] <- "NACAP1"

ahub$gene[grep("^PLAC8", ahub$gene)] |> sort() |> unique()
df_long$gene_v93[df_long$gene_symbol == "LOC107987169"] <- "PLAC8"
df_long$gene_v93[df_long$gene_symbol == "LOC124905457"] <- "PLAC8"

ahub$gene[grep("^PPIA", ahub$gene)] |> sort() |> unique()
df_long$gene_v93[df_long$gene_symbol == "PPIAL4H"] <- "PPIAL4C"

ahub$gene[grep("^SCG5", ahub$gene)] |> sort() |> unique()
df_long$gene_v93[df_long$gene_symbol == "ARHGAP11A-SCG5"] <- "SCG5"

ahub$gene[grep("^VEGF", ahub$gene)] |> sort() |> unique()
df_long$gene_v93[df_long$gene_symbol == "PIR-FIGF"] <- "VEGFD"

df_long$gene_v93[which(!(df_long$gene_v93 %in% ahub$gene))]


### v98

ahub <- read.delim("data/AnnotHub_databases/hs_EnsDb_v98.tsv")
df_long$gene_symbol[which(!(df_long$gene_symbol %in% ahub$gene))]

df_long$gene_v98 <- df_long$gene_symbol

ahub$gene[grep("^ATP5F1B|^ATP5F1E|^ATP5E|^ATPSE", ahub$gene)] |> sort() |> unique()
df_long$gene_v98[df_long$gene_symbol == "SLMO2-ATP5E"] <- "ATP5F1E"

"C9orf16" %in% ahub$gene
df_long$gene_v98[df_long$gene_symbol == "BBLN"] <- "C9orf16"

ahub$gene[grep("^SCYA4L|^CCL4L", ahub$gene)] |> sort() |> unique()
df_long$gene_v98[df_long$gene_symbol == "CCL4L1"] <- NA  # probe capturing CCL4L2 as well

ahub$gene[grep("CUZD1", ahub$gene)] |> sort() |> unique()
df_long$gene_v98[df_long$gene_symbol == "FAM24B-CUZD1"] <- "CUZD1"

"DDX58" %in% ahub$gene
df_long$gene_v98[df_long$gene_symbol == "RIGI"] <- "DDX58"

ahub$gene[grep("^FCGR3", ahub$gene)] |> sort() |> unique()
df_long$gene_v98[df_long$gene_symbol == "LOC124905743"] <- "FCGR3B"

ahub$gene[grep("^H2AZ|H2AFZ", ahub$gene)] |> sort() |> unique()
df_long$gene_v98[df_long$gene_symbol == "H2AZ1"] <- "H2AFZ"

ahub$gene[grep("^H4FG|HIST1H4C", ahub$gene)] |> sort() |> unique()
df_long$gene_v98[df_long$gene_symbol == "H4C3"] <- "HIST1H4C"

ahub$gene[grep("^HMGN2$|^TRMT5|C15orf21", ahub$gene)] |> sort() |> unique()
df_long$gene_v98[df_long$gene_symbol == "LOC107986383"] <- "HMGN2"

"ICOSLG" %in% ahub$gene
df_long$gene_v98[df_long$gene_symbol == "LOC102723996"] <- "ICOSLG"

ahub$gene[grep("^LGALS9", ahub$gene)] |> sort() |> unique()
df_long$gene_v98[df_long$gene_symbol == "LOC124907928"] <- "LGALS9C"

ahub$gene[grep("^PLAC8", ahub$gene)] |> sort() |> unique()
df_long$gene_v98[df_long$gene_symbol == "LOC107987169"] <- "PLAC8"
df_long$gene_v98[df_long$gene_symbol == "LOC124905457"] <- "PLAC8"

ahub$gene[grep("^SCG5", ahub$gene)] |> sort() |> unique()
df_long$gene_v98[df_long$gene_symbol == "ARHGAP11A-SCG5"] <- "SCG5"

ahub$gene[grep("^VEGF", ahub$gene)] |> sort() |> unique()
df_long$gene_v98[df_long$gene_symbol == "PIR-FIGF"] <- "VEGFD"

df_long$gene_v98[which(!(df_long$gene_v98 %in% ahub$gene))]


### v110

ahub <- read.delim("data/AnnotHub_databases/hs_EnsDb_v110.tsv")
df_long$gene_symbol[which(!(df_long$gene_symbol %in% ahub$gene))]

df_long$gene_v110 <- df_long$gene_symbol

ahub$gene[grep("^ATP5F1B|^ATP5F1E|^ATP5E|^ATPSE", ahub$gene)] |> sort() |> unique()
df_long$gene_v110[df_long$gene_symbol == "SLMO2-ATP5E"] <- "ATP5F1E"

ahub$gene[grep("^SCYA4L|^CCL4L", ahub$gene)] |> sort() |> unique()
df_long$gene_v110[df_long$gene_symbol == "CCL4L1"] <- NA  # probe capturing CCL4L2 as well

ahub$gene[grep("CUZD1", ahub$gene)] |> sort() |> unique()
df_long$gene_v110[df_long$gene_symbol == "FAM24B-CUZD1"] <- "CUZD1"

ahub$gene[grep("^FCGR3", ahub$gene)] |> sort() |> unique()
df_long$gene_v110[df_long$gene_symbol == "LOC124905743"] <- "FCGR3B"

ahub$gene[grep("^HMGN2$|^TRMT5|C15orf21", ahub$gene)] |> sort() |> unique()
df_long$gene_v110[df_long$gene_symbol == "LOC107986383"] <- "HMGN2"

"ICOSLG" %in% ahub$gene
df_long$gene_v110[df_long$gene_symbol == "LOC102723996"] <- "ICOSLG"

ahub$gene[grep("^LGALS9", ahub$gene)] |> sort() |> unique()
df_long$gene_v110[df_long$gene_symbol == "LOC124907928"] <- "LGALS9C"

ahub$gene[grep("^PLAC8", ahub$gene)] |> sort() |> unique()
df_long$gene_v110[df_long$gene_symbol == "LOC107987169"] <- "PLAC8"
df_long$gene_v110[df_long$gene_symbol == "LOC124905457"] <- "PLAC8"

ahub$gene[grep("^VEGF", ahub$gene)] |> sort() |> unique()
df_long$gene_v110[df_long$gene_symbol == "PIR-FIGF"] <- "VEGFD"

df_long$gene_v110[which(!(df_long$gene_v110 %in% ahub$gene))]


# Save mapping keys

df <- df_long |> dplyr::select(probe_name, gene_v87, gene_v91, 
                               gene_v93, gene_v98, gene_v110,
                               panel_type)
df <- df[!duplicated(df), ]
write.csv(df, row.names = FALSE, quote = FALSE, 
          file = "data/cosmx_human_ucc_1k_ensembl_mapping.csv")
