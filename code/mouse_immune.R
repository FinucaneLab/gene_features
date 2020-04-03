#' ---
#' title: "Create gene features from mouse_immune"
#' author: "Jacob Ulirsch"
#' date: "July 12, 2018"
#' ---

#' 
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = FALSE, eval = TRUE
.libPaths(c(.libPaths(), "/PHShome/cl322/R/x86_64-pc-linux-gnu-library/3.4"))
library(tidyverse)
library(data.table)
#library(MUDAN)
library(BuenColors)
library(Seurat)
library(irlba)
source("code/utils.R")
#library(reticulate)
#use_virtualenv("/data/aryee/julirsch/python/venv3/bin/activate", required = TRUE)

#' Read in data and annotations
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
file <- "data/mouse_immune/17aug_commonRNA_normalized.txt.gz"
mouse_immune <- data.frame(fread(paste0("zcat < ", file), skip = 0, sep = "\t")) %>% 
  data.matrix() %>% 
  Matrix()
row.names(mouse_immune) <- read.table("data/mouse_immune/rnaRowLabels.txt") %>% .$V1
#mouse_immune.annot <- as.data.frame(mouse_immune.annot)
#row.names(mouse_immune.annot) <- mouse_immune.annot$tissue_subcluster

#' Identify human homologs and subset to concise gene set
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
# Convert mouse symbols to human symbols
s2m <- read.table("resources/symbol2ensmusg.txt", sep = "\t", header = F, stringsAsFactors = F, col.names = c("msymbol", "ENSMUSG"))
m2h <- read.table("resources/ensmusg2ensg.txt", sep = "\t", header = F, stringsAsFactors = F, col.names = c("ENSMUSG", "ENSG"))
s2h <- merge(s2m, m2h, by = "ENSMUSG")
# Remove one-to-many homologs
one2many <- unique(s2h[duplicated(s2h$ENSG),]$ENSG)
s2h.one2one <- s2h %>%
  dplyr::filter(ENSG %ni% one2many)
# Subset to concise gene set 
keep <- read.table("resources/gene_annot.txt", sep = "\t", header = T, stringsAsFactors = F, col.names = c("ENSG", "chr", "start", "end"))
s2h.keep <- merge(s2h.one2one, keep, by = "ENSG")
# Apply to expression matrix
rowkeep <- row.names(mouse_immune) %in% s2h.keep$msymbol
mouse_immune <- mouse_immune[rowkeep,]
row.names(mouse_immune) <- s2h.keep[match(row.names(mouse_immune), s2h.keep$msymbol),]$ENSG

#' Filter, normalize, and scale data
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
mouse_immune.so <- CreateSeuratObject(raw.data = mouse_immune, project = "mouse_immune", min.cells = 5) #, meta.data = mouse_immune.annot)
mouse_immune.so <- FilterCells(mouse_immune.so, subset.names = "nGene", low.thresholds = 500, high.thresholds = Inf)
mouse_immune.so <- NormalizeData(mouse_immune.so) # log normalize, scale by library size
mouse_immune.so <- ScaleData(mouse_immune.so, min.cells.to.block = 1, block.size = 500) # center, scale variance within expression bins

#' Identify overdispersed genes
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
pdf("features/mouse_immune/variablegenes.pdf")
mouse_immune.so <- FindVariableGenes(mouse_immune.so, do.plot = TRUE, do.text = FALSE, do.contour = FALSE, cex.use = 0.1)
dev.off()

#' Perform PCA by way of partial SVD
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
mouse_immune.so <- RunPCA(mouse_immune.so, do.print = FALSE, pcs.compute = 20)
mouse_immune.so <- ProjectPCA(mouse_immune.so, do.print = FALSE)
#RunUMAP(mouse_immune.so, reduction.use = "pca", n_neighbors = 30L, min_dist = 0.3)

#' Cluster cells in PC space
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
mouse_immune.so <- FindClusters(object = mouse_immune.so, reduction.type = "pca", k.param = 5, dims.use = 1:5, save.SNN = TRUE, resolution = 1.0, force.recalc = TRUE)

#' Plot PC space to see clusters and other annotations
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
pdf("features/mouse_immune/PCA.pdf")
DimPlot(object = mouse_immune.so, reduction.use = "pca", cols.use = jdb_palette("lawhoops"), pt.size = 0.5, group.by = "ident")
DimPlot(object = mouse_immune.so, reduction.use = "pca", dim.1 = 3, dim.2 = 4, cols.use = jdb_palette("lawhoops"), pt.size = 0.5, group.by = "ident")
dev.off()

#' Write out projected gene loadings
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
mouse_immune.u <- GetGeneLoadings(mouse_immune.so, use.full = TRUE)
write.table(mouse_immune.u, "features/mouse_immune/u_matrix.txt", quote = F, row.names = T, col.names = T, sep = "\t")

#' Write out normalized expression across specified clusters
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
mouse_immune.so@meta.data$celltype <- row.names(mouse_immune.so@meta.data)
mouse_immune.so <- SetAllIdent(mouse_immune.so, id = "celltype")
mouse_immune.ae <- AverageExpression(mouse_immune.so, use.scale = TRUE)
write.table(mouse_immune.ae, "features/mouse_immune/ave_expr.txt", quote = F, row.names = T, col.names = T, sep = "\t")
