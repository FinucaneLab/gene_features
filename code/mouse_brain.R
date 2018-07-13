#' ---
#' title: "Create gene features from mouse_brain"
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
# Saved as RDS objects already
mouse_brain <- readRDS("data/mouse_brain/metacells.BrainCellAtlas_Saunders_version_2018.04.01.RDS")
mouse_brain.annot <- readRDS("data/mouse_brain/annotation.BrainCellAtlas_Saunders_version_2018.04.01.RDS")
mouse_brain.annot <- as.data.frame(mouse_brain.annot)
row.names(mouse_brain.annot) <- mouse_brain.annot$tissue_subcluster

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
rowkeep <- row.names(mouse_brain) %in% s2h.keep$msymbol
mouse_brain <- mouse_brain[rowkeep,]
row.names(mouse_brain) <- s2h.keep[match(row.names(mouse_brain), s2h.keep$msymbol),]$ENSG

#' Filter, normalize, and scale data
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
mouse_brain.so <- CreateSeuratObject(raw.data = mouse_brain, project = "mouse_brain", min.cells = 5, meta.data = mouse_brain.annot)
mouse_brain.so <- FilterCells(mouse_brain.so, subset.names = "nGene", low.thresholds = 500, high.thresholds = Inf)
mouse_brain.so <- NormalizeData(mouse_brain.so) # log normalize, scale by library size
mouse_brain.so <- ScaleData(mouse_brain.so, min.cells.to.block = 1, block.size = 500) # center, scale variance within expression bins

#' Identify overdispersed genes
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
pdf("features/mouse_brain/variablegenes.pdf")
mouse_brain.so <- FindVariableGenes(mouse_brain.so, do.plot = TRUE, do.text = FALSE, do.contour = FALSE, cex.use = 0.1)
dev.off()

#' Perform PCA by way of partial SVD
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
mouse_brain.so <- RunPCA(mouse_brain.so, do.print = FALSE, pcs.compute = 100)
mouse_brain.so <- ProjectPCA(mouse_brain.so, do.print = FALSE)
#RunUMAP(mouse_brain.so, reduction.use = "pca", n_neighbors = 30L, min_dist = 0.3)

#' Cluster cells in PC space
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
mouse_brain.so <- FindClusters(object = mouse_brain.so, reduction.type = "pca", k.param = 20, dims.use = 1:20, save.SNN = TRUE, resolution = 0.6, force.recalc = TRUE)

#' Plot PC space to see clusters and other annotations
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
pdf("features/mouse_brain/PCA.pdf")
DimPlot(object = mouse_brain.so, reduction.use = "pca", cols.use = jdb_palette("lawhoops"), pt.size = 0.5, group.by = "ident")
DimPlot(object = mouse_brain.so, reduction.use = "pca", dim.1 = 3, dim.2 = 4, cols.use = jdb_palette("lawhoops"), pt.size = 0.5, group.by = "ident")
DimPlot(object = mouse_brain.so, reduction.use = "pca", cols.use = jdb_palette("lawhoops"), pt.size = 0.5, group.by = "class")
DimPlot(object = mouse_brain.so, reduction.use = "pca", dim.1 = 3, dim.2 = 4, cols.use = jdb_palette("lawhoops"), pt.size = 0.5, group.by = "class")
DimPlot(object = mouse_brain.so, reduction.use = "pca", cols.use = jdb_palette("lawhoops"), pt.size = 0.5, group.by = "tissue")
DimPlot(object = mouse_brain.so, reduction.use = "pca", dim.1 = 3, dim.2 = 4, cols.use = jdb_palette("lawhoops"), pt.size = 0.5, group.by = "tissue")
dev.off()

#' Write out projected gene loadings
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
mouse_brain.u <- GetGeneLoadings(mouse_brain.so, use.full = TRUE)
write.table(mouse_brain.u, "features/mouse_brain/u_matrix.txt", quote = F, row.names = T, col.names = T, sep = "\t")

#' Write out normalized expression across specified clusters
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
mouse_brain.so <- SetAllIdent(mouse_brain.so, id = "class")
mouse_brain.ae <- AverageExpression(mouse_brain.so, use.scale = TRUE)
write.table(mouse_brain.ae, "features/mouse_brain/ave_expr.txt", quote = F, row.names = T, col.names = T, sep = "\t")


