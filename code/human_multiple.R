#' ---
#' title: "Create gene features from human_multiple"
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
human_multiple <- readRDS("data/human_multiple/GTEx.rds")
human_multiple.annot <- read.table("data/human_multiple/GTEx_samples.txt", header = T, sep = "\t")
human_multiple.annot <- as.data.frame(human_multiple.annot)
row.names(human_multiple.annot) <- gsub("-", ".", human_multiple.annot$SAMPID)

#' Subset to concise gene set
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
# Subset to concise gene set 
keep <- read.table("resources/gene_annot.txt", sep = "\t", header = T, stringsAsFactors = F, col.names = c("ENSG", "chr", "start", "end"))
# Apply to expression matrix
rowkeep <- row.names(human_multiple) %in% keep$ENSG
human_multiple <- human_multiple[rowkeep,]
row.names(human_multiple) <- keep[match(row.names(human_multiple), keep$ENSG),]$ENSG

#' Filter, normalize, and scale data
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
human_multiple.so <- CreateSeuratObject(raw.data = human_multiple, project = "human_multiple", min.cells = 5, meta.data = human_multiple.annot)
human_multiple.so <- FilterCells(human_multiple.so, subset.names = "nGene", low.thresholds = 500, high.thresholds = Inf)
human_multiple.so <- NormalizeData(human_multiple.so) # log normalize, scale by library size
human_multiple.so <- ScaleData(human_multiple.so, min.cells.to.block = 1, block.size = 500) # center, scale variance within expression bins

#' Identify overdispersed genes
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
pdf("features/human_multiple/variablegenes.pdf")
human_multiple.so <- FindVariableGenes(human_multiple.so, do.plot = TRUE, do.text = FALSE, do.contour = FALSE, cex.use = 0.1)
dev.off()

#' Perform PCA by way of partial SVD
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
human_multiple.so <- RunPCA(human_multiple.so, do.print = FALSE, pcs.compute = 100)
human_multiple.so <- ProjectPCA(human_multiple.so, do.print = FALSE)
#RunUMAP(human_multiple.so, reduction.use = "pca", n_neighbors = 30L, min_dist = 0.3)

#' Cluster cells in PC space
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
human_multiple.so <- FindClusters(object = human_multiple.so, reduction.type = "pca", k.param = 30, dims.use = 1:20, save.SNN = TRUE, resolution = 0.6, force.recalc = TRUE)

#' Plot PC space to see clusters and other annotations
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
pdf("features/human_multiple/PCA.pdf")
DimPlot(object = human_multiple.so, reduction.use = "pca", pt.size = 0.5, group.by = "ident")
DimPlot(object = human_multiple.so, reduction.use = "pca", dim.1 = 3, dim.2 = 4, pt.size = 0.5, group.by = "ident")
DimPlot(object = human_multiple.so, reduction.use = "pca", pt.size = 0.5, group.by = "SMTS")
DimPlot(object = human_multiple.so, reduction.use = "pca", dim.1 = 3, dim.2 = 4, pt.size = 0.5, group.by = "SMTS")
dev.off()

#' Write out projected gene loadings
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
human_multiple.u <- GetGeneLoadings(human_multiple.so, use.full = TRUE)
write.table(human_multiple.u, "features/human_multiple/u_matrix.txt", quote = F, row.names = T, col.names = T, sep = "\t")

#' Write out normalized expression across specified clusters
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
human_multiple.so <- SetAllIdent(human_multiple.so, id = "SMTS")
human_multiple.ae <- AverageExpression(human_multiple.so, use.scale = TRUE)
write.table(human_multiple.ae, "features/human_multiple/ave_expr.txt", quote = F, row.names = T, col.names = T, sep = "\t")
