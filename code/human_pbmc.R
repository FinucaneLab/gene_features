#' ---
#' title: "Create gene features from human_pbmc"
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
human_pbmc <- Read10X("data/human_pbmc")
colnames(human_pbmc) <- gsub("_Cer_", "_", colnames(human_pbmc))
human_pbmc.annot <- read.table(text = colnames(human_pbmc), header = FALSE, sep = "_", col.names = c("annot1", "annot2", "barcode"))
row.names(human_pbmc.annot) <- colnames(human_pbmc)

#' Subset to concise gene set
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
# Subset to concise gene set 
s2h <- read.table("resources/ensg2symbol.txt", sep = "\t", header = F, stringsAsFactors = F, col.names = c("ENSG", "hsymbol")) 
keep <- read.table("resources/gene_annot.txt", sep = "\t", header = T, stringsAsFactors = F, col.names = c("ENSG", "chr", "start", "end"))
s2h.keep <- merge(s2h, keep, by = "ENSG")
# Apply to expression matrix
rowkeep <- row.names(human_pbmc) %in% s2h.keep$hsymbol
human_pbmc <- human_pbmc[rowkeep,]
row.names(human_pbmc) <- keep[match(row.names(human_pbmc), s2h.keep$hsymbol),]$ENSG

#' Filter, normalize, and scale data
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
human_pbmc.so <- CreateSeuratObject(raw.data = human_pbmc, project = "human_pbmc", min.cells = 5, meta.data = human_pbmc.annot)
human_pbmc.so <- FilterCells(human_pbmc.so, subset.names = "nGene", low.thresholds = 500, high.thresholds = Inf)
human_pbmc.so <- NormalizeData(human_pbmc.so) # log normalize, scale by library size
human_pbmc.so <- ScaleData(human_pbmc.so, min.cells.to.block = 1, block.size = 500) # center, scale variance within expression bins

#' Identify overdispersed genes
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
pdf("features/human_pbmc/variablegenes.pdf")
human_pbmc.so <- FindVariableGenes(human_pbmc.so, do.plot = TRUE, do.text = FALSE, do.contour = FALSE, cex.use = 0.1, num.bin = 40, y.cutoff = 0.5, x.high.cutoff = 10)
dev.off()

#' Perform PCA by way of partial SVD
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
human_pbmc.so <- RunPCA(human_pbmc.so, do.print = FALSE, pcs.compute = 100)
human_pbmc.so <- ProjectPCA(human_pbmc.so, do.print = FALSE)
#RunUMAP(human_pbmc.so, reduction.use = "pca", n_neighbors = 30L, min_dist = 0.3)

#' Cluster cells in PC space
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
human_pbmc.so <- FindClusters(object = human_pbmc.so, reduction.type = "pca", k.param = 20, dims.use = 1:20, save.SNN = TRUE, resolution = 0.6, force.recalc = TRUE)

#' Plot PC space to see clusters and other annotations
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
pdf("features/human_pbmc/PCA.pdf")
DimPlot(object = human_pbmc.so, reduction.use = "pca", cols.use = jdb_palette("lawhoops"), pt.size = 0.5, group.by = "ident")
DimPlot(object = human_pbmc.so, reduction.use = "pca", cols.use = jdb_palette("lawhoops"), dim.1 = 3, dim.2 = 4, pt.size = 0.5, group.by = "ident")
DimPlot(object = human_pbmc.so, reduction.use = "pca", pt.size = 0.5, group.by = "annot1")
DimPlot(object = human_pbmc.so, reduction.use = "pca", dim.1 = 3, dim.2 = 4, pt.size = 0.5, group.by = "annot1")
DimPlot(object = human_pbmc.so, reduction.use = "pca", pt.size = 0.5, group.by = "annot2")
DimPlot(object = human_pbmc.so, reduction.use = "pca", dim.1 = 3, dim.2 = 4, pt.size = 0.5, group.by = "annot2")
dev.off()

#' Write out projected gene loadings
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
human_pbmc.u <- GetGeneLoadings(human_pbmc.so, use.full = TRUE)
write.table(human_pbmc.u, "features/human_pbmc/u_matrix.txt", quote = F, row.names = T, col.names = T, sep = "\t")

#' Write out normalized expression across specified clusters
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
human_pbmc.so <- SetAllIdent(human_pbmc.so, id = "annot1")
human_pbmc.ae <- AverageExpression(human_pbmc.so, use.scale = TRUE)
write.table(human_pbmc.ae, "features/human_pbmc/ave_expr.txt", quote = F, row.names = T, col.names = T, sep = "\t")