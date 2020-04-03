#' ---
#' title: "Create gene features from human_immune"
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
library(cellrangerRkit)
source("code/utils.R")
#library(reticulate)
#use_virtualenv("/data/aryee/julirsch/python/venv3/bin/activate", required = TRUE)

#' Read in data and annotations
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
human_immune <- get_matrix_from_h5("data/human_immume/ica_bone_marrow_h5.h5", genome = "GRCh38")
human_immune.annot <- cbind(donor = gsub("Manton", "", gsub("_.*", "", colnames(human_immune))), library = gsub("-.*", "", gsub(".*HiSeq_", "", colnames(human_immune))))
row.names(human_immune.annot) <- colnames(human_immune)

#' Subset to concise gene set
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
# Subset to concise gene set 
s2h <- read.table("resources/ensg2symbol.txt", sep = "\t", header = F, stringsAsFactors = F, col.names = c("ENSG", "hsymbol")) 
keep <- read.table("resources/gene_annot.txt", sep = "\t", header = T, stringsAsFactors = F, col.names = c("ENSG", "chr", "start", "end"))
s2h.keep <- merge(s2h, keep, by = "ENSG")
# Apply to expression matrix
rowkeep <- row.names(human_immune) %in% s2h.keep$ENSG
human_immune <- human_immune[rowkeep,]
row.names(human_immune) <- keep[match(row.names(human_immune), s2h.keep$hsymbol),]$ENSG

#' Filter, normalize, and scale data
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
human_immune.so <- CreateSeuratObject(raw.data = exprs(human_immune), project = "human_immune", min.cells = 5)
human_immune.so <- FilterCells(human_immune.so, subset.names = "nGene", low.thresholds = 500, high.thresholds = Inf)
human_immune.so <- NormalizeData(human_immune.so) # log normalize, scale by library size
human_immune.so <- ScaleData(human_immune.so, min.cells.to.block = 1, block.size = 500) # center, scale variance within expression bins

#' Identify overdispersed genes
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
pdf("features/human_immune/variablegenes.pdf")
human_immune.so <- FindVariableGenes(human_immune.so, do.plot = TRUE, do.text = FALSE, do.contour = FALSE, cex.use = 0.1, num.bin = 40, y.cutoff = 0.25, x.high.cutoff = 30, x.low.cutoff = 0.05)
dev.off()

#' Perform PCA by way of partial SVD
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
human_immune.so <- RunPCA(human_immune.so, do.print = FALSE, pcs.compute = 100)
human_immune.so <- ProjectPCA(human_immune.so, do.print = FALSE)
#RunUMAP(human_immune.so, reduction.use = "pca", n_neighbors = 30L, min_dist = 0.3)

#' Cluster cells in PC space
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
human_immune.so <- FindClusters(object = human_immune.so, reduction.type = "pca", k.param = 20, dims.use = 1:20, save.SNN = FALSE, resolution = 0.6, force.recalc = TRUE)

#' Plot PC space to see clusters and other annotations
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
pdf("features/human_immune/PCA.pdf")
DimPlot(object = human_immune.so, reduction.use = "pca", pt.size = 0.5, group.by = "ident")
DimPlot(object = human_immune.so, reduction.use = "pca", dim.1 = 3, dim.2 = 4, pt.size = 0.5, group.by = "ident")
DimPlot(object = human_immune.so, reduction.use = "pca", cols.use = jdb_palette("lawhoops"), pt.size = 0.5, group.by = "orig.ident")
DimPlot(object = human_immune.so, reduction.use = "pca", cols.use = jdb_palette("lawhoops"), dim.1 = 3, dim.2 = 4, pt.size = 0.5, group.by = "orig.ident")
dev.off()

#' Write out projected gene loadings
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
human_immune.u <- GetGeneLoadings(human_immune.so, use.full = TRUE)
write.table(human_immune.u, "features/human_immune/u_matrix.txt", quote = F, row.names = T, col.names = T, sep = "\t")

#' Write out normalized expression across specified clusters
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
human_immune.so <- SetAllIdent(human_immune.so, id = "celltype")
human_immune.ae <- AverageExpression(human_immune.so, use.scale = TRUE)
write.table(human_immune.ae, "features/human_immune/ave_expr.txt", quote = F, row.names = T, col.names = T, sep = "\t")
