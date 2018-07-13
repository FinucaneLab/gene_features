#' ---
#' title: "Create gene features from human_pancreas"
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
file <- "data/human_pancreas/GSE85241_cellsystems_dataset_4donors_updated.csv.gz"
human_pancreas <- data.frame(fread(paste0("zcat < ", file), skip = 1, sep = "\t"))[,-1] %>% 
  data.matrix() %>% 
  Matrix()
rownames(human_pancreas) <- data.frame(fread(paste0("zcat < ", file), select = 1, skip = 1, sep = "\t"))[,1] %>%
  gsub("__.*", "", .)
file <- "data/human_pancreas/cellnames.txt"
colnames(human_pancreas) <- read.table(file, sep = "\t")[,1]
human_pancreas.annot <- as.data.frame(cbind(cell = colnames(human_pancreas), donor = gsub("-.*", "", colnames(human_pancreas)), library = gsub(".*-|_.*", "", colnames(human_pancreas))))
row.names(human_pancreas.annot) <- human_pancreas.annot$cell

#' Subset to concise gene set
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
# Subset to concise gene set 
s2h <- read.table("resources/ensg2symbol.txt", sep = "\t", header = F, stringsAsFactors = F, col.names = c("ENSG", "hsymbol")) 
keep <- read.table("resources/gene_annot.txt", sep = "\t", header = T, stringsAsFactors = F, col.names = c("ENSG", "chr", "start", "end"))
s2h.keep <- merge(s2h, keep, by = "ENSG")
# Apply to expression matrix
rowkeep <- row.names(human_pancreas) %in% s2h.keep$hsymbol
human_pancreas <- human_pancreas[rowkeep,]
row.names(human_pancreas) <- keep[match(row.names(human_pancreas), s2h.keep$hsymbol),]$ENSG

#' Filter, normalize, and scale data
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
human_pancreas.so <- CreateSeuratObject(raw.data = human_pancreas, project = "human_pancreas", min.cells = 5, meta.data = human_pancreas.annot)
human_pancreas.so <- FilterCells(human_pancreas.so, subset.names = "nGene", low.thresholds = 500, high.thresholds = Inf)
human_pancreas.so <- NormalizeData(human_pancreas.so) # log normalize, scale by library size
human_pancreas.so <- ScaleData(human_pancreas.so, min.cells.to.block = 1, block.size = 500) # center, scale variance within expression bins

#' Identify overdispersed genes
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
pdf("features/human_pancreas/variablegenes.pdf")
human_pancreas.so <- FindVariableGenes(human_pancreas.so, do.plot = TRUE, do.text = FALSE, do.contour = FALSE, cex.use = 0.1)
dev.off()

#' Perform PCA by way of partial SVD
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
human_pancreas.so <- RunPCA(human_pancreas.so, do.print = FALSE, pcs.compute = 100)
human_pancreas.so <- ProjectPCA(human_pancreas.so, do.print = FALSE)
#RunUMAP(human_pancreas.so, reduction.use = "pca", n_neighbors = 30L, min_dist = 0.3)

#' Cluster cells in PC space
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
human_pancreas.so <- FindClusters(object = human_pancreas.so, reduction.type = "pca", k.param = 20, dims.use = 1:20, save.SNN = TRUE, resolution = 0.6, force.recalc = TRUE)

#' Plot PC space to see clusters and other annotations
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
pdf("features/human_pancreas/PCA.pdf")
DimPlot(object = human_pancreas.so, reduction.use = "pca", pt.size = 0.5, group.by = "ident")
DimPlot(object = human_pancreas.so, reduction.use = "pca", dim.1 = 3, dim.2 = 4, pt.size = 0.5, group.by = "ident")
DimPlot(object = human_pancreas.so, reduction.use = "pca", pt.size = 0.5, group.by = "donor")
DimPlot(object = human_pancreas.so, reduction.use = "pca", dim.1 = 3, dim.2 = 4, pt.size = 0.5, group.by = "donor")
DimPlot(object = human_pancreas.so, reduction.use = "pca", pt.size = 0.5, group.by = "library")
DimPlot(object = human_pancreas.so, reduction.use = "pca", dim.1 = 3, dim.2 = 4, pt.size = 0.5, group.by = "library")
dev.off()

#' Write out projected gene loadings
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
human_pancreas.u <- GetGeneLoadings(human_pancreas.so, use.full = TRUE)
write.table(human_pancreas.u, "features/human_pancreas/u_matrix.txt", quote = F, row.names = T, col.names = T, sep = "\t")

#' Write out normalized expression across specified clusters
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
human_pancreas.ae <- AverageExpression(human_pancreas.so, use.scale = TRUE)
write.table(human_pancreas.ae, "features/human_pancreas/ave_expr.txt", quote = F, row.names = T, col.names = T, sep = "\t")
