#' ---
#' title: "Create gene features from human_heme"
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
file <- "data/human_heme/16populations_RNAcounts.txt"
human_heme <- data.frame(fread(file, sep = "\t"))[,-1] %>% 
  data.matrix() %>% 
  Matrix()
rownames(human_heme) <- data.frame(fread(file), select = 1, skip = 1, sep = "\t")[,1]

#' Subset to concise gene set
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
# Subset to concise gene set 
s2h <- read.table("resources/ensg2symbol.txt", sep = "\t", header = F, stringsAsFactors = F, col.names = c("ENSG", "hsymbol")) 
keep <- read.table("resources/gene_annot.txt", sep = "\t", header = T, stringsAsFactors = F, col.names = c("ENSG", "chr", "start", "end"))
s2h.keep <- merge(s2h, keep, by = "ENSG")
# Apply to expression matrix
rowkeep <- row.names(human_heme) %in% s2h.keep$hsymbol
human_heme <- human_heme[rowkeep,]
row.names(human_heme) <- keep[match(row.names(human_heme), s2h.keep$hsymbol),]$ENSG
# Remove duplicates
rowkeep <- rownames(human_heme) %ni% rownames(human_heme)[duplicated(rownames(human_heme))]
human_heme <- human_heme[rowkeep,]

#' Filter, normalize, and scale data
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
human_heme.so <- CreateSeuratObject(raw.data = human_heme, project = "human_heme", min.cells = 5)
human_heme.so <- FilterCells(human_heme.so, subset.names = "nGene", low.thresholds = 500, high.thresholds = Inf)
human_heme.so <- NormalizeData(human_heme.so) # log normalize, scale by library size
human_heme.so <- ScaleData(human_heme.so, min.cells.to.block = 1, block.size = 500) # center, scale variance within expression bins

#' Identify overdispersed genes
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
pdf("features/human_heme/variablegenes.pdf")
human_heme.so <- FindVariableGenes(human_heme.so, do.plot = TRUE, do.text = FALSE, do.contour = FALSE, cex.use = 0.1)
dev.off()

#' Perform PCA by way of partial SVD
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
human_heme.so <- RunPCA(human_heme.so, do.print = FALSE, pcs.compute = 6)
human_heme.so <- ProjectPCA(human_heme.so, do.print = FALSE)
#RunUMAP(human_heme.so, reduction.use = "pca", n_neighbors = 30L, min_dist = 0.3)

#' Cluster cells in PC space
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
human_heme.so <- FindClusters(object = human_heme.so, reduction.type = "pca", k.param = 5, dims.use = 1:5, save.SNN = TRUE, resolution = 2.0, force.recalc = TRUE)

#' Plot PC space to see clusters and other annotations
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
pdf("features/human_heme/PCA.pdf")
human_heme.so@meta.data$type <- row.names(human_heme.so@meta.data)
DimPlot(object = human_heme.so, reduction.use = "pca", cols.use = jdb_palette("lawhoops"), pt.size = 0.5, group.by = "ident")
DimPlot(object = human_heme.so, reduction.use = "pca", cols.use = jdb_palette("lawhoops"), dim.1 = 3, dim.2 = 4, pt.size = 0.5, group.by = "ident")
DimPlot(object = human_heme.so, reduction.use = "pca", cols.use = jdb_palette("lawhoops"), pt.size = 0.5, group.by = "type")
DimPlot(object = human_heme.so, reduction.use = "pca", cols.use = jdb_palette("lawhoops"), dim.1 = 3, dim.2 = 4, pt.size = 0.5, group.by = "type")
dev.off()

#' Write out projected gene loadings
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
human_heme.u <- GetGeneLoadings(human_heme.so, use.full = TRUE)
write.table(human_heme.u, "features/human_heme/u_matrix.txt", quote = F, row.names = T, col.names = T, sep = "\t")

#' Write out normalized expression across specified clusters
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
human_heme.so <- SetAllIdent(human_heme.so, id = "type")
human_heme.ae <- AverageExpression(human_heme.so, use.scale = TRUE)
write.table(human_heme.ae, "features/human_heme/ave_expr.txt", quote = F, row.names = T, col.names = T, sep = "\t")
