#' ---
#' title: "Create gene features from human_gut"
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
file <- "data/human_gut/Epi.250.dge.txt.gz"
human_gut <- data.frame(fread(paste0("zcat < ", file), skip = 0, sep = "\t"))[,-1] %>% 
  data.matrix() %>% 
  Matrix()
rownames(human_gut) <- data.frame(fread(paste0("zcat < ", file), select = 1, skip = 1, sep = "\t"))[,1]
human_gut.annot.meta <- fread("data/human_gut/all.anno.txt", sep = "\t", skip = 1, header = F, col.names = c("ident", "name", "description", "group"))
human_gut.annot.ident <- fread("data/human_gut/all.ident.txt", sep = "\t", skip = 1, header = F, col.names = c("cell", "ident"))
human_gut.annot <- merge(human_gut.annot.ident, human_gut.annot.meta, by = "ident")
human_gut.annot <- as.data.frame(human_gut.annot)
row.names(human_gut.annot) <- human_gut.annot$cell
human_gut <- human_gut[,colnames(human_gut) %in% human_gut.annot$cell]
human_gut.annot <- human_gut.annot[colnames(human_gut),]

#' Subset to concise gene set
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
# Subset to concise gene set 
s2h <- read.table("resources/ensg2symbol.txt", sep = "\t", header = F, stringsAsFactors = F, col.names = c("ENSG", "hsymbol")) 
keep <- read.table("resources/gene_annot.txt", sep = "\t", header = T, stringsAsFactors = F, col.names = c("ENSG", "chr", "start", "end"))
s2h.keep <- merge(s2h, keep, by = "ENSG")
# Apply to expression matrix
rowkeep <- row.names(human_gut) %in% s2h.keep$hsymbol
human_gut <- human_gut[rowkeep,]
row.names(human_gut) <- keep[match(row.names(human_gut), s2h.keep$hsymbol),]$ENSG

#' Filter, normalize, and scale data
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
human_gut.so <- CreateSeuratObject(raw.data = human_gut, project = "human_gut", min.cells = 5, meta.data = human_gut.annot)
human_gut.so <- FilterCells(human_gut.so, subset.names = "nGene", low.thresholds = 500, high.thresholds = Inf)
human_gut.so <- NormalizeData(human_gut.so) # log normalize, scale by library size
human_gut.so <- ScaleData(human_gut.so, min.cells.to.block = 1, block.size = 500) # center, scale variance within expression bins

#' Identify overdispersed genes
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
pdf("features/human_gut/variablegenes.pdf")
human_gut.so <- FindVariableGenes(human_gut.so, do.plot = TRUE, do.text = FALSE, do.contour = FALSE, cex.use = 0.1, num.bin = 40, y.cutoff = 0.5, x.high.cutoff = 20)
dev.off()

#' Perform PCA by way of partial SVD
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
human_gut.so <- RunPCA(human_gut.so, do.print = FALSE, pcs.compute = 100)
human_gut.so <- ProjectPCA(human_gut.so, do.print = FALSE)
#RunUMAP(human_gut.so, reduction.use = "pca", n_neighbors = 30L, min_dist = 0.3)

#' Cluster cells in PC space
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
human_gut.so <- FindClusters(object = human_gut.so, reduction.type = "pca", k.param = 20, dims.use = 1:20, save.SNN = TRUE, resolution = 0.6, force.recalc = TRUE)

#' Plot PC space to see clusters and other annotations
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
pdf("features/human_gut/PCA.pdf")
DimPlot(object = human_gut.so, reduction.use = "pca", cols.use = jdb_palette("lawhoops"), pt.size = 0.5, group.by = "ident")
DimPlot(object = human_gut.so, reduction.use = "pca", cols.use = jdb_palette("lawhoops"), dim.1 = 3, dim.2 = 4, pt.size = 0.5, group.by = "ident")
DimPlot(object = human_gut.so, reduction.use = "pca", cols.use = jdb_palette("lawhoops"), pt.size = 0.5, group.by = "celltype")
DimPlot(object = human_gut.so, reduction.use = "pca", cols.use = jdb_palette("lawhoops"), dim.1 = 3, dim.2 = 4, pt.size = 0.5, group.by = "celltype")
dev.off()

#' Write out projected gene loadings
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
human_gut.u <- GetGeneLoadings(human_gut.so, use.full = TRUE)
write.table(human_gut.u, "features/human_gut/u_matrix.txt", quote = F, row.names = T, col.names = T, sep = "\t")

#' Write out normalized expression across specified clusters
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
human_gut.so <- SetAllIdent(human_gut.so, id = "celltype")
human_gut.ae <- AverageExpression(human_gut.so, use.scale = TRUE)
write.table(human_gut.ae, "features/human_gut/ave_expr.txt", quote = F, row.names = T, col.names = T, sep = "\t")

