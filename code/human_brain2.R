#' ---
#' title: "Create gene features from human_brain2"
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
file <- "data/human_brain2/RNAseq/expression_matrix.csv"
human_brain2 <- data.frame(fread(file, sep = ","))[,-1] %>% 
  data.matrix() %>% 
  Matrix()
file <- "data/human_brain2/RNAseq/rows_metadata.csv"
rownames(human_brain2) <- data.frame(fread(file, select = 3, skip = 1, sep = ","))[,1]
file <- "data/human_brain2/RNAseq/columns_metadata.csv"
human_brain2.annot <- read.table(file, sep = ",", header = TRUE)
colnames(human_brain2) <- paste0("V", seq(1:dim(human_brain2)[2]))
row.names(human_brain2.annot) <- colnames(human_brain2)

#' Subset to concise gene set
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
# Subset to concise gene set 
s2h <- read.table("resources/ensg2symbol.txt", sep = "\t", header = F, stringsAsFactors = F, col.names = c("ENSG", "hsymbol")) 
keep <- read.table("resources/gene_annot.txt", sep = "\t", header = T, stringsAsFactors = F, col.names = c("ENSG", "chr", "start", "end"))
s2h.keep <- merge(s2h, keep, by = "ENSG")
# Apply to expression matrix
rowkeep <- row.names(human_brain2) %in% s2h.keep$ENSG
human_brain2 <- human_brain2[rowkeep,]

#' Filter, normalize, and scale data
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
human_brain2.so <- CreateSeuratObject(raw.data = human_brain2, project = "human_brain2", min.cells = 5, meta.data = human_brain2.annot)
human_brain2.so <- FilterCells(human_brain2.so, subset.names = "nGene", low.thresholds = 500, high.thresholds = Inf)
human_brain2.so <- NormalizeData(human_brain2.so) # log normalize, scale by library size
human_brain2.so <- ScaleData(human_brain2.so, min.cells.to.block = 1, block.size = 500) # center, scale variance within expression bins

#' Identify overdispersed genes
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
pdf("features/human_brain2/variablegenes.pdf")
human_brain2.so <- FindVariableGenes(human_brain2.so, do.plot = TRUE, do.text = FALSE, do.contour = FALSE, cex.use = 0.1)
dev.off()

#' Perform PCA by way of partial SVD
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
human_brain2.so <- RunPCA(human_brain2.so, do.print = FALSE, pcs.compute = 100)
human_brain2.so <- ProjectPCA(human_brain2.so, do.print = FALSE)
#RunUMAP(human_brain2.so, reduction.use = "pca", n_neighbors = 30L, min_dist = 0.3)

#' Cluster cells in PC space
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
human_brain2.so <- FindClusters(object = human_brain2.so, reduction.type = "pca", k.param = 20, dims.use = 1:20, save.SNN = TRUE, resolution = 1.0, force.recalc = TRUE)

#' Plot PC space to see clusters and other annotations
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
pdf("features/human_brain2/PCA.pdf")
DimPlot(object = human_brain2.so, reduction.use = "pca", cols.use = jdb_palette("lawhoops"), pt.size = 0.5, group.by = "ident")
DimPlot(object = human_brain2.so, reduction.use = "pca", cols.use = jdb_palette("lawhoops"), dim.1 = 3, dim.2 = 4, pt.size = 0.5, group.by = "ident")
DimPlot(object = human_brain2.so, reduction.use = "pca", pt.size = 0.5, group.by = "donor_id")
DimPlot(object = human_brain2.so, reduction.use = "pca", dim.1 = 3, dim.2 = 4, pt.size = 0.5, group.by = "donor_id")
DimPlot(object = human_brain2.so, reduction.use = "pca", pt.size = 0.5, group.by = "structure_acronym")
DimPlot(object = human_brain2.so, reduction.use = "pca", dim.1 = 3, dim.2 = 4, pt.size = 0.5, group.by = "structure_acronym")
dev.off()

#' Write out projected gene loadings
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
human_brain2.u <- GetGeneLoadings(human_brain2.so, use.full = TRUE)
write.table(human_brain2.u, "features/human_brain2/u_matrix.txt", quote = F, row.names = T, col.names = T, sep = "\t")

#' Write out normalized expression across specified clusters
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
human_brain2.so <- SetAllIdent(human_brain2.so, id = "structure_acronym")
human_brain2.ae <- AverageExpression(human_brain2.so, use.scale = TRUE)
write.table(human_brain2.ae, "features/human_brain2/ave_expr.txt", quote = F, row.names = T, col.names = T, sep = "\t")
