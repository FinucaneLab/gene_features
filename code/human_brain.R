#' ---
#' title: "Create gene features from human_brain"
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
file <- "data/human_brain/GSE97930_CerebellarHem_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt.gz"
ch <- data.frame(fread(paste0("zcat < ", file), skip = 1, sep = "\t"))[,-1] %>% 
  data.matrix() %>% 
  Matrix()
rownames(ch) <- data.frame(fread(paste0("zcat < ", file), select = 1, skip = 1, sep = "\t"))[,1]
colnames(ch) <- colnames(fread(paste0("zcat < ", file), skip = 0, nrows = 1, sep = "\t", header = TRUE, fill= TRUE))[-(dim(ch)[2] + 1)]
file <- "data/human_brain/GSE97930_FrontalCortex_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt.gz"
fc <- data.frame(fread(paste0("zcat < ", file), skip = 1, sep = "\t"))[,-1] %>% 
  data.matrix() %>% 
  Matrix()
rownames(fc) <- data.frame(fread(paste0("zcat < ", file), select = 1, skip = 1, sep = "\t"))[,1]
colnames(fc) <- colnames(fread(paste0("zcat < ", file), skip = 0, nrows = 1, sep = "\t", header = TRUE, fill= TRUE))[-(dim(fc)[2] + 1)]
file <- "data/human_brain/GSE97930_VisualCortex_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt.gz"
vc <- data.frame(fread(paste0("zcat < ", file), skip = 1, sep = "\t"))[,-1] %>% 
  data.matrix() %>% 
  Matrix()
rownames(vc) <- data.frame(fread(paste0("zcat < ", file), select = 1, skip = 1, sep = "\t"))[,1]
colnames(vc) <- colnames(fread(paste0("zcat < ", file), skip = 0, nrows = 1, sep = "\t", header = TRUE, fill= TRUE))[-(dim(vc)[2] + 1)]
genes.int <- intersect(intersect(rownames(ch), rownames(fc)), rownames(vc))
length(genes.int)
human_brain <- cbind(ch[genes.int,], fc[genes.int,], vc[genes.int,])
rm(ch); rm(fc); rm(vc)
colnames(human_brain) <- gsub("_Cer_", "_", colnames(human_brain))
human_brain.annot <- read.table(text = colnames(human_brain), header = FALSE, sep = "_", col.names = c("annot1", "annot2", "barcode"))
row.names(human_brain.annot) <- colnames(human_brain)

#' Subset to concise gene set
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
# Subset to concise gene set 
s2h <- read.table("resources/ensg2symbol.txt", sep = "\t", header = F, stringsAsFactors = F, col.names = c("ENSG", "hsymbol")) 
keep <- read.table("resources/gene_annot.txt", sep = "\t", header = T, stringsAsFactors = F, col.names = c("ENSG", "chr", "start", "end"))
s2h.keep <- merge(s2h, keep, by = "ENSG")
# Apply to expression matrix
rowkeep <- row.names(human_brain) %in% s2h.keep$hsymbol
human_brain <- human_brain[rowkeep,]
row.names(human_brain) <- keep[match(row.names(human_brain), s2h.keep$hsymbol),]$ENSG

#' Filter, normalize, and scale data
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
human_brain.so <- CreateSeuratObject(raw.data = human_brain, project = "human_brain", min.cells = 5, meta.data = human_brain.annot)
human_brain.so <- FilterCells(human_brain.so, subset.names = "nGene", low.thresholds = 500, high.thresholds = Inf)
human_brain.so <- NormalizeData(human_brain.so) # log normalize, scale by library size
human_brain.so <- ScaleData(human_brain.so, min.cells.to.block = 1, block.size = 500) # center, scale variance within expression bins

#' Identify overdispersed genes
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
pdf("features/human_brain/variablegenes.pdf")
human_brain.so <- FindVariableGenes(human_brain.so, do.plot = TRUE, do.text = FALSE, do.contour = FALSE, cex.use = 0.1, num.bin = 40, y.cutoff = 0.5, x.high.cutoff = 20)
dev.off()

#' Perform PCA by way of partial SVD
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
human_brain.so <- RunPCA(human_brain.so, do.print = FALSE, pcs.compute = 100)
human_brain.so <- ProjectPCA(human_brain.so, do.print = FALSE)
#RunUMAP(human_brain.so, reduction.use = "pca", n_neighbors = 30L, min_dist = 0.3)

#' Cluster cells in PC space
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
human_brain.so <- FindClusters(object = human_brain.so, reduction.type = "pca", k.param = 20, dims.use = 1:20, save.SNN = TRUE, resolution = 0.6, force.recalc = TRUE)

#' Plot PC space to see clusters and other annotations
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
pdf("features/human_brain/PCA.pdf")
DimPlot(object = human_brain.so, reduction.use = "pca", cols.use = jdb_palette("lawhoops"), pt.size = 0.5, group.by = "ident")
DimPlot(object = human_brain.so, reduction.use = "pca", cols.use = jdb_palette("lawhoops"), dim.1 = 3, dim.2 = 4, pt.size = 0.5, group.by = "ident")
DimPlot(object = human_brain.so, reduction.use = "pca", cols.use = jdb_palette("lawhoops"), pt.size = 0.5, group.by = "celltype")
DimPlot(object = human_brain.so, reduction.use = "pca", cols.use = jdb_palette("lawhoops"), dim.1 = 3, dim.2 = 4, pt.size = 0.5, group.by = "celltype")
dev.off()

#' Write out projected gene loadings
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
human_brain.u <- GetGeneLoadings(human_brain.so, use.full = TRUE)
write.table(human_brain.u, "features/human_brain/u_matrix.txt", quote = F, row.names = T, col.names = T, sep = "\t")

#' Write out normalized expression across specified clusters
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
human_brain.so <- SetAllIdent(human_brain.so, id = "celltype")
human_brain.ae <- AverageExpression(human_brain.so, use.scale = TRUE)
write.table(human_brain.ae, "features/human_brain/ave_expr.txt", quote = F, row.names = T, col.names = T, sep = "\t")
