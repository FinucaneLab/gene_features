#' ---
#' title: "Create gene features from mouse_aorta"
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
file <- "data/mouse_aorta/scRNAseq_aorta_counts.csv.gz"
mouse_aorta <- data.frame(fread(paste0("zcat < ", file), skip = 0, sep = ","))[,-1] %>% 
  data.matrix() %>% 
  Matrix()
rownames(mouse_aorta) <- data.frame(fread(paste0("zcat < ", file), select = 1:2, skip = 1, sep = ","))[,1] %>%
  gsub("\\..*", "", .)
#mouse_aorta.annot <- readRDS("data/mouse_aorta/")
#mouse_aorta.annot <- as.data.frame(mouse_aorta.annot)
#row.names(mouse_aorta.annot) <- mouse_aorta.annot$tissue_subcluster

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
rowkeep <- row.names(mouse_aorta) %in% s2h.keep$msymbol
mouse_aorta <- mouse_aorta[rowkeep,]
row.names(mouse_aorta) <- s2h.keep[match(row.names(mouse_aorta), s2h.keep$msymbol),]$ENSG
# Remove duplicates
rowkeep <- rownames(mouse_aorta) %ni% rownames(mouse_aorta)[duplicated(rownames(mouse_aorta))]
mouse_aorta <- mouse_aorta[rowkeep,]

#' Filter, normalize, and scale data
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
mouse_aorta.so <- CreateSeuratObject(raw.data = mouse_aorta, project = "mouse_aorta", min.cells = 5) #, meta.data = mouse_aorta.annot)
mouse_aorta.so <- FilterCells(mouse_aorta.so, subset.names = "nGene", low.thresholds = 500, high.thresholds = Inf)
mouse_aorta.so <- NormalizeData(mouse_aorta.so) # log normalize, scale by library size
mouse_aorta.so <- ScaleData(mouse_aorta.so, min.cells.to.block = 1, block.size = 500) # center, scale variance within expression bins

#' Identify overdispersed genes
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
pdf("features/mouse_aorta/variablegenes.pdf")
mouse_aorta.so <- FindVariableGenes(mouse_aorta.so, do.plot = TRUE, do.text = FALSE, do.contour = FALSE, cex.use = 0.1)
dev.off()

#' Perform PCA by way of partial SVD
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
mouse_aorta.so <- RunPCA(mouse_aorta.so, do.print = FALSE, pcs.compute = 100)
mouse_aorta.so <- ProjectPCA(mouse_aorta.so, do.print = FALSE)
#RunUMAP(mouse_aorta.so, reduction.use = "pca", n_neighbors = 30L, min_dist = 0.3)

#' Cluster cells in PC space
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
mouse_aorta.so <- FindClusters(object = mouse_aorta.so, reduction.type = "pca", k.param = 20, dims.use = 1:20, save.SNN = TRUE, resolution = 0.6, force.recalc = TRUE)

#' Plot PC space to see clusters and other annotations
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
pdf("features/mouse_aorta/PCA.pdf")
DimPlot(object = mouse_aorta.so, reduction.use = "pca", cols.use = jdb_palette("lawhoops"), pt.size = 0.5, group.by = "ident")
DimPlot(object = mouse_aorta.so, reduction.use = "pca", dim.1 = 3, dim.2 = 4, cols.use = jdb_palette("lawhoops"), pt.size = 0.5, group.by = "ident")
dev.off()

#' Write out projected gene loadings
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
mouse_aorta.u <- GetGeneLoadings(mouse_aorta.so, use.full = TRUE)
write.table(mouse_aorta.u, "features/mouse_aorta/u_matrix.txt", quote = F, row.names = T, col.names = T, sep = "\t")

#' Write out normalized expression across specified clusters
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = TRUE, eval = TRUE
#mouse_aorta.so <- SetAllIdent(mouse_aorta.so, id = "class")
mouse_aorta.ae <- AverageExpression(mouse_aorta.so, use.scale = TRUE)
write.table(mouse_aorta.ae, "features/mouse_aorta/ave_expr.txt", quote = F, row.names = T, col.names = T, sep = "\t")
