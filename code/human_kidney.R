#------------------------------------------------------SETUP-----------------------------------------------------#

# Load libraries
library(tidyverse)
library(data.table)
library(BuenColors)
library(Seurat)
library(irlba)
library(Matrix)
library(future)
library(reticulate)
library(ggrastr)
library(tidytext)
library(matrixTests)
source("utils.R")

# Set up parallelization
# Remember to use htop to delete forgotten forks
Sys.setenv(R_FUTURE_FORK_ENABLE = T)
options(future.globals.maxSize = 4096 * 1024^2)
plan(strategy = "multicore", workers = 32)

# Parameters
name <- "human_kidney"
number_pcs <- 30
vargenes <- 1600
clus_res <- 0.2

# Setup
dir.create(paste0("../plots/", name))
dir.create(paste0("../features/", name))

# Notes on data:
# We observe strong batch effects, so we run the integration pipeline.

#------------------------------------------------LOAD AND FORMAT DATA-----------------------------------------------#

# Load this in in case we need it later
keep <- read.table("../resources/gene_annot_jun10.txt", sep = "\t", header = T, stringsAsFactors = F, col.names = c("ENSG", "symbol", "chr", "start", "end", "TSS"))

# Read in data and annotations
# We'll have to read in each batch and cbind
read_sparse_mat <- function(x) {
  file <- paste0("../data/", name, "/", x, "matrix.mtx.gz")
  mat <- readMM(file) %>%
    data.matrix() %>%
    Matrix(sparse = TRUE)
  file <- paste0("../data/", name, "/", x, "features.tsv.gz")
  rows <- data.frame(fread(file, header=F), row.names=1)
  file <- paste0("../data/", name, "/", x, "barcodes.tsv.gz")
  cols <- data.frame(fread(file, header=F), row.names=1)
  rownames(mat) <- rownames(rows)
  colnames(mat) <- rownames(cols)
  mat <- ConvertToENSGAndProcessMatrix(mat, "ENSG")
  ### Reorder according to keep so that we can cbind later
  mat <- mat[match(keep$ENSG, rownames(mat)),]
  return(mat)
}

filelist = c("GSM4145204_kidney1_", "GSM4145205_kidney2_", "GSM4145206_kidney3_")
# Read all matrices into a list
datalist <- lapply(filelist, read_sparse_mat)

# Bind assuming same row order
mat <- do.call("cbind", datalist)
colnames(mat) <- make.names(colnames(mat), unique=T)

# Make batch indicators
mat.annot <- matrix(0, nrow = dim(mat)[2], ncol = 1)
curr_batch_start_ind = 1
for (i in 1:length(datalist)) {
  num_cells <- dim(datalist[[i]])[2]
  mat.annot[curr_batch_start_ind:(curr_batch_start_ind + num_cells - 1), 1] = toString(i)
  curr_batch_start_ind <- curr_batch_start_ind + num_cells
}
mat.annot <- data.frame(mat.annot)
rownames(mat.annot) <- colnames(mat)
colnames(mat.annot) = c("BATCH_ID")

#--------------------------------------------------COMPUTE FEATURES-------------------------------------------------#

# Create Seurat object
# min.features determined for each dataset
so <- CreateSeuratObject(counts = mat, project = name, min.features = 200, meta.data = mat.annot)

# Clean up
rm(mat)
rm(datalist)

# QC
so <- subset(so, 
             subset = nFeature_RNA > quantile(so$nFeature_RNA, 0.05) & 
               nFeature_RNA < quantile(so$nFeature_RNA, 0.95))

# Regress out batch effects
so.list <- SplitObject(so, split.by = "BATCH_ID")
for (i in 1:length(so.list)) {
  so.list[[i]] <- NormalizeData(so.list[[i]])
  so.list[[i]] <- FindVariableFeatures(so.list[[i]], nfeatures = vargenes)
}
so.anchors <- FindIntegrationAnchors(object.list = so.list, dims = 1:30, k.filter=100, anchor.features = vargenes)
so.integrated <- IntegrateData(anchorset = so.anchors, dims = 1:30, features.to.integrate = rownames(so))
# This is a bit of a hack, which makes it so that we don't have to change our pipeline
so@assays$RNA@data <- so.integrated@assays$integrated@data

# No need to normalize an integrated dataset
# so <- NormalizeData(so, normalization.method = "LogNormalize", scale.factor = 1000000)
so <- ScaleData(so, min.cells.to.block = 1, block.size = 500)

# Variables genes are pre-identified during integration; we don't recompute or plot anything
so@assays$RNA@var.features <- so.integrated@assays$integrated@var.features
# so <- FindVariableFeatures(so, nfeatures = vargenes)
# Plot variable genes with and without labels
# PlotAndSaveHVG(so, name)

# Run PCA
so <- RunPCA(so, npcs = 100)
# Project PCA to all genes
so <- ProjectDim(so, do.center = T)
# Plot Elbow
PlotAndSavePCAElbow(so, 100, name)

# Run ICA
so <- RunICA(so, nics = number_pcs)
# Project ICA to all genes
so <- ProjectDim(so, reduction = "ica", do.center = T)

# Cluster cells
so <- FindNeighbors(so, k.param = 15, dims = 1:number_pcs, nn.eps = 0)
so <- FindClusters(so, resolution = clus_res, n.start = 100)

# UMAP dim reduction
so <- RunUMAP(so, dims = 1:number_pcs, min.dist = 0.2, n.epochs = 500,
              n.neighbors = 15, learning.rate = 0.1, spread = 2)

# Check for batch effects
PlotAndSaveUMAPClusters(so, so@meta.data$BATCH_ID, name, suffix = "_batch_effects")
# Batch effects have been regressed out

# Plot UMAP clusters
PlotAndSaveUMAPClusters(so, so@meta.data$seurat_clusters, name)

# Plot PCs on UMAP
PlotAndSavePCsOnUMAP(so, name)
# Plot ICs on UMAP
PlotAndSaveICsOnUMAP(so, name)
# Plot known marker genes on UMAP 
# marker_genes <- c()
# PlotAndSaveKnownMarkerGenesOnUMAP(so, keep, marker_genes, name)

# Save global features
SaveGlobalFeatures(so, name)

# Compute any cluster dependent features (DE genes, within-cluster PCs, etc.) and save them
# Seurat clusters
Idents(object=so) <- "seurat_clusters"
clus <- levels(so@meta.data$seurat_clusters)
demarkers <- WithinClusterFeatures(so, "seurat_clusters", clus, name)

# Plot DE genes on UMAP
PlotAndSaveDEGenesOnUMAP(so, demarkers, name, height = 10, rank_by_tstat = TRUE)

# Save Seurat object
saveRDS(so, paste0("../data/", name, "/so.rds"))

