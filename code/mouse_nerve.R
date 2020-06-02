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
options(future.globals.maxSize = 2048 * 1024^2)
plan(strategy = "multicore", workers = 32)

# Parameters
name <- "mouse_nerve"
number_pcs <- 40
vargenes <- 1800
clus_res <- 1.0

# Setup
dir.create(paste0("../plots/", name))
dir.create(paste0("../features/", name))

# Notes on data:
# Cluster labels are provided, so we use these.

#------------------------------------------------LOAD AND FORMAT DATA-----------------------------------------------#

# Read in data and annotations
file <- paste0("../data/", name, "/GSE142541_mouse_matrix.mtx.gz")
mat <- readMM(file) %>%
  data.matrix() %>%
  Matrix(sparse = TRUE)
file <- paste0("../data/", name, "/GSE142541_mouse_barcodes.tsv.gz")
col.annot <- data.frame(fread(file, header=F), row.names=1)
file <- paste0("../data/", name, "/GSE142541_mouse_features.tsv.gz")
row.annot <- data.frame(fread(file, header=F), row.names=1)
colnames(mat) <- rownames(col.annot)
rownames(mat) <- rownames(row.annot)
file <- paste0("../data/", name, "/GSE142541_mouse_meta.csv.gz")
mat.annot <- data.frame(fread(file, header=T))
mat.annot <- mat.annot[!duplicated(mat.annot$V1), ]
rownames(mat.annot) <- mat.annot[,1]
mat <- mat[,colnames(mat) %in% rownames(mat.annot)]

# Convert to ENSG, drop duplicates, and fill in missing genes
mat <- ConvertToENSGAndProcessMatrix(mat, "ENSMUSG")

# Load this in in case we need it later
keep <- read.table("../resources/gene_annot_jun10.txt", sep = "\t", header = T, stringsAsFactors = F, col.names = c("ENSG", "symbol", "chr", "start", "end", "TSS"))

#--------------------------------------------------COMPUTE FEATURES-------------------------------------------------#

# Create Seurat object
# min.features determined for each dataset
so <- CreateSeuratObject(counts = mat, project = name, min.features = 200, meta.data = mat.annot)

# Clean up
rm(mat)

# QC
so <- subset(so, 
             subset = nFeature_RNA > quantile(so$nFeature_RNA, 0.05) & 
               nFeature_RNA < quantile(so$nFeature_RNA, 0.95))
so <- NormalizeData(so, normalization.method = "LogNormalize", scale.factor = 1000000)
so <- ScaleData(so, min.cells.to.block = 1, block.size = 500)

# Identify variable genes
so <- FindVariableFeatures(so, nfeatures = vargenes)
# Plot variable genes with and without labels
PlotAndSaveHVG(so, name)

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
so <- FindNeighbors(so, dims = 1:number_pcs, nn.eps = 0)
so <- FindClusters(so, resolution = clus_res, n.start = 100)

# UMAP dim reduction
so <- RunUMAP(so, dims = 1:number_pcs, min.dist = 0.4, n.epochs = 500,
              n.neighbors = 10, learning.rate = 0.1, spread = 2)

# Plot UMAP clusters
PlotAndSaveUMAPClusters(so, so@meta.data$seurat_clusters, name)
# Plot known clusters on UMAP (if applicable)
PlotAndSaveUMAPClusters(so, so@meta.data$x, name, suffix = "_pre_def")

# Plot PCs on UMAP
PlotAndSavePCsOnUMAP(so, name)
# Plot ICs on UMAP
PlotAndSaveICsOnUMAP(so, name)
# Plot known marker genes on UMAP 
marker_genes <- c("CD3E", "CD68", "LYZ2", "PTPRC", "SLC16A1", "CLDN1", "PROX1", "LYVE1", "RGS5", "ACTA2", "SFRP4", "ECM1", "DPT", "SMOC2", "APOD", "PLP1", "MBP", "SOX10", "NGFR")
PlotAndSaveKnownMarkerGenesOnUMAP(so, keep, marker_genes, name)

# Save global features
SaveGlobalFeatures(so, name)

# Compute any cluster dependent features (DE genes, within-cluster PCs, etc.) and save them
# Seurat clusters
Idents(object=so) <- "seurat_clusters"
clus <- levels(so@meta.data$seurat_clusters)
demarkers <- WithinClusterFeatures(so, "seurat_clusters", clus, name)
# Pre-defined cluster dependent features (if applicable)
Idents(object=so) <- "x"
clus <- unique(so@meta.data$x)
demarkers_pre_def <- WithinClusterFeatures(so, "x", clus, name, suffix = "_pre_def")

# Plot DE genes on UMAP
PlotAndSaveDEGenesOnUMAP(so, demarkers, name, height = 10, rank_by_tstat = TRUE)
# Plot DE genes from pre-defined clusters on UMAP (if applicable)
PlotAndSaveDEGenesOnUMAP(so, demarkers_pre_def, name, suffix = "_pre_def", height = 10, rank_by_tstat = TRUE)

# Save Seurat object
saveRDS(so, paste0("../data/", name, "/so.rds"))

