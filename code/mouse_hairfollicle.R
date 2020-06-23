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
name <- "mouse_hairfollicle"
number_pcs <- 60
vargenes <- 2500
clus_res <- 0.6

# Setup
dir.create(paste0("../plots/", name))
dir.create(paste0("../features/", name))

# Notes on data:
# No annotations provided

#------------------------------------------------LOAD AND FORMAT DATA-----------------------------------------------#

# Load this in in case we need it later
keep <- read.table("../resources/gene_annot_jun10.txt", sep = "\t", header = T, stringsAsFactors = F, col.names = c("ENSG", "symbol", "chr", "start", "end", "TSS"))

# Read in data and annotations
# We'll have to read in each batch and cbind
read_sparse_mat <- function(x) {
  mat <- readMM(paste0(x, "matrix.mtx.gz")) %>%
    data.matrix() %>%
    Matrix(sparse = TRUE)
  rows <- data.frame(fread(paste0(x, "genes.tsv.gz"), header=F), row.names=1)
  cols <- data.frame(fread(paste0(x, "barcodes.tsv.gz"), header=F), row.names=1)
  rownames(mat) <- rownames(rows)
  colnames(mat) <- rownames(cols)
  mat <- ConvertToENSGAndProcessMatrix(mat, "ENSMUSG")
  ### Reorder according to keep so that we can cbind later
  mat <- mat[match(keep$ENSG, rownames(mat)),]
  return(mat)
}
read_sparse_mat_2 <- function(x) {
  mat <- readMM(paste0(x, "matrix.mtx.gz")) %>%
    data.matrix() %>%
    Matrix(sparse = TRUE)
  rows <- data.frame(fread(paste0(x, "features.tsv.gz"), header=F), row.names=1)
  cols <- data.frame(fread(paste0(x, "barcodes.tsv.gz"), header=F), row.names=1)
  rownames(mat) <- rownames(rows)
  colnames(mat) <- rownames(cols)
  mat <- ConvertToENSGAndProcessMatrix(mat, "ENSMUSG")
  ### Reorder according to keep so that we can cbind later
  mat <- mat[match(keep$ENSG, rownames(mat)),]
  return(mat)
}

filelist = c(paste0("../data/", name, "/GSM3177991_sample_1_"),
             paste0("../data/", name, "/GSM3177992_sample_2_"),
             paste0("../data/", name, "/GSM4155928_sample_3_Wisoo_Young_Rev_Count_"),
             paste0("../data/", name, "/GSM4155929_sample_4_Wisoo_Aged_Rev_Count_"))
# Read all matrices into a list
datalist <- lapply(filelist[1:2], read_sparse_mat)
datalist_2 <- lapply(filelist[3:4], read_sparse_mat_2)

# Bind assuming same row order
mat <- do.call("cbind", datalist)
mat2 <- do.call("cbind", datalist_2)
mat <- do.call("cbind", c(mat, mat2))
colnames(mat) <- make.names(colnames(mat), unique=T)

#--------------------------------------------------COMPUTE FEATURES-------------------------------------------------#

# Create Seurat object
# min.features determined for each dataset
so <- CreateSeuratObject(counts = mat, project = name, min.features = 200)

# Clean up
rm(mat)
rm(datalist)

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
              n.neighbors = 20, learning.rate = 0.1, spread = 2)

# Plot UMAP clusters
PlotAndSaveUMAPClusters(so, so@meta.data$seurat_clusters, name)

# Plot PCs on UMAP
PlotAndSavePCsOnUMAP(so, name)
# Plot ICs on UMAP
PlotAndSaveICsOnUMAP(so, name)
# Plot known marker genes on UMAP 
marker_genes <- c("Mylk", "Des", "Tpm2", "Myh11", "Cnn1", "Actg2", "Ppp1r14a", "Ckb", "Myl9", "Acta2", "Flna", "Csrp1", "Msrb1", "Tagln", "Synpo2", "Ckm", "Cd59a")
marker_genes <- lapply(marker_genes, toupper)
PlotAndSaveKnownMarkerGenesOnUMAP(so, keep, marker_genes, name)

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




