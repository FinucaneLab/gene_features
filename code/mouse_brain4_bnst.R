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
options(future.globals.maxSize = 4 * 2048 * 1024^2)
plan(strategy = "multicore", workers = 32)

# Parameters
name <- "mouse_brain4_bnst"
number_pcs <- 70
vargenes <- 2000
clus_res <- 0.2

# Setup
dir.create(paste0("../plots/", name))
dir.create(paste0("../features/", name))

# Notes on data:
# No annotations provided.

#------------------------------------------------LOAD AND FORMAT DATA-----------------------------------------------#

# Load this in in case we need it later
keep <- read.table("../resources/gene_annot_jun10.txt", sep = "\t", header = T, stringsAsFactors = F, col.names = c("ENSG", "symbol", "chr", "start", "end", "TSS"))

# Read in data and annotations
# We'll have to read in each batch and cbind
read_sparse_mat <- function(x) {
  mat <- readMM(x) %>%
    Matrix(sparse = TRUE) %>%
    as("dgCMatrix")
  file <- sub("matrix", "genes", x)
  file <- sub("mtx", "csv", file)
  rows <- data.frame(fread(file, header=T), row.names=1)
  file <- sub("matrix", "barcodes", x)
  file <- sub("mtx", "csv", file)
  cols <- data.frame(fread(file, header=T), row.names=1)
  rownames(mat) <- rownames(rows)
  colnames(mat) <- rownames(cols)
  mat <- ConvertToENSGAndProcessMatrix(mat, "mouse_symbol")
  ### Reorder according to keep so that we can cbind later
  mat <- mat[match(keep$ENSG, rownames(mat)),]
  return(mat)
}

filelist <- Sys.glob(paste0("../data/", name, "/*.mtx.gz"))
# Read all matrices into a list
datalist <- lapply(filelist, read_sparse_mat)

# Bind assuming same row order
mat <- do.call("cbind", datalist)
colnames(mat) <- make.names(colnames(mat), unique=T)

#--------------------------------------------------COMPUTE FEATURES-------------------------------------------------#

# Create Seurat object
# min.features determined for each dataset
so <- CreateSeuratObject(counts = mat, project = name, min.features = 200)

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

# Plot PCs on UMAP
PlotAndSavePCsOnUMAP(so, name)
# Plot ICs on UMAP
PlotAndSaveICsOnUMAP(so, name)
# Plot known marker genes on UMAP 
marker_genes <- c("SLC17A6", "GAD2", "GAD1", "SLN", "CHAT", "IMPG1", "SEMA3E", "TRP73", "WIF1", "ARHGEF38", "ESR1", "GIL3", "LRRC9", "NPFFR2", "PAX6", "TH", "ZEB2", "ESR2", "ANGPT1", "BNC2", "EBF1", "EPSTI1", "HAUS4", "ROR1", "SST", "TAC2", "AVP", "CDH23", "PIEZO2", "VIP")
PlotAndSaveKnownMarkerGenesOnUMAP(so, keep, marker_genes, name)

# Save global features
SaveGlobalFeatures(so, name)

# Compute any cluster dependent features (DE genes, within-cluster PCs, etc.) and save them
# Seurat clusters
Idents(object=so) <- "seurat_clusters"
clus <- levels(so@meta.data$seurat_clusters)
demarkers <- WithinClusterFeatures(so, "seurat_clusters", clus, name)

# Plot DE genes on UMAP
PlotAndSaveDEGenesOnUMAP(so, demarkers, name, height = 30, rank_by_tstat = TRUE)

# Save Seurat object
saveRDS(so, paste0("../data/", name, "/so.rds"))

