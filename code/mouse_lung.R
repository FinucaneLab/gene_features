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
options(future.globals.maxSize = 1024 * 8 * 1024^2)
plan(strategy = "multicore", workers = 32)

# Parameters
name <- "mouse_lung"
number_pcs <- 50
vargenes <- 2000
clus_res <- 0.2

# Setup
dir.create(paste0("../plots/", name))
dir.create(paste0("../features/", name))

# Notes on data:
# Data is provided in batches so we regress out batch effects (we choose sequencing batch instead of amplification batch)
# Tissues annotations are provided but some are NaN, so we ignore them
# We see batch effects so we use Seurat's integration workflow

#------------------------------------------------LOAD AND FORMAT DATA-----------------------------------------------#

# Load this in in case we need it later
keep <- read.table("../resources/gene_annot_jun10.txt", sep = "\t", header = T, stringsAsFactors = F, col.names = c("ENSG", "symbol", "chr", "start", "end", "TSS"))

# Read in data and annotations
# We'll have to read in each batch and cbind
read_sparse_mat <- function(x) {
  mat <- data.frame(fread(paste0("../data/", name, "/", x)), row.names=1) %>%
    data.matrix() %>%
    Matrix(sparse = TRUE)
  mat <- ConvertToENSGAndProcessMatrix(mat, "mouse_symbol")
  ### Reorder according to keep so that we can cbind later
  mat <- mat[match(keep$ENSG, rownames(mat)),]
  return(mat)
}

# Get file names
filelist = list.files(path=paste0("../data/", name, "/"), pattern = "^GSM")
# Read all matrices into a list
datalist <- lapply(filelist, read_sparse_mat)

# Bind assuming same row order
mat <- do.call("cbind", datalist)

# Read in annotations
file <- paste0("../data/", name, "/GSE119228_metadata.txt.gz")
mat.annot <- data.frame(fread(file, skip=12), row.names=1)

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
so.list <- SplitObject(so, split.by = "Seq_batch_ID")
for (i in 1:length(so.list)) {
  print(i)
  so.list[[i]] <- NormalizeData(so.list[[i]])
  so.list[[i]] <- FindVariableFeatures(so.list[[i]], nfeatures = vargenes)
}
sample_size_list = lapply(so.list, FUN = function(x) dim(x)[2])
reference_dataset = which.max(sample_size_list)
so.anchors <- FindIntegrationAnchors(object.list = so.list, dims = 1:30, k.filter=100, anchor.features = vargenes, reference = reference_dataset)
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
so <- FindNeighbors(so, dims = 1:number_pcs, nn.eps = 0)
so <- FindClusters(so, resolution = clus_res, n.start = 100)

# UMAP dim reduction
so <- RunUMAP(so, dims = 1:number_pcs, min.dist = 0.4, n.epochs = 500,
              n.neighbors = 20, learning.rate = 0.1, spread = 2)

# Check for batch effects
PlotAndSaveUMAPClusters(so, so@meta.data$Seq_batch_ID, name, suffix = "_batch_effects")
# Batch effects appear mostly controlled for

# Plot UMAP clusters
PlotAndSaveUMAPClusters(so, so@meta.data$seurat_clusters, name)

# Plot PCs on UMAP
PlotAndSavePCsOnUMAP(so, name)
# Plot ICs on UMAP
PlotAndSaveICsOnUMAP(so, name)
# Plot known marker genes on UMAP 
marker_genes <- c("IL1RL1", "IL33", "MCPT8", "IL6", "IL13", "IL1B", "TNF", "CXCL2", "CSF2", "POU2F2")
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

