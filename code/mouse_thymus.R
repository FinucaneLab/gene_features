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
name <- "mouse_thymus"
number_pcs <- 40
vargenes <- 2500
clus_res <- 0.8

# Setup
dir.create(paste0("../plots/", name))
dir.create(paste0("../features/", name))

# Notes on data:
# There appear to be mild batch effects, but trying to correct for them seems to make a lot of structure disappear
# It also appears that the datasets might come from different sources; i.e. whole thymus vs thymic lobe?
# We don't try to regress out batch effects, hoping that they are mild enough and/or correspond to actual biology
# No annotations are provided

#------------------------------------------------LOAD AND FORMAT DATA-----------------------------------------------#

# Load this in in case we need it later
keep <- read.table("../resources/gene_annot_jun10.txt", sep = "\t", header = T, stringsAsFactors = F, col.names = c("ENSG", "symbol", "chr", "start", "end", "TSS"))

# Read in data and annotations
# We'll have to read in each batch and cbind
read_sparse_mat <- function(x) {
  mat <- data.frame(fread(x), row.names=1) %>%
    data.matrix() %>%
    Matrix(sparse = TRUE)
  mat <- ConvertToENSGAndProcessMatrix(mat, "mouse_symbol")
  ### Reorder according to keep so that we can cbind later
  mat <- mat[match(keep$ENSG, rownames(mat)),]
  return(mat)
}

filelist = Sys.glob(paste0("../data/", name, "/*"))
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

# Check for batch effects
PlotAndSaveUMAPClusters(so, so@meta.data$BATCH_ID, name, suffix = "_batch_effects")
# Is ok

# Plot UMAP clusters
PlotAndSaveUMAPClusters(so, so@meta.data$seurat_clusters, name)

# Plot PCs on UMAP
PlotAndSavePCsOnUMAP(so, name)
# Plot ICs on UMAP
PlotAndSaveICsOnUMAP(so, name)
# Plot known marker genes on UMAP 
marker_genes <- c("DLK1", "COL3A1", "VIM", "COL4A1", "IGF2", "SPARC", "KRT8", "PTN", "MDK", "SFN", "KRT5", "CRIP1", "LCP1", "CD8A", "TCF7", "CTSS")
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

