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
options(future.globals.maxSize = 3 * 4096 * 1024^2)
plan(strategy = "multicore", workers = 32)

# Parameters
name <- "human_thymus"
number_pcs <- 40
vargenes <- 2000
clus_res <- 0.4

# Setup
dir.create(paste0("../plots/", name))
dir.create(paste0("../features/", name))

# Notes on data:
# h5ad file was written using ScanPy/AnnData in Python, and was not compatible with R version of Seurat
# Read in the file in Python and wrote it out into matrix market format
# There are batch effects (across different donors) so we run the integration pipeline.

#------------------------------------------------LOAD AND FORMAT DATA-----------------------------------------------#

# Read in data and annotations
file <- paste0("../data/", name, "/human_thymus_data_processed.mtx")
mat <- readMM(file) %>%
  data.matrix() %>%
  Matrix(sparse = TRUE) %>%
  t()
file <- paste0("../data/", name, "/human_thymus_data_processed_genes.tsv")
rows <- data.frame(fread(file, header=F))
file <- paste0("../data/", name, "/human_thymus_data_processed_cells.tsv")
cols <- data.frame(fread(file, header=F), row.names=1)
rownames(rows) <- make.names(rows[,1], unique=T)
rownames(mat) <- rownames(rows)
colnames(mat) <- rownames(cols)

# Convert to ENSG
mat <- ConvertToENSGAndProcessMatrix(mat, "human_symbol")

# Read in annotations
file <- paste0("../data/", name, "/HTA08.v01.A06.Science_human_tcells.csv")
mat.annot <- data.frame(fread(file), row.names=1)

# Load this in in case we need it later
keep <- read.table("../resources/gene_annot_jun10.txt", sep = "\t", header = T, stringsAsFactors = F, col.names = c("ENSG", "symbol", "chr", "start", "end", "TSS"))

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
so.list <- SplitObject(so, split.by = "donor")
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
PlotAndSaveUMAPClusters(so, so@meta.data$donor, name, suffix = "_batch_effects")
# Batch effects appear mostly controlled for

# Plot UMAP clusters
PlotAndSaveUMAPClusters(so, so@meta.data$seurat_clusters, name)
# Plot known clusters on UMAP (if applicable)
PlotAndSaveUMAPClusters(so, so@meta.data$cell.types, name, suffix = "_pre_def")

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
# Pre-defined cluster dependent features (if applicable)
Idents(object=so) <- "cell.types"
clus <- unique(so@meta.data$cell.types)
demarkers_pre_def <- WithinClusterFeatures(so, "cell.types", clus, name, suffix = "_pre_def")

# Plot DE genes on UMAP
PlotAndSaveDEGenesOnUMAP(so, demarkers, name, height = 30, rank_by_tstat = TRUE)
# Plot DE genes from pre-defined clusters on UMAP (if applicable)
PlotAndSaveDEGenesOnUMAP(so, demarkers_pre_def, name, suffix = "_pre_def", height = 30, rank_by_tstat = TRUE)

# Save Seurat object
saveRDS(so, paste0("../data/", name, "/so.rds"))

