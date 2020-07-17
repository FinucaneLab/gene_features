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
options(future.globals.maxSize = 8 * 2048 * 1024^2)
plan(strategy = "multicore", workers = 32)

# Parameters
name <- "mouse_multiple2_facs"
number_pcs <- 65
vargenes <- 1500
clus_res <- 0.2

# Setup
dir.create(paste0("../plots/", name))
dir.create(paste0("../features/", name))

# Notes on data:
# Batch labels are provided, though I don't see any batch effects on the UMAP, so I ignore them.
# Cell-type annotations provided, but there are a ton of these,
# so we skip plotting differentially expressed genes

#------------------------------------------------LOAD AND FORMAT DATA-----------------------------------------------#

# Load this in in case we need it later
keep <- read.table("../resources/gene_annot_jun10.txt", sep = "\t", header = T, stringsAsFactors = F, col.names = c("ENSG", "symbol", "chr", "start", "end", "TSS"))

file <- paste0("../data/", name, "/facs_filtered_matrix.mtx")
mat <- readMM(file) %>%
  Matrix(sparse = TRUE) %>%
  as("dgCMatrix")
file <- paste0("../data/", name, "/facs_filtered_genes.tsv")
rows <- data.frame(fread(file, header=T), row.names=1)
file <- paste0("../data/", name, "/facs_filtered_cells.tsv")
mat.annot <- data.frame(fread(file, header=T), row.names=1)
rownames(mat) <- rownames(rows)
colnames(mat) <- rownames(mat.annot)
mat <- ConvertToENSGAndProcessMatrix(mat, "mouse_symbol")

#--------------------------------------------------COMPUTE FEATURES-------------------------------------------------#

# Subset columns because data is giant
sample_col_inds <- sample(dim(mat)[2], 60000)
mat <- mat[,sample_col_inds]
mat.annot <- mat.annot[sample_col_inds,,drop=F]

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
PlotAndSaveUMAPClusters(so, so@meta.data$cell_ontology_class, name, suffix = "_pre_def")

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
Idents(object=so) <- "cell_ontology_class"
clus <- unique(so@meta.data$cell_ontology_class)
demarkers_pre_def <- WithinClusterFeatures(so, "cell_ontology_class", clus, name, suffix = "_pre_def")

# Plot DE genes on UMAP
# PlotAndSaveDEGenesOnUMAP(so, demarkers, name, height = 20, rank_by_tstat = TRUE)
# Plot DE genes from pre-defined clusters on UMAP (if applicable)
# PlotAndSaveDEGenesOnUMAP(so, demarkers_pre_def, name, suffix = "_pre_def", height = 20, rank_by_tstat = TRUE)

# Save Seurat object
saveRDS(so, paste0("../data/", name, "/so.rds"))

