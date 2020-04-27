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
options(future.globals.maxSize = 2048 * 1024^2 * 2)
plan(strategy = "multicore", workers = 32)

# Parameters
name <- "human_multiple"
number_pcs <- 60
vargenes <- 5000
clus_res <- 0.6

# Setup
dir.create(paste0("../plots/", name))

# Notes on data:
# Bulk data in TPMs

#------------------------------------------------LOAD AND FORMAT DATA-----------------------------------------------#

# Read in data and annotations
file <- paste0("../data/", name, "/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz")
mat <- data.frame(fread(cmd = paste0("zcat < ", file), skip = 2, sep = "\t"), row.names = 1)[-1] %>%
  data.matrix() %>%
  Matrix(sparse = TRUE)
row.names(mat) <- gsub("\\..*", "", row.names(mat))
file <- paste0("../data/", name, "/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")
mat.annot <- data.frame(fread(file), skip=1) %>%
  dplyr::select(SAMPID, SMTSD) %>%
  dplyr::mutate(SAMPID = gsub("-", ".", SAMPID)) 
colnames(mat.annot) <- c("NAME", "Cluster")
row.names(mat.annot) <- make.names(mat.annot$NAME)

# Read in annotations
keep <- read.table("../resources/gene_annot_jun10.txt", sep = "\t", header = T, stringsAsFactors = F, col.names = c("ENSG", "symbol", "chr", "start", "end", "TSS"))

# Filter expression matrix
rowkeep <- row.names(mat) %in% keep$ENSG
mat <- mat[rowkeep,]

# Add missing genes as sparse rows
notkeep <- keep %>%
  dplyr::filter(ENSG %ni% row.names(mat))
missing <- sparseMatrix(dims = c(dim(notkeep)[1], length(colnames(mat))), i={}, j={})
row.names(missing) <- notkeep$ENSG
colnames(missing) <- colnames(mat)
mat <- rbind(mat, missing)

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
so@assays$RNA@data <- log1p(so@assays$RNA@counts)
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
so <- RunICA(so, nics = 60)
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
PlotAndSaveUMAPClusters(so, so@meta.data$Cluster, name, suffix = "_pre_def", width = 8, height = 7)

# Plot PCs on UMAP
PlotAndSavePCsOnUMAP(so, name)
# Plot ICs on UMAP
PlotAndSaveICsOnUMAP(so, name)
# Plot known marker genes on UMAP 
#marker_genes <- c("")
#PlotAndSaveKnownMarkerGenesOnUMAP(so, keep, marker_genes, name)

# Save global features
SaveGlobalFeatures(so, name)

# Compute any cluster dependent features (DE genes, within-cluster PCs, etc.) and save them
# Seurat clusters
Idents(object=so) <- "seurat_clusters"
clus <- levels(so@meta.data$seurat_clusters)
demarkers <- WithinClusterFeatures(so, "seurat_clusters", clus, name)

# Pre-defined cluster dependent features (if applicable)
Idents(object=so) <- "Cluster"
clus <- unique(so@meta.data$Cluster)
demarkers_pre_def <- WithinClusterFeatures(so, "Cluster", clus, name, suffix = "_pre_def")

# Plot DE genes on UMAP
PlotAndSaveDEGenesOnUMAP(so, demarkers, name, height = 30)
# Plot DE genes from pre-defined clusters on UMAP (if applicable)
PlotAndSaveDEGenesOnUMAP(so, demarkers_pre_def, name, suffix = "_pre_def", height = 30)

# Save Seurat object
saveRDS(so, paste0("../data/", name, "/so.rds"))




