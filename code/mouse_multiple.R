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
name <- "mouse_multiple"
number_pcs <- 60
vargenes <- 2500
clus_res <- 0.6

# Setup
dir.create(paste0("../plots/", name))
dir.create(paste0("../features/", name))

# Notes on data:
# 

#------------------------------------------------LOAD AND FORMAT DATA-----------------------------------------------#

# Read in data and annotations
file = paste0("../data/", name, "/Figure2-batch-removed.txt")
mat <- data.frame(fread(file), row.names=1) %>%
  data.matrix() %>%
  Matrix(sparse = TRUE)
file <- paste0("../data/", name, "/MCA_CellAssignments.csv")
mat.annot <- data.frame(fread(file, header=T), row.names=2)

# Convert to ENSG, drop duplicates, and fill in missing genes
mat <- ConvertToENSGAndProcessMatrix(mat, "mouse_symbol")

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
PlotAndSaveUMAPClusters(so, so@meta.data$Tissue, name, suffix = "_pre_def")

# Plot PCs on UMAP
PlotAndSavePCsOnUMAP(so, name)
# Plot ICs on UMAP
PlotAndSaveICsOnUMAP(so, name)

# Save global features
SaveGlobalFeatures(so, name)

# Compute any cluster dependent features (DE genes, within-cluster PCs, etc.) and save them
# Seurat clusters
Idents(object=so) <- "seurat_clusters"
clus <- levels(so@meta.data$seurat_clusters)
demarkers <- WithinClusterFeatures(so, "seurat_clusters", clus, name)
# Pre-defined cluster dependent features (if applicable)
Idents(object=so) <- "Tissue"
clus <- unique(so@meta.data$Tissue)
clus <- clus[!is.na(clus)]
so@meta.data$Tissue[is.na(so@meta.data$Tissue)] = "NA"
demarkers_pre_def <- WithinClusterFeatures(so, "Tissue", clus, name, suffix = "_pre_def")

# Plot DE genes on UMAP
PlotAndSaveDEGenesOnUMAP(so, demarkers, name, height = 30, rank_by_tstat = TRUE)
# Plot DE genes from pre-defined clusters on UMAP (if applicable)
PlotAndSaveDEGenesOnUMAP(so, demarkers_pre_def, name, suffix = "_pre_def", height = 30, rank_by_tstat = TRUE)

# Write out normalized expression within sub-clusters and across all cells
Idents(object=so) <- "ClusterID"
so.ae <- AverageExpression(so, slot = "scale.data")$RNA
colnames(so.ae) <- paste0("Cluster", colnames(so.ae))
so.ae$Allcells <- apply(so@assays$RNA@scale.data, 1, mean)
so.ae %>%
  data.frame() %>%
  rownames_to_column(., var = "ENSG") %>%
  as.tibble() %>%
  arrange(factor(ENSG, levels = keep$ENSG)) %>%
  fwrite(., 
         paste0("../features/", name, "/average_expression", "_pre_def_2", ".txt"),
         quote = F, row.names = F, col.names = T, sep = "\t")
system(paste0("gzip ../features/", name, "/average_expression", "_pre_def_2", ".txt"))

# Save Seurat object
saveRDS(so, paste0("../data/", name, "/so.rds"))

