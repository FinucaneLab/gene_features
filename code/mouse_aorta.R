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
source("utils.R")

# Set up parallelization
# Remember to use htop to delete forgotten forks
Sys.setenv(R_FUTURE_FORK_ENABLE = T)
options(future.globals.maxSize = 2048 * 1024^2)
plan(strategy = "multicore", workers = 32)

# Parameters
name <- "mouse_aorta"
number_pcs <- 40
vargenes <- 2000
clus_res <- 0.6

# Setup
dir.create(paste0("../plots/", name))

# Notes on data:
# Authors provide cluster labels and what appears to be the first 2 PCs
# We'll just use the cluster labels

#------------------------------------------------LOAD AND FORMAT DATA-----------------------------------------------#

# Read in data and annotations
file <- paste0("../data/", name, "/Seurat_Chow_12PCs_outfile.mtx")
mat <- data.frame(fread(file), row.names=1) %>%
  data.matrix() %>%
  Matrix(sparse = TRUE)
file <- paste0("../data/", name, "/META_DATA_Chow_12PCs_outfile.txt")
mat.annot <- data.frame(fread(cmd = paste0("zcat < ", file), skip=2))
colnames(mat.annot) <- c("NAME", "Cluster", "Sub.Cluster", "Average.Intensity")
row.names(mat.annot) <- make.names(mat.annot$NAME)

# Convert to ENSG
sym2ensmusg <- read.table("../resources/symbol2ensmusg.txt", header = F, stringsAsFactors = F, col.names = c("symbol", "ENSMUSG"))
ensmusg2ensg <- read.table("../resources/ensmusg2ensg.txt", header = F, stringsAsFactors = F, col.names = c("ENSMUSG", "ENSG"))
### Rename and drop duplicates
ensmusg2ensg <- ensmusg2ensg[!duplicated(ensmusg2ensg$ENSMUSG),]
row.names(sym2ensmusg) <- sym2ensmusg$symbol
row.names(ensmusg2ensg) <- ensmusg2ensg$ENSMUSG
### First Sym -> ENSG
mat <- merge(mat, sym2ensmusg, by="row.names")
row.names(mat) <- mat$ENSMUSG
mat <- mat[, !(names(mat) %in% c("Row.names", "symbol", "ENSMUSG"))]
### Now ENSG -> ENSMUSG
mat <- merge(mat, ensmusg2ensg, by="row.names")
### Drop duplicates again
mat <- mat[!duplicated(mat$ENSG),]
row.names(mat) <- mat$ENSG
mat <- mat[, !(names(mat) %in% c("Row.names", "ENSMUSG", "ENSG"))]

# Read in annotations
keep <- read.table("../resources/gene_annot_jun10.txt", sep = "\t", header = T, stringsAsFactors = F, col.names = c("ENSG", "symbol", "chr", "start", "end", "TSS"))

# Filter expression matrix
rowkeep <- row.names(mat) %in% keep$ENSG
mat <- mat[rowkeep,]

# Add missing genes as sparse rows
notkeep <- keep %>%
  dplyr::filter(ENSG %ni% row.names(mat))
missing <- sparseMatrix(dims = c(dim(notkeep)[1], length(colnames(mat))), i={}, j={}) %>%
  data.frame()
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
PlotAndSaveUMAPClusters(so, so@meta.data$Sub.Cluster, name, suffix = "_pre_def")

# Plot PCs on UMAP
PlotAndSavePCsOnUMAP(so, name)
# Plot ICs on UMAP
PlotAndSaveICsOnUMAP(so, name)
# Plot known marker genes on UMAP 
marker_genes <- c("MYH11", "TPM2", "MYL9", "TAGLN", "ACTA2", "C1QA", "PECAM1", "GPIHBP1", "CD52", "RAC2")
PlotAndSaveKnownMarkerGenesOnUMAP(so, keep, marker_genes, name)

# Save global features
SaveGlobalFeatures(so, name)

# Compute any cluster dependent features (DE genes, within-cluster PCs, etc.) and save them
# Seurat clusters
Idents(object=so) <- "seurat_clusters"
clus <- levels(so@meta.data$seurat_clusters)
demarkers <- WithinClusterFeatures(so, clus, name)
# Pre-defined cluster dependent features (if applicable)
Idents(object=so) <- "Sub.Cluster"
clus <- unique(so@meta.data$Sub.Cluster)
demarkers_pre_def <- WithinClusterFeatures(so, clus, name, suffix = "_pre_def")

# Plot DE genes on UMAP
PlotAndSaveDEGenesOnUMAP(so, demarkers, name)
# Plot DE genes from pre-defined clusters on UMAP (if applicable)
PlotAndSaveDEGenesOnUMAP(so, demarkers_pre_def, name, suffix = "_pre_def")

# Save Seurat object
saveRDS(so, paste0("../data/", name, "/so.rds"))

