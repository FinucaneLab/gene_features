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
name <- "mouse_gastrulation"
number_pcs <- 50
vargenes <- 2000
clus_res <- 0.2

# Setup
dir.create(paste0("../plots/", name))
dir.create(paste0("../features/", name))

# Notes on data:
# I plotted the batch labels and saw no batch effects. Trying to integrate batches breaks clustering.
# So I ignore them. I use the cell type annotations. From these, it appears we are capturing the biology without batch correction.

#------------------------------------------------LOAD AND FORMAT DATA-----------------------------------------------#

# Load this in in case we need it later
keep <- read.table("../resources/gene_annot_jun10.txt", sep = "\t", header = T, stringsAsFactors = F, col.names = c("ENSG", "symbol", "chr", "start", "end", "TSS"))

file <- paste0("../data/", name, "/atlas/raw_counts.mtx")
mat <- readMM(file) %>%
  Matrix(sparse = TRUE) %>%
  as("dgCMatrix")
file <- paste0("../data/", name, "/atlas/genes.tsv")
rows <- data.frame(fread(file, header=F), row.names=1)
file <- paste0("../data/", name, "/atlas/barcodes.tsv")
cols <- data.frame(fread(file, header=F), row.names=1)
rownames(mat) <- rownames(rows)
colnames(mat) <- rownames(cols)
file <- paste0("../data/", name, "/atlas/meta.tab")
mat.annot <- data.frame(fread(file, header=T), row.names=1)
mat.annot$sequencing.batch <- as(mat.annot$sequencing.batch, "character")

mat <- ConvertToENSGAndProcessMatrix(mat, "ENSMUSG")

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

# Check for batch effects
PlotAndSaveUMAPClusters(so, so@meta.data$sequencing.batch, name, suffix = "_batch_effects")

# Plot UMAP clusters
PlotAndSaveUMAPClusters(so, so@meta.data$seurat_clusters, name)
# Plot known clusters on UMAP (if applicable)
PlotAndSaveUMAPClusters(so, so@meta.data$celltype, name, suffix = "_pre_def")

# Plot PCs on UMAP
PlotAndSavePCsOnUMAP(so, name)
# Plot ICs on UMAP
PlotAndSaveICsOnUMAP(so, name)
# Plot known marker genes on UMAP 
marker_genes <- c("ISL1", "GFPT2", "TNNT1", "HAS2", "FOXE1", "HHEX", "PHLDA2", "HHEX", "TTR", "SFRP5", "UPP1", "FZD7", "OSR1", "RIPPLY3", "HOXA1", "IRX1", "NEPN", "TFPI", "TRF", "XLR3A", "RHOX5", "TRAP1A", "CITED1", "CDX1", "HOXB2", "STRA6", "MNX1", "TLX2", "HOXC9", "AXIN2", "WNT5B", "SMIM3", "CDX2")
PlotAndSaveKnownMarkerGenesOnUMAP(so, keep, marker_genes, name)

# Save global features
SaveGlobalFeatures(so, name)

# Compute any cluster dependent features (DE genes, within-cluster PCs, etc.) and save them
# Seurat clusters
Idents(object=so) <- "seurat_clusters"
clus <- levels(so@meta.data$seurat_clusters)
demarkers <- WithinClusterFeatures(so, "seurat_clusters", clus, name)
# Pre-defined cluster dependent features (if applicable)
Idents(object=so) <- "celltype"
clus <- unique(so@meta.data$celltype)
clus <- clus[!is.na(clus)]
demarkers_pre_def <- WithinClusterFeatures(so, "celltype", clus, name, suffix = "_pre_def")

# Plot DE genes on UMAP
PlotAndSaveDEGenesOnUMAP(so, demarkers, name, height = 30, rank_by_tstat = TRUE)
# Plot DE genes from pre-defined clusters on UMAP (if applicable)
PlotAndSaveDEGenesOnUMAP(so, demarkers_pre_def, name, suffix = "_pre_def", height = 30, rank_by_tstat = TRUE)

# Save Seurat object
saveRDS(so, paste0("../data/", name, "/so.rds"))

