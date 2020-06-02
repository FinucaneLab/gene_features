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
name <- "human_coloncancer"
number_pcs <- 40
vargenes <- 2000
clus_res <- 0.2

# Setup
dir.create(paste0("../plots/", name))
dir.create(paste0("../features/", name))

# Notes on data:
# Data comes from two platforms, 10X and Smart-seq. We use the Seurat integration pipeline.
# Also, data appears to come in log(1 + TPM) form. This is because the files are named TPM but the columns don't add up to a million.
# Sub-cluster annotations provided, so we use those

#------------------------------------------------LOAD AND FORMAT DATA-----------------------------------------------#

# Load this in in case we need it later
keep <- read.table("../resources/gene_annot_jun10.txt", sep = "\t", header = T, stringsAsFactors = F, col.names = c("ENSG", "symbol", "chr", "start", "end", "TSS"))

# Read in data and annotations
# We'll have to read in each batch and cbind
read_sparse_mat <- function(x) {
  mat <- data.frame(fread(paste0("../data/", name, "/", x)), row.names=1) %>%
    data.matrix() %>%
    Matrix(sparse = TRUE)
  mat <- ConvertToENSGAndProcessMatrix(mat, "human_symbol")
  ### Reorder according to keep so that we can cbind later
  mat <- mat[match(keep$ENSG, rownames(mat)),]
  return(mat)
}

# Get file names
filelist = c("GSE146771_CRC.Leukocyte.10x.TPM.txt.gz", "GSE146771_CRC.Leukocyte.Smart-seq2.TPM.txt.gz")
# Read all matrices into a list
datalist <- lapply(filelist, read_sparse_mat)

# Bind assuming same row order
mat <- do.call("cbind", datalist)

# Read in annotations
file <- paste0("../data/", name, "/GSE146771_CRC.Leukocyte.10x.Metadata.txt.gz")
mat.annot.10x <- data.frame(fread(file), row.names=1)
file <- paste0("../data/", name, "/GSE146771_CRC.Leukocyte.Smart-seq2.Metadata.txt.gz")
mat.annot.smartseq <- data.frame(fread(file), row.names=1)
mat.annot <- rbind(mat.annot.10x, mat.annot.smartseq)
rownames(mat.annot) <- mat.annot[,"CellName"]

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
so.list <- SplitObject(so, split.by = "Platform")
for (i in 1:length(so.list)) {
  print(i)
  so.list[[i]] <- FindVariableFeatures(so.list[[i]], nfeatures = vargenes)
}
so.anchors <- FindIntegrationAnchors(object.list = so.list, dims = 1:30, k.filter=100, anchor.features = vargenes)
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
PlotAndSaveUMAPClusters(so, so@meta.data$Platform, name, suffix = "_batch_effects")
# Batch effects appear mostly controlled for

# Plot UMAP clusters
PlotAndSaveUMAPClusters(so, so@meta.data$seurat_clusters, name)
# Plot known clusters on UMAP (if applicable)
PlotAndSaveUMAPClusters(so, so@meta.data$Sub_Cluster, name, suffix = "_pre_def")

# Plot PCs on UMAP
PlotAndSavePCsOnUMAP(so, name)
# Plot ICs on UMAP
PlotAndSaveICsOnUMAP(so, name)
# Plot known marker genes on UMAP 
marker_genes <- c("CPA3", "TPSD1", "SLC18A2", "TCF4", "MZB1", "IL3RA", "FCER1A", "FCGR2B", "CD1A", "DAPP1", "C1ORF54", "IDO1", "VCAN", "CD36", "CD14", "IFITM2", "FCGR3A", "SIGLEC10", "EREG", "IL1B", "CCL3", "C1QB", "C1QC", "TREM2", "SPP1", "ABL2", "SDC4")
PlotAndSaveKnownMarkerGenesOnUMAP(so, keep, marker_genes, name)

# Save global features
SaveGlobalFeatures(so, name)

# Compute any cluster dependent features (DE genes, within-cluster PCs, etc.) and save them
# Seurat clusters
Idents(object=so) <- "seurat_clusters"
clus <- levels(so@meta.data$seurat_clusters)
demarkers <- WithinClusterFeatures(so, "seurat_clusters", clus, name)
# Pre-defined cluster dependent features (if applicable)
Idents(object=so) <- "Sub_Cluster"
clus <- unique(so@meta.data$Sub_Cluster)
demarkers_pre_def <- WithinClusterFeatures(so, "Sub_Cluster", clus, name, suffix = "_pre_def")

# Plot DE genes on UMAP
PlotAndSaveDEGenesOnUMAP(so, demarkers, name, height = 10, rank_by_tstat = TRUE)
# Plot DE genes from pre-defined clusters on UMAP (if applicable)
PlotAndSaveDEGenesOnUMAP(so, demarkers_pre_def, name, suffix = "_pre_def", height = 10, rank_by_tstat = TRUE)

# Save Seurat object
saveRDS(so, paste0("../data/", name, "/so.rds"))

