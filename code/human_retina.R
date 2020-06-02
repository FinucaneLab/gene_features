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
name <- "human_retina"
number_pcs <- 50
vargenes <- 2000
clus_res <- 0.2

# Setup
dir.create(paste0("../plots/", name))
dir.create(paste0("../features/", name))

# Notes on data:
# No annotations provided
# There appear to be two donor-types, so we first check for batch effects
# There are strong batch effects so run the integration pipeline

#------------------------------------------------LOAD AND FORMAT DATA-----------------------------------------------#

# Load this in in case we need it later
keep <- read.table("../resources/gene_annot_jun10.txt", sep = "\t", header = T, stringsAsFactors = F, col.names = c("ENSG", "symbol", "chr", "start", "end", "TSS"))

# Load and process first mat
file <- paste0("../data/", name, "/GSE116106_retina_aggr_10_matrix.mtx.gz")
mat1 <- readMM(file) %>%
  data.matrix() %>%
  Matrix(sparse = TRUE)
file <- paste0("../data/", name, "/GSE116106_retina_aggr_10_barcodes.tsv.gz")
mat1.colnames <- data.frame(fread(file, header=F), row.names=1)
file <- paste0("../data/", name, "/GSE116106_retina_aggr_10_genes.tsv.gz")
mat1.rownames <- data.frame(fread(file, header=F), row.names=1)
rownames(mat1) <- rownames(mat1.rownames)
colnames(mat1) <- rownames(mat1.colnames)
mat1 <- ConvertToENSGAndProcessMatrix(mat1, "ENSG")
mat1 <- mat1[match(keep$ENSG, rownames(mat1)),]

# Load and process second mat
file <- paste0("../data/", name, "/GSE122970_retina_aggr_day8_matrix.mtx.gz")
mat2 <- readMM(file) %>%
  data.matrix() %>%
  Matrix(sparse = TRUE)
file <- paste0("../data/", name, "/GSE122970_retina_aggr_day8_barcodes.tsv.gz")
mat2.colnames <- data.frame(fread(file, header=F), row.names=1)
file <- paste0("../data/", name, "/GSE122970_retina_aggr_day8_genes.tsv.gz")
mat2.rownames <- data.frame(fread(file, header=F), row.names=1)
rownames(mat2) <- rownames(mat2.rownames)
colnames(mat2) <- rownames(mat2.colnames)
mat2 <- ConvertToENSGAndProcessMatrix(mat2, "ENSG")
mat2 <- mat2[match(keep$ENSG, rownames(mat2)),]

# Bind assuming same row order
mat <- cbind(mat1, mat2)
colnames(mat) <- make.names(colnames(mat), unique=T)

# Make batch indicators
mat.annot <- matrix(0, nrow = dim(mat)[2], ncol = 1)
mat.annot[1:dim(mat1)[2], 1] = "1"
mat.annot[(1 + dim(mat1)[2]):(dim(mat)[2]), 1] = "2"
mat.annot <- data.frame(mat.annot)
rownames(mat.annot) <- colnames(mat)
colnames(mat.annot) = c("BATCH_ID")

#--------------------------------------------------COMPUTE FEATURES-------------------------------------------------#

# Create Seurat object
# min.features determined for each dataset
so <- CreateSeuratObject(counts = mat, project = name, min.features = 200, meta.data = mat.annot)

# Clean up
rm(mat)
rm(mat1)
rm(mat2)

# QC
so <- subset(so, 
             subset = nFeature_RNA > quantile(so$nFeature_RNA, 0.05) & 
               nFeature_RNA < quantile(so$nFeature_RNA, 0.95))

# Regress out batch effects
so.list <- SplitObject(so, split.by = "BATCH_ID")
for (i in 1:length(so.list)) {
  print(i)
  so.list[[i]] <- NormalizeData(so.list[[i]])
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
PlotAndSaveUMAPClusters(so, so@meta.data$BATCH_ID, name, suffix = "_batch_effects")
# This isn't perfect but is maybe alright, given that cell-type proportions in the gestation period versus 8 days postnatal could conceivably be very different

# Plot UMAP clusters
PlotAndSaveUMAPClusters(so, so@meta.data$seurat_clusters, name)

# Plot PCs on UMAP
PlotAndSavePCsOnUMAP(so, name)
# Plot ICs on UMAP
PlotAndSaveICsOnUMAP(so, name)
# Plot known marker genes on UMAP 
marker_genes <- c("EBF3", "POU4F2", "NHLH1", "DLX2", "NHLH2", "TFAP2B", "TFAP2A", "SCRT2", "PRDM13", "MEIS2", "POU2F2", "NEUROD4", "PTF1A", "THRB", "RAX2", "PRDM8", "RUNX1T1", "RAX2")
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

