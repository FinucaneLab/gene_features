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
options(future.globals.maxSize = 12 * 2048 * 1024^2)
plan(strategy = "multicore", workers = 16)

# Parameters
name <- "human_fetalblood"
number_pcs <- 50
vargenes <- 2000
clus_res <- 0.2

# Setup
dir.create(paste0("../plots/", name))
dir.create(paste0("../features/", name))

# Notes on data:
# Annotations provided. There don't appear to be serious batch effects, and we capture biology regardless
# (looking at pre-defined annotations), so we don't do batch correction. We use the pre-defined cell type labels.

#------------------------------------------------LOAD AND FORMAT DATA-----------------------------------------------#

# Load this in in case we need it later
keep <- read.table("../resources/gene_annot_jun10.txt", sep = "\t", header = T, stringsAsFactors = F, col.names = c("ENSG", "symbol", "chr", "start", "end", "TSS"))

# Read in data and annotations
# We'll have to read in each batch and cbind
read_sparse_mat <- function(x) {
  file <- paste0("../data/", name, "/", x, "/GRCh38/matrix.mtx")
  mat <- readMM(file) %>%
    Matrix(sparse = TRUE) %>%
    as("dgCMatrix")
  file <- paste0("../data/", name, "/", x, "/GRCh38/genes.tsv")
  rows <- data.frame(fread(file, header=F), row.names=1)
  file <- paste0("../data/", name, "/", x, "/GRCh38/barcodes.tsv")
  cols <- data.frame(fread(file, header=F), row.names=1)
  rownames(mat) <- rownames(rows)
  colnames(mat) <- paste0(x, "_", lapply(rownames(cols), FUN = function(s) substr(s, 1, 16)))
  mat <- ConvertToENSGAndProcessMatrix(mat, "ENSG")
  ### Reorder according to keep so that we can cbind later
  mat <- mat[match(keep$ENSG, rownames(mat)),]
  return(mat)
}
read_annotations <- function(x) {
  file <- paste0("../data/", name, "/", x, "/", x, ".csv")
  mat.annot <- data.frame(fread(file, header=T), row.names=1)
  rownames(mat.annot) <- paste0(x, "_", rownames(mat.annot))
  return(mat.annot)
}

pathlist <- list.dirs(path = paste0("../data/", name, "/"), full.names = FALSE, recursive = FALSE)
# Read all matrices into a list
datalist <- lapply(pathlist, read_sparse_mat)
# Read all annotations into a list
ann_list <- lapply(pathlist, read_annotations)

# Bind assuming same row order
mat <- do.call("cbind", datalist)
mat.annot <- do.call("rbind", ann_list)

# Read in experiment annotations
file <- paste0("../data/", name, "/E-MTAB-7407.sdrf.txt")
exp_annot <- data.frame(fread(file, header=T), row.names=1)

# Make batch indicators
mat.annot$BATCH_ID <- NA
curr_batch_start_ind <- 1
for (i in 1:length(pathlist)) {
  batch_name <- exp_annot[pathlist[[i]],"Characteristics.age."]
  batch_name <- strsplit(batch_name," ")[[1]][1]
  if (batch_name %in% c("7", "8")) {
    batch_name <- "7-8"
  } else if (batch_name %in% c("9", "10", "11")) {
    batch_name <- "9-11"
  } else if (batch_name %in% c("12", "13", "14")) {
    batch_name <- "12-14"
  } else if (batch_name %in% c("15", "16", "17")) {
    batch_name <- "15-17"
  }
  num_cells <- dim(ann_list[[i]])[1]
  mat.annot[curr_batch_start_ind:(curr_batch_start_ind + num_cells - 1), "BATCH_ID"] = batch_name
  curr_batch_start_ind <- curr_batch_start_ind + num_cells
}

# Drop un-annotated cells
mat.annot <- mat.annot[mat.annot$Cell.Labels != "NA",]
mat <- mat[,colnames(mat) %in% rownames(mat.annot)]

#--------------------------------------------------COMPUTE FEATURES-------------------------------------------------#

# Subset columns because data is giant
sample_col_inds <- sample(dim(mat)[2], 40000)
mat <- mat[,sample_col_inds]
mat.annot <- mat.annot[sample_col_inds,,drop=F] 

# Create Seurat object
# min.features determined for each dataset
so <- CreateSeuratObject(counts = mat, project = name, min.features = 200, meta.data = mat.annot)

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
PlotAndSaveUMAPClusters(so, so@meta.data$BATCH_ID, name, suffix = "_batch_effects")

# Plot UMAP clusters
PlotAndSaveUMAPClusters(so, so@meta.data$seurat_clusters, name)
# Plot known clusters on UMAP (if applicable)
PlotAndSaveUMAPClusters(so, so@meta.data$Cell.Labels, name, suffix = "_pre_def")

# Plot PCs on UMAP
PlotAndSavePCsOnUMAP(so, name)
# Plot ICs on UMAP
PlotAndSaveICsOnUMAP(so, name)
# Plot known marker genes on UMAP 
marker_genes <- c("CD34", "SPINK2", "JCHAIN", "IGLL1", "CD79B", "TCL1A", "IGKC", "MS4A1", "CD19", "LTB", "KLRB1", "CD7", "IL32", "NKG7", "XCL2", "NCAM1", "MPO", "LYZ", "CD1C", "CD4", "C1QA", "VCAM1", "KIT", "PF4", "ITGA2B", "KLF1", "ALAS2", "HBA1", "BPGM", "ESAM", "ECM1", "APOA1")
PlotAndSaveKnownMarkerGenesOnUMAP(so, keep, marker_genes, name)

# Save global features
SaveGlobalFeatures(so, name)

# Compute any cluster dependent features (DE genes, within-cluster PCs, etc.) and save them
# Seurat clusters
Idents(object=so) <- "seurat_clusters"
clus <- levels(so@meta.data$seurat_clusters)
demarkers <- WithinClusterFeatures(so, "seurat_clusters", clus, name)
# Pre-defined cluster dependent features (if applicable)
Idents(object=so) <- "Cell.Labels"
clus <- unique(so@meta.data$Cell.Labels)
demarkers_pre_def <- WithinClusterFeatures(so, "Cell.Labels", clus, name, suffix = "_pre_def")

# Plot DE genes on UMAP
PlotAndSaveDEGenesOnUMAP(so, demarkers, name, height = 30, rank_by_tstat = TRUE)
# Plot DE genes from pre-defined clusters on UMAP (if applicable)
PlotAndSaveDEGenesOnUMAP(so, demarkers_pre_def, name, suffix = "_pre_def", height = 30, rank_by_tstat = TRUE)

# Save Seurat object
saveRDS(so, paste0("../data/", name, "/so.rds"))

