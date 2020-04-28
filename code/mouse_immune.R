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
options(future.globals.maxSize = 2048 * 1024^2)
plan(strategy = "multicore", workers = 32)

# Parameters
name <- "mouse_immune"

# Setup
dir.create(paste0("../plots/", name))
dir.create(paste0("../features/", name))

# Notes on data:
# Data is already collapsed across donors
# Thus, we just report the normalized samples themselves, one v rest logFC, and a few PCs/ICs

#------------------------------------------------LOAD AND FORMAT DATA-----------------------------------------------#

# Read in data and annotations
file <- paste0("../data/", name, "/17aug_commonRNA_normalized.txt.gz")
mat <- data.frame(fread(file, skip = 0, sep = "\t")) %>% 
  data.matrix() %>% 
  Matrix(sparse = TRUE)
row.names(mat) <- read.table(paste0("../data/", name, "/rnaRowLabels.txt")) %>% .$V1

# Convert to ENSG, drop duplicates, and fill in missing genes
mat <- ConvertToENSGAndProcessMatrix(mat, "mouse_symbol")

# Load this in in case we need it later
keep <- read.table("../resources/gene_annot_jun10.txt", sep = "\t", header = T, stringsAsFactors = F, col.names = c("ENSG", "symbol", "chr", "start", "end", "TSS"))

#--------------------------------------------------COMPUTE FEATURES-------------------------------------------------#

so <- CreateSeuratObject(counts = mat, project = name)
so <- NormalizeData(so, normalization.method = "LogNormalize", scale.factor = 1000000)
so <- ScaleData(so)

# Run PCA (we use all features for PCA)
so <- RunPCA(so, npcs = 85, features = row.names(so))
# Project PCA to all genes
so <- ProjectDim(so, do.center = T)
# Plot Elbow
PlotAndSavePCAElbow(so, 85, name)

# Identify variable genes
so <- FindVariableFeatures(so, nfeatures = 4000)
# Plot variable genes with and without labels
PlotAndSaveHVG(so, name)

# Run ICA
so <- RunICA(so, nics = 30)
# Project ICA to all genes
so <- ProjectDim(so, reduction = "ica", do.center = T)

# Write out projected gene loadings across all cells
so@reductions$pca@feature.loadings.projected %>%
  data.frame() %>%
  rownames_to_column(., var = "ENSG") %>%
  as.tibble() %>%
  arrange(factor(ENSG, levels = keep$ENSG)) %>%
  fwrite(., 
         paste0("../features/", name, "/projected_pcaloadings.txt"),
         quote = F, row.names = F, col.names = T, sep = "\t")
system(paste0("gzip ../features/", name, "/projected_pcaloadings.txt"))

# Write out projected IC gene loadings across all cells
so@reductions$ica@feature.loadings.projected %>%
  data.frame() %>%
  rownames_to_column(., var = "ENSG") %>%
  as.tibble() %>%
  arrange(factor(ENSG, levels = keep$ENSG)) %>%
  fwrite(., 
         paste0("../features/", name, "/projected_icaloadings.txt"),
         quote = F, row.names = F, col.names = T, sep = "\t")
system(paste0("gzip ../features/", name, "/projected_icaloadings.txt"))

# Write out the samples (we name it average_expression for consistency with rest of features)
so@assays$RNA@scale.data %>%
  data.frame() %>%
  rownames_to_column(., var="ENSG") %>%
  as.tibble() %>%
  arrange(factor(ENSG, levels = keep$ENSG)) %>%
  fwrite(.,
         paste0("../features/", name, "/average_expression.txt"),
         quote = F, row.names = F, col.names = T, sep = "\t")
system(paste0("gzip ../features/", name, "/average_expression.txt"))

# Compute LogFC as just difference between log transformed scaled samples. We do this to avoid infinities.
lfc <- lapply(
  X = 1:86,
  FUN = function(x) {
    print(x)
    set1 <- so@assays$RNA@scale.data[,x] %>%
      data.matrix()
    set2 <- so@assays$RNA@scale.data[,c(1:86)[-x]] %>%
      data.matrix() %>%
      rowSums() %>%
      data.matrix()
    diff = set1 - set2
    return(diff)
  }
) %>%
  set_names(c(1:86)) %>%
  bind_rows()
row.names(lfc) <- row.names(so)

# Write out LFC (we name it as tstat from pre-defined clusters for consistency with the rest of the features)
lfc %>%
  rownames_to_column(., var="ENSG") %>%
  as.tibble() %>%
  arrange(factor(ENSG, levels = keep$ENSG)) %>%
  fwrite(.,
         paste0("../features/", name, "/diffexprs_tstat_clusters_pre_def.txt"),
         quote = F, row.names = F, col.names = T, sep = "\t")
system(paste0("gzip ../features/", name, "/diffexprs_tstat_clusters_pre_def.txt"))

# Save Seurat object
saveRDS(so, paste0("../data/", name, "/so.rds"))

