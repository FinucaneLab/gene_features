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
source("utils.R")

# Parameters
name <- "human_pancreas"
number_pcs <- 40
vargenes <- 3000
clus_res <- 1

# Setup
dir.create(paste0("../plots/", name))

# Set up parallelization
# Remember to use htop to delete forgotten forks
Sys.setenv(R_FUTURE_FORK_ENABLE = T)
options(future.globals.maxSize = 2048 * 1024^2)
plan(strategy = "multicore", workers = 32)

# Notes on data:
# Authors DID NOT provide cell type / cluster IDs for cells
# Authors provided gene symbols, some are duplicated

# Read in data and annotations
file <- paste0("../data/", name, "/GSE85241_cellsystems_dataset_4donors_updated.csv.gz")
mat <- data.frame(fread(cmd = paste0("zcat < ", file), skip = 1, sep = "\t"))[,-1] %>%
  data.matrix() %>% 
  Matrix()
rownames(mat) <- data.frame(cmd = fread(paste0("zcat < ", file), select = 1, skip = 1, sep = "\t"))[,1] %>%
  gsub("__.*", "", .)
file <- paste0("../data/", name, "/cellnames.txt")
mat.annot <-  read_delim(file, delim = " ", col_names = c("cell", "donor", "library", "id")) %>% 
  data.frame()
row.names(mat.annot) <- mat.annot$cell
colnames(mat) <- mat.annot$cell

# Read in annotations
keep <- read.table("../resources/gene_annot_jun10.txt", sep = "\t", header = T, stringsAsFactors = F, col.names = c("ENSG", "symbol", "chr", "start", "end", "TSS"))

# Filter expression matrix
rowkeep <- row.names(mat) %in% keep$symbol
mat <- mat[rowkeep,]

# Deal with duplicate gene symbols
mat <- mat[!duplicated(row.names(mat)),]
dups <- keep %>%
  dplyr::filter(symbol %in% row.names(mat)) %>%
  count(symbol) %>%
  dplyr::filter(n > 1) %>%
  .$symbol
keep <- keep %>%
  group_by(symbol) %>%
  dplyr::mutate(size = abs(start - end),
                rank = rank(size))
matdups <- mat[dups,]
keepdups <- keep %>%
  dplyr::filter(rank == 2)
row.names(matdups) <- keepdups[match(row.names(matdups), keepdups$symbol),]$ENSG
keepnotdups <- keep %>%
  dplyr::filter(rank == 1)
row.names(mat) <- keepnotdups[match(row.names(mat), keepnotdups$symbol),]$ENSG
mat <- rbind(mat, matdups)

# Add missing genes as sparse rows
notkeep <- keep %>%
  dplyr::filter(ENSG %ni% row.names(mat))
missing <- sparseMatrix(dims = c(dim(notkeep)[1], length(colnames(mat))), i={}, j={})
row.names(missing) <- notkeep$ENSG
colnames(missing) <- colnames(mat)
mat <- rbind(mat, missing)

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
p <- so@assays$RNA@meta.features %>%
  dplyr::rename("variable_gene" = vst.variable) %>%
  ggplot(.) +
  geom_histogram(aes(x = log10(vst.variance.standardized), fill = variable_gene), bins = 50) +
  pretty_plot() + 
  scale_fill_manual(values = jdb_palette("FantasticFox")[c(3,5)]) +
  xlab("log10 mean-standardized variance")
p
ggsave(p, filename = paste0("../plots/", name, "/variablegenes.pdf"), device = cairo_pdf, width = 6, height = 4, family = "Helvetica")

# Run PCA
so <- RunPCA(so, npcs = 100)

# Project PCA to all genes
so <- ProjectDim(so, do.center = T)

# Plot Elbow
p <- cbind("PC" = seq(1:100), "stdev" = so@reductions$pca@stdev) %>%
  as.tibble() %>%
  ggplot(., aes(x =  PC, y = stdev)) +
  geom_point() + 
  pretty_plot() +
  ylab("standard deviation") + 
  xlab("principal components")
p
ggsave(p, filename = paste0("../plots/", name, "/pcaelbow.pdf"), device = cairo_pdf, width = 6, height = 4, family = "Helvetica")

# Run ICA
so <- RunICA(so, nics = number_pcs)

# Project ICA to all genes
so <- ProjectDim(so, reduction = "ica", do.center = T)

# Cluster cells
# higher k.param?
# plot with PCA to verify cluster number? some quantitative metric?
so <- FindNeighbors(so, dims = 1:number_pcs, nn.eps = 0)
so <- FindClusters(so, resolution = clus_res, n.start = 100)

# UMAP dim reduction
so <- RunUMAP(so, dims = 1:number_pcs, min.dist = 0.4, n.epochs = 500,
                                     n.neighbors = 10, learning.rate = 0.1, spread = 2)

# Plot UMAP clusters
clusters.df <- bind_cols(data.frame("cluster" = so@meta.data$seurat_clusters),
                        data.frame(so@reductions$umap@cell.embeddings)) %>%
  as.tibble()
p <- clusters.df %>%
  ggplot(., aes(x = UMAP_1, y = UMAP_2)) + 
  geom_point_rast(aes(color = cluster), size = 0.5) +
  scale_color_manual(values = jdb_palette("corona")) +
  pretty_plot() 
ggsave(p + theme(legend.position = "none"), filename = paste0("../plots/", name, "/umap_clusters.pdf"), device = cairo_pdf, width = 6, height = 5, family = "Helvetica")

# Plot PCs on UMAP
pcs.df <- bind_cols(data.frame(so@reductions$pca@cell.embeddings[,1:24]),
                        data.frame(so@reductions$umap@cell.embeddings)) %>%
  as.tibble()
p <- pcs.df %>%
  reshape2::melt(id.vars = c("UMAP_1", "UMAP_2")) %>%
  dplyr::mutate(value = case_when(value > 3 ~ 3,
                                  value < -3 ~ -3,
                                  T ~ value)) %>%
  dplyr::rename("scores" = value) %>%
  ggplot(., aes(x = UMAP_1, y = UMAP_2)) + 
  geom_point_rast(aes(color = scores), size = 1) +
  scale_color_gradientn(colors = jdb_palette("solar_extra")) +
  pretty_plot() +
  facet_wrap(~variable, ncol = 4)
ggsave(p + theme(legend.position = "none"), filename = paste0("../plots/", name, "/umap_pcs.pdf"), device = cairo_pdf, width = 7, height = 10, family = "Helvetica")

# Plot ICs on UMAP
ics.df <- bind_cols(data.frame(so@reductions$ica@cell.embeddings[,1:24]),
                    data.frame(so@reductions$umap@cell.embeddings)) %>%
  as.tibble()
p <- ics.df %>%
  reshape2::melt(id.vars = c("UMAP_1", "UMAP_2")) %>%
  dplyr::mutate(value = case_when(value > 3 ~ 3,
                                  value < -3 ~ -3,
                                  T ~ value)) %>%
  dplyr::rename("scores" = value) %>%
  ggplot(., aes(x = UMAP_1, y = UMAP_2)) + 
  geom_point_rast(aes(color = scores), size = 1) +
  scale_color_gradientn(colors = jdb_palette("solar_extra")) +
  pretty_plot() +
  facet_wrap(~variable, ncol = 4)
ggsave(p + theme(legend.position = "none"), filename = paste0("../plots/", name, "/umap_ics.pdf"), device = cairo_pdf, width = 7, height = 10, family = "Helvetica")

# Plot known marker genes on UMAP  
knownmarker_genes <- keep %>%
  dplyr::select(ENSG, symbol) %>%
  dplyr::filter(symbol %in% c("GCG", "INS", "SST", "MAFA", "PPY", "PRSS1", "KRT19", "MAFB", "GSTA2", "GHRL", "ESAM"))
knownmarkers.df <- bind_cols(data.frame(t(so@assays$RNA@scale.data[knownmarker_genes$ENSG,])),
                        data.frame(so@reductions$umap@cell.embeddings)) %>%
  as.tibble()
p <- knownmarkers.df %>%
  reshape2::melt(id.vars = c("UMAP_1", "UMAP_2")) %>%
  merge(., knownmarker_genes, by.x = "variable", by.y = "ENSG") %>%
  as.tibble() %>%
  dplyr::mutate(value = case_when(value > 3 ~ 3,
                          value < -3 ~ -3,
                          T ~ value)) %>%
  dplyr::rename("exprs" = value) %>%
  ggplot(., aes(x = UMAP_1, y = UMAP_2)) + 
  geom_point_rast(aes(color = exprs), size = 1) +
  scale_color_gradientn(colors = jdb_palette("solar_extra")) +
  pretty_plot() +
  facet_wrap(~symbol, ncol = 4)
ggsave(p + theme(legend.position = "none"), filename = paste0("../plots/", name, "/umap_knownmarkers.pdf"), device = cairo_pdf, width = 7, height = 10, family = "Helvetica")

# Calculate differentially expressed genes
# is t test in log space?
# include down regulated as separate feature?
markers <- FindAllMarkers(so, test.use = "t", logfc.threshold = 0)
markers <- markers %>%
  merge(., keep  %>% dplyr::select(ENSG, symbol), by.x = "gene", by.y = "ENSG") %>%
  dplyr::mutate(p_val = ifelse(p_val < 10^-200, 10^-200, p_val),
                tstat = qt(p_val * 2, dim(so)[2] - 2, lower.tail = F) * sign(avg_logFC)) %>%
  as.tibble()
demarkers <- markers %>% 
  dplyr::filter(p_val_adj < 0.05, avg_logFC > log(2)) 
demarkers %>%
  group_by(cluster) %>% 
  count()
demarkers.mat <- demarkers %>%
  dplyr::mutate(value = 1) %>%
  dplyr::select(cluster, gene, value) %>%
  dplyr::mutate(cluster = paste0("Cluster", cluster)) %>%
  cast_sparse(gene, cluster, value)
notkeep <- keep %>%
  dplyr::filter(ENSG %ni% row.names(demarkers.mat))
missing <- sparseMatrix(dims = c(dim(notkeep)[1], length(colnames(demarkers.mat))), i={}, j={})
row.names(missing) <- notkeep$ENSG
colnames(missing) <- colnames(demarkers.mat)
demarkers.df <- rbind(demarkers.mat, missing) %>%
  data.frame()
markers.mat <- markers %>%
  dplyr::select(cluster, gene, tstat) %>%
  dplyr::mutate(cluster = paste0("Cluster", cluster)) %>%
  cast_sparse(gene, cluster, tstat)
notkeep <- keep %>%
  dplyr::filter(ENSG %ni% row.names(markers.mat))
missing <- sparseMatrix(dims = c(dim(notkeep)[1], length(colnames(markers.mat))), i={}, j={})
row.names(missing) <- notkeep$ENSG
colnames(missing) <- colnames(markers.mat)
markers.df <- rbind(markers.mat, missing) %>%
  data.frame()

# Plot DE genes on UMAP
topdegenes <- demarkers %>%
  group_by(cluster) %>% 
  top_n(n = 2, wt = avg_logFC)
topdegenes.df <- bind_cols(data.frame(t(so@assays$RNA@scale.data[topdegenes$gene,])),
                             data.frame(so@reductions$umap@cell.embeddings)) %>%
  as.tibble()
p <- topdegenes.df %>%
  reshape2::melt(id.vars = c("UMAP_1", "UMAP_2")) %>%
  merge(., keep, by.x = "variable", by.y = "ENSG") %>%
  as.tibble() %>%
  dplyr::mutate(value = case_when(value > 3 ~ 3,
                                  value < -3 ~ -3,
                                  T ~ value)) %>%
  dplyr::rename("exprs" = value) %>%
  ggplot(., aes(x = UMAP_1, y = UMAP_2)) + 
  geom_point_rast(aes(color = exprs), size = 1) +
  scale_color_gradientn(colors = jdb_palette("solar_extra")) +
  pretty_plot() +
  facet_wrap(~symbol, ncol = 4)
ggsave(p + theme(legend.position = "none"), filename = paste0("../plots/", name, "/umap_degenes.pdf"), device = cairo_pdf, width = 7, height = 10, family = "Helvetica")

# Calculate PCs within cluster
PC_within_cluster <- function(so, clus) {
  so.clus <- subset(so, idents = clus) 
  if(dim(so.clus)[2] > 100) {
    so.clus <- RunPCA(so.clus, npcs = 10)
    so.clus <- ProjectDim(so.clus, do.center = T)
    so.gl <- so.clus@reductions$pca@feature.loadings.projected
    colnames(so.gl) <- paste0("Cluster", clus, "_", colnames(so.gl))
    return(so.gl)
  }
}
clus <- levels(so@meta.data$seurat_clusters)
so.clus.pcs <- lapply(1:length(clus), function(x) {PC_within_cluster(so, clus[x])})
so.clus.pcs[sapply(so.clus.pcs, is.null)] <- NULL
so.clus.pcs <- so.clus.pcs %>%
  do.call(cbind, .)

# Save Seurat object
saveRDS(so, paste0("../data/", name, "/so.rds"))

# Write out projected gene loadings across all cells
so@reductions$pca@feature.loadings.projected %>%
  data.frame() %>%
  rownames_to_column(., var = "ENSG") %>%
  as.tibble() %>%
  arrange(factor(ENSG, levels = keep$ENSG)) %>%
fwrite(., 
       "../features/human_pancreas/projected_pcaloadings.txt", 
       quote = F, row.names = F, col.names = T, sep = "\t")

# Write out projected gene loadings within clusters
so.clus.pcs %>%
  data.frame() %>%
  rownames_to_column(., var = "ENSG") %>%
  as.tibble() %>%
  arrange(factor(ENSG, levels = keep$ENSG)) %>%
  fwrite(., 
         "../features/human_pancreas/projected_pcaloadings_clusters.txt", 
         quote = F, row.names = F, col.names = T, sep = "\t")

# Write out projected gene loadings within pre-defined clusters (where available)

# Write out normalized expression within clusters and across all cells
human_pancreas.ae <- AverageExpression(so, slot = "scale.data")$RNA
colnames(human_pancreas.ae) <- paste0("Cluster", colnames(human_pancreas.ae))
human_pancreas.ae$Allcells <- apply(so@assays$RNA@scale.data, 1, mean)
human_pancreas.ae %>%
  data.frame() %>%
  rownames_to_column(., var = "ENSG") %>%
  as.tibble() %>%
  arrange(factor(ENSG, levels = keep$ENSG)) %>%
  fwrite(., 
         "../features/human_pancreas/average_expression.txt", 
         quote = F, row.names = F, col.names = T, sep = "\t")

# Write out normalized expression within pre-defined clusters (where available)

# Write differential expression (DE genes) between clusters
demarkers.df %>%
  rownames_to_column(., var = "ENSG") %>%
  as.tibble() %>%
  arrange(factor(ENSG, levels = keep$ENSG)) %>%
  fwrite(., 
         "../features/human_pancreas/diffexprs_genes_clusters.txt", 
         quote = F, row.names = F, col.names = T, sep = "\t")

# Write differential expression (t-stat) between clusters
markers.df %>%
  rownames_to_column(., var = "ENSG") %>%
  as.tibble() %>%
  arrange(factor(ENSG, levels = keep$ENSG)) %>%
  fwrite(., 
         "../features/human_pancreas/diffexprs_tstat_clusters.txt", 
         quote = F, row.names = F, col.names = T, sep = "\t")

# Write out ICA across all cells
so@reductions$ica@feature.loadings.projected %>%
  data.frame() %>%
  rownames_to_column(., var = "ENSG") %>%
  as.tibble() %>%
  arrange(factor(ENSG, levels = keep$ENSG)) %>%
  fwrite(., 
         "../features/human_pancreas/projected_icaloadings.txt", 
         quote = F, row.names = F, col.names = T, sep = "\t")

