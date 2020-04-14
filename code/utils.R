# Not %in% function
"%ni%" <- Negate("%in%")

# t-test
DiffTTest <- function (object, cells.1, cells.2, genes.use = NULL, print.bar = TRUE, 
          assay.type = "RNA") 
{
  data.test <- GetAssayData(object = object, assay.type = assay.type, 
                            slot = "data")
  genes.use <- SetIfNull(x = genes.use, default = rownames(x = data.test))
  if (print.bar) {
    iterate.fxn = pblapply
  }
  else {
    iterate.fxn <- lapply
  }
  p_val <- unlist(x = iterate.fxn(X = genes.use, FUN = function(x) {
    as.numeric(t.test(x = data.test[x, cells.1], y = data.test[x, cells.2])$statistic)
  }))
  to.return <- data.frame(p_val, row.names = genes.use)
  return(to.return)
}

# T stat
my.t.test <- function(c){
  n <- sqrt(length(c))
  mean(c)*n/sd(c)
}

# Plot and save variable genes
PlotAndSaveHVG <- function(so, name, display = T) {
  p <- so@assays$RNA@meta.features %>%
    dplyr::rename("variable_gene" = vst.variable) %>%
    ggplot(.) +
    geom_histogram(aes(x = log10(vst.variance.standardized), fill = variable_gene), bins = 50) +
    pretty_plot() + 
    scale_fill_manual(values = jdb_palette("FantasticFox")[c(3,5)]) +
    xlab("log10 mean-standardized variance")
  if (display) {
    plot(p)
  }
  ggsave(p, filename = paste0("../plots/", name, "/variablegenes.pdf"), device = cairo_pdf, width = 6, height = 4, family = "Helvetica")
}

# Plot and save PCA Elbow
PlotAndSavePCAElbow <- function(so, num_pcs, name, display = T) {
  p <- cbind("PC" = seq(1:num_pcs), "stdev" = so@reductions$pca@stdev) %>%
    as.tibble() %>%
    ggplot(., aes(x =  PC, y = stdev)) +
    geom_point() + 
    pretty_plot() +
    ylab("standard deviation") + 
    xlab("principal components")
  if(display) {
    plot(p)
  }
  ggsave(p, filename = paste0("../plots/", name, "/pcaelbow.pdf"), device = cairo_pdf, width = 6, height = 4, family = "Helvetica")
}

# Plot and save UMAP Clusters
PlotAndSaveUMAPClusters <- function(so, clust_col, name, suffix = "", display = T) {
  clusters.df <- bind_cols(data.frame("cluster" = clust_col),
                           data.frame(so@reductions$umap@cell.embeddings)) %>%
    as.tibble()
  p <- clusters.df %>%
    ggplot(., aes(x = UMAP_1, y = UMAP_2)) + 
    geom_point_rast(aes(color = cluster), size = 0.5) +
    scale_color_manual(values = jdb_palette("corona")) +
    pretty_plot()
  if(display) {
    plot(p)
  }
  ggsave(p + theme(legend.position = "right"), filename = paste0("../plots/", name, "/umap_clusters", suffix, ".pdf"), device = cairo_pdf, width = 6, height = 5, family = "Helvetica")
}

# Plot and save PCs on UMAP
# No option to display because figure is usually large
PlotAndSavePCsOnUMAP <- function(so, name) {
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
}

# Plot and save ICs on UMAP
# No option to display because figure is usually large
PlotAndSaveICsOnUMAP <- function(so, name) {
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
}

# Plot and save known marker genes on UMAP
# No option to display because figure is usually large
PlotAndSaveKnownMarkerGenesOnUMAP <- function(so, keep, marker_genes, name) {
  knownmarker_genes <- keep %>%
    dplyr::select(ENSG, symbol) %>%
    dplyr::filter(symbol %in% marker_genes)
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
}

# Save global features
SaveGlobalFeatures <- function(so, name) {
  # Write out projected gene loadings across all cells
  so@reductions$pca@feature.loadings.projected %>%
    data.frame() %>%
    rownames_to_column(., var = "ENSG") %>%
    as.tibble() %>%
    arrange(factor(ENSG, levels = keep$ENSG)) %>%
    fwrite(., 
           paste0("../features/", name, "/projected_pcaloadings.txt"),
           quote = F, row.names = F, col.names = T, sep = "\t")
  
  # Write out ICA across all cells
  so@reductions$ica@feature.loadings.projected %>%
    data.frame() %>%
    rownames_to_column(., var = "ENSG") %>%
    as.tibble() %>%
    arrange(factor(ENSG, levels = keep$ENSG)) %>%
    fwrite(., 
           paste0("../features/", name, "/projected_icaloadings.txt"),
           quote = F, row.names = F, col.names = T, sep = "\t")
}

# Big within-cluster features function
WithinClusterFeatures <- function(so, clus, name, suffix="") {
  # Calculate differentially expressed genes
  markers <- FindAllMarkers(so, test.use = "t", logfc.threshold = 0)
  markers <- markers %>%
    merge(., keep  %>% dplyr::select(ENSG, symbol), by.x = "gene", by.y = "ENSG") %>%
    dplyr::mutate(p_val = ifelse(p_val < 10^-200, 10^-200, p_val),
                  tstat = qt(p_val * 2, dim(so)[2] - 2, lower.tail = F) * sign(avg_logFC)) %>%
    as.tibble()
  #Get upregulated
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
  #Get downregulated
  demarkers_down <- markers %>%
    dplyr::filter(p_val_adj < 0.05, avg_logFC < -log(2))
  demarkers_down %>%
    group_by(cluster) %>%
    count()
  demarkers_down.mat <- demarkers_down %>%
    dplyr::mutate(value = 1) %>%
    dplyr::select(cluster, gene, value) %>%
    dplyr::mutate(cluster = paste0("Cluster", cluster)) %>%
    cast_sparse(gene, cluster, value)
  notkeep <- keep %>%
    dplyr::filter(ENSG %ni% row.names(demarkers_down.mat))
  missing <- sparseMatrix(dims = c(dim(notkeep)[1], length(colnames(demarkers_down.mat))), i={}, j={})
  row.names(missing) <- notkeep$ENSG
  colnames(missing) <- colnames(demarkers_down.mat)
  demarkers_down.df <- rbind(demarkers_down.mat, missing) %>%
    data.frame()
  #Format continuous
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
  so.clus.pcs <- lapply(1:length(clus), function(x) {PC_within_cluster(so, clus[x])})
  so.clus.pcs[sapply(so.clus.pcs, is.null)] <- NULL
  so.clus.pcs <- so.clus.pcs %>%
    do.call(cbind, .)
  
  # Write out projected gene loadings within clusters
  so.clus.pcs %>%
    data.frame() %>%
    rownames_to_column(., var = "ENSG") %>%
    as.tibble() %>%
    arrange(factor(ENSG, levels = keep$ENSG)) %>%
    fwrite(., 
           paste0("../features/", name, "/projected_pcaloadings_clusters", suffix, ".txt"),
           quote = F, row.names = F, col.names = T, sep = "\t")
  
  # Write out normalized expression within clusters and across all cells
  so.ae <- AverageExpression(so, slot = "scale.data")$RNA
  colnames(so.ae) <- paste0("Cluster", colnames(so.ae))
  so.ae$Allcells <- apply(so@assays$RNA@scale.data, 1, mean)
  so.ae %>%
    data.frame() %>%
    rownames_to_column(., var = "ENSG") %>%
    as.tibble() %>%
    arrange(factor(ENSG, levels = keep$ENSG)) %>%
    fwrite(., 
           paste0("../features/", name, "/average_expression", suffix, ".txt"),
           quote = F, row.names = F, col.names = T, sep = "\t")
  
  # Write differential expression (DE genes) between clusters
  demarkers.df %>%
    rownames_to_column(., var = "ENSG") %>%
    as.tibble() %>%
    arrange(factor(ENSG, levels = keep$ENSG)) %>%
    fwrite(.,
           paste0("../features/", name, "/diffexprs_genes_clusters", suffix, ".txt"),
           quote = F, row.names = F, col.names = T, sep = "\t")

  # Write differential expression (DE genes) between clusters (downregulated)
  demarkers_down.df %>%
    rownames_to_column(., var = "ENSG") %>%
    as.tibble() %>%
    arrange(factor(ENSG, levels = keep$ENSG)) %>%
    fwrite(.,
           paste0("../features/", name, "/diffexprs_down_genes_clusters", suffix, ".txt"),
           quote = F, row.names = F, col.names = T, sep = "\t")

  # Write differential expression (t-stat) between clusters
  markers.df %>%
    rownames_to_column(., var = "ENSG") %>%
    as.tibble() %>%
    arrange(factor(ENSG, levels = keep$ENSG)) %>%
    fwrite(.,
           paste0("../features/", name, "/diffexprs_tstat_clusters", suffix, ".txt"),
           quote = F, row.names = F, col.names = T, sep = "\t")

  return(demarkers)
}

# Plot top DE genes on UMAP
# No option to display because figure is usually large
PlotAndSaveDEGenesOnUMAP <- function(so, demarkers, name, suffix = "") {
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
  ggsave(p + theme(legend.position = "none"), filename = paste0("../plots/", name, "/umap_degenes", suffix, ".pdf"), device = cairo_pdf, width = 7, height = 10, family = "Helvetica")
}

