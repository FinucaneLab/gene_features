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

### Assumes a sparse dgCMatrix as input
### Accepts row_id_type = ENSG, ENSMUSG, human_symbol, mouse_symbol
ConvertToENSGAndProcessMatrix <- function(mat, row_id_type) {
  ### Convert to ENSG
  if (row_id_type == "mouse_symbol") {
    sym2ensmusg <- read.table("../resources/symbol2ensmusg.txt", header = F, stringsAsFactors = F, col.names = c("symbol", "ENSMUSG"))
    ### Drop rows that don't exist in column
    mat <- mat[row.names(mat) %in% sym2ensmusg$symbol,]
    ### Rename
    idx <- match(row.names(mat), sym2ensmusg$symbol)
    row.names(mat) <- sym2ensmusg$ENSMUSG[idx]
    ### Update row type variable
    row_id_type <- "ENSMUSG"
  }
  if (row_id_type == "human_symbol") {
    sym2ensg <- read.table("../resources/ensg2symbol.txt", header = F, stringsAsFactors = F, col.names = c("ENSG", "symbol"))
    ### Drop rows that don't exist in column
    mat <- mat[row.names(mat) %in% sym2ensg$symbol,]
    ### Rename
    idx <- match(row.names(mat), sym2ensg$symbol)
    row.names(mat) <- sym2ensg$ENSG[idx]
    ### Update row type variable
    row_id_type <- "ENSG"
  }
  if (row_id_type == "ENSMUSG") {
    ensmusg2ensg <- read.table("../resources/ensmusg2ensg.txt", header = F, stringsAsFactors = F, col.names = c("ENSMUSG", "ENSG"))
    ### Drop rows that don't exist in column
    mat <- mat[row.names(mat) %in% ensmusg2ensg$ENSMUSG,]
    ### Rename
    idx <- match(row.names(mat), ensmusg2ensg$ENSMUSG)
    row.names(mat) <- ensmusg2ensg$ENSG[idx]
    ### Update row type variable
    row_id_type <- "ENSG"
  }
  ### Drop duplicate rows
  mat <- mat[!duplicated(row.names(mat)),]
  ### Read in annotations
  keep <- read.table("../resources/gene_annot_jun10.txt", sep = "\t", header = T, stringsAsFactors = F, col.names = c("ENSG", "symbol", "chr", "start", "end", "TSS"))
  ### Filter to genes of interest
  rowkeep <- row.names(mat) %in% keep$ENSG
  mat <- mat[rowkeep,]
  ### Add missing genes
  notkeep <- keep %>%
    dplyr::filter(ENSG %ni% row.names(mat))
  missing <- sparseMatrix(dims = c(dim(notkeep)[1], length(colnames(mat))), i={}, j={})
  row.names(missing) <- notkeep$ENSG
  colnames(missing) <- colnames(mat)
  mat <- rbind(mat, missing)
  ### Return
  return(mat)
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
PlotAndSaveUMAPClusters <- function(so, clust_col, name, suffix = "", display = T, raster_dpi = 300, width = 6, height = 5) {
  clusters.df <- bind_cols(data.frame("cluster" = clust_col),
                           data.frame(so@reductions$umap@cell.embeddings)) %>%
    as.tibble()
  p <- clusters.df %>%
    ggplot(., aes(x = UMAP_1, y = UMAP_2)) + 
    geom_point_rast(aes(color = cluster), size = 0.5, raster.dpi = raster_dpi) +
    scale_color_manual(values = c(jdb_palette("corona"), jdb_palette("corona"))) +
    pretty_plot() +
    guides(col = guide_legend(ncol = 1)) +
    theme(legend.key.size = unit(0.1, 'lines'))
  if(display) {
    plot(p)
  }
  ggsave(p + theme(legend.position = "right"), filename = paste0("../plots/", name, "/umap_clusters", suffix, ".pdf"), device = cairo_pdf, width = width, height = height, family = "Helvetica")
}

# Plot and save PCs on UMAP
PlotAndSavePCsOnUMAP <- function(so, name, display=T, raster_dpi=100) {
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
    geom_point_rast(aes(color = scores), size = 1, raster.dpi = raster_dpi) +
    scale_color_gradientn(colors = jdb_palette("solar_extra")) +
    pretty_plot() +
    facet_wrap(~variable, ncol = 4)
  if (display) {
    plot(p)
  }
  ggsave(p + theme(legend.position = "none"), filename = paste0("../plots/", name, "/umap_pcs.pdf"), device = cairo_pdf, width = 7, height = 10, family = "Helvetica")
}

# Plot and save ICs on UMAP
PlotAndSaveICsOnUMAP <- function(so, name, display=T, raster_dpi=100) {
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
    geom_point_rast(aes(color = scores), size = 1, raster.dpi = raster_dpi) +
    scale_color_gradientn(colors = jdb_palette("solar_extra")) +
    pretty_plot() +
    facet_wrap(~variable, ncol = 4)
  if (display) {
    plot(p)
  }
  ggsave(p + theme(legend.position = "none"), filename = paste0("../plots/", name, "/umap_ics.pdf"), device = cairo_pdf, width = 7, height = 10, family = "Helvetica")
}

# Plot and save known marker genes on UMAP
PlotAndSaveKnownMarkerGenesOnUMAP <- function(so, keep, marker_genes, name, display=T, raster_dpi=100) {
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
    geom_point_rast(aes(color = exprs), size = 1, raster.dpi = raster_dpi) +
    scale_color_gradientn(colors = jdb_palette("solar_extra")) +
    pretty_plot() +
    facet_wrap(~symbol, ncol = 4)
  if (display) {
    plot(p)
  }
  ggsave(p + theme(legend.position = "none"), filename = paste0("../plots/", name, "/umap_knownmarkers.pdf"), device = cairo_pdf, width = 7, height = 10, family = "Helvetica")
}

# Save global features
SaveGlobalFeatures <- function(so, name, compress=T) {
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
  
  # Compress if directed
  if (compress) {
    system(paste0("gzip ../features/", name, "/projected_pcaloadings.txt"))
    system(paste0("gzip ../features/", name, "/projected_icaloadings.txt"))
  }
}

# Define fast t-test function for sparse matrices
row_t_welch2 <- function (x, y, alternative = "two.sided", mu = 0, conf.level = 0.95) 
{
  force(x)
  force(y)
  if (is.vector(x)) 
    x <- matrix(x, nrow = 1)
  if (is.vector(y)) 
    y <- matrix(y, nrow = 1)
  if (is.data.frame(x) && all(sapply(x, is.numeric))) 
    x <- data.matrix(x)
  if (is.data.frame(y) && all(sapply(y, is.numeric))) 
    y <- data.matrix(y)
  if (nrow(y) == 1 & nrow(x) > 1) {
    y <- matrix(y, nrow = nrow(x), ncol = ncol(y), byrow = TRUE)
  }
  matrixTests:::assert_equal_nrow(x, y)
  if (length(alternative) == 1) 
    alternative <- rep(alternative, length.out = nrow(x))
  matrixTests:::assert_character_vec_length(alternative, 1, nrow(x))
  choices <- c("two.sided", "less", "greater")
  alternative <- choices[pmatch(alternative, choices, duplicates.ok = TRUE)]
  matrixTests:::assert_all_in_set(alternative, choices)
  if (length(mu) == 1) 
    mu <- rep(mu, length.out = nrow(x))
  matrixTests:::assert_numeric_vec_length(mu, 1, nrow(x))
  matrixTests:::assert_all_in_closed_interval(mu, -Inf, Inf)
  if (length(conf.level) == 1) 
    conf.level <- rep(conf.level, length.out = nrow(x))
  matrixTests:::assert_numeric_vec_length(conf.level, 1, nrow(x))
  matrixTests:::assert_all_in_closed_interval(conf.level, 0, 1)
  mxs <- rowMeans(x, na.rm = F)
  mys <- rowMeans(y, na.rm = F)
  mxys <- mxs - mys
  nxs <- rep.int(ncol(x), nrow(x))
  nys <- rep.int(ncol(y), nrow(y))
  nxys <- nxs + nys
  ####
  # Rewrite to work with sparse matrices
  x@x <- x@x - mxs[x@i+1]
  empty_sum_sq <- (dim(x)[2] - tabulate(x@i + 1, dim(x)[1])) * mxs^2
  vxs <- (rowSums(x^2, na.rm = F) + empty_sum_sq) / (nxs - 1)
  y@x <- y@x - mxs[y@i+1]
  empty_sum_sq <- (dim(y)[2] - tabulate(y@i + 1, dim(y)[1])) * mys^2
  vys <- (rowSums(y^2, na.rm = F) + empty_sum_sq) / (nys - 1)
  #vxs <- rowSums((x - mxs)^2, na.rm = F)/(nxs - 1)
  #vys <- rowSums((y - mys)^2, na.rm = F)/(nys - 1)
  ####
  stderxs <- vxs/nxs
  stderys <- vys/nys
  stders <- stderxs + stderys
  dfs <- stders * stders/(stderxs * stderxs/(nxs - 1) + stderys * 
                            stderys/(nys - 1))
  stders <- sqrt(stders)
  tres <- matrixTests:::do_ttest(mxys, mu, stders, alternative, dfs, conf.level)
  w1 <- nxs < 2
  matrixTests:::showWarning(w1, "had less than 2 \"x\" observations")
  w2 <- !w1 & nys < 2
  matrixTests:::showWarning(w2, "had less than 2 \"y\" observations")
  w3 <- stders <= 10 * .Machine$double.eps * pmax(abs(mxs), 
                                                  abs(mys))
  matrixTests:::showWarning(w3, "were essentially constant")
  tres[w1 | w2 | w3, ] <- NA
  rnames <- rownames(x)
  if (!is.null(rnames)) 
    rnames <- make.unique(rnames)
  data.frame(obs.x = nxs, obs.y = nys, obs.tot = nxys, mean.x = mxs, 
             mean.y = mys, mean.diff = mxys, var.x = vxs, var.y = vys, 
             stderr = stders, df = dfs, statistic = tres[, 1], pvalue = tres[, 
                                                                             2], conf.low = tres[, 3], conf.high = tres[, 4], 
             alternative = alternative, mean.null = mu, conf.level = conf.level, 
             stringsAsFactors = FALSE, row.names = rnames)
}

# Define fast t-test 1 vs all function
FindAllMarkers2 <- function(so, clusters, type) {
  cluster_levels <- unique(so@meta.data[,clusters])
  cluster_levels <- cluster_levels[!is.na(cluster_levels)]
  lapply(
    X = 1:length(cluster_levels),
    #X = 1:3,
    FUN = function(x) {
      print(x)
      ind1 <- so@meta.data[,clusters] == cluster_levels[x]
      ind2 <- so@meta.data[,clusters] != cluster_levels[x]
      ind1[is.na(ind1)] = F
      ind2[is.na(ind2)] = F
      set1 <- so@assays$RNA@data[,ind1]
      set2 <- so@assays$RNA@data[,ind2]
      if(type == "sparse") {
        out <- row_t_welch2(set1, set2)
      }
      else {
        out <- row_t_welch(set1 %>% data.matrix(), set2 %>% data.matrix())
      }
      out <- out %>%
        rownames_to_column(var = "ENSG") %>%
        dplyr::mutate(cluster = as.numeric(x),
                      avg_logFC = log2(mean.x) - log2(mean.y)) %>%
        dplyr::rename("p_val" = pvalue, "tstat" = statistic) %>%
        dplyr::select(cluster, ENSG, avg_logFC, tstat, p_val) %>%
        dplyr::mutate(avg_logFC = ifelse(is.na(avg_logFC), 0, avg_logFC),
                      tstat = ifelse(is.na(tstat), 0, tstat),
                      p_val = ifelse(is.na(p_val), 0, p_val))
      return(out)
    }
  ) %>%
    bind_rows() %>%
    dplyr::mutate(p_val = ifelse(p_val < 10^-300, 10^-300, p_val),
                  p_val_adj = p.adjust(p_val, method = "BH"))
}

# Big within-cluster features function
WithinClusterFeatures <- function(so, clusters, clus, name, suffix="", compress=T, type = "dense") {
  # Calculate differentially expressed genes
  #markers <- FindAllMarkers(so, test.use = "t", logfc.threshold = 0)
  markers <- FindAllMarkers2(so, clusters, type)
  markers <- markers %>%
    merge(., keep  %>% dplyr::select(ENSG, symbol), by = "ENSG") %>%
    #dplyr::mutate(p_val = ifelse(p_val < 10^-200, 10^-200, p_val),
    #              tstat = qt(p_val * 2, dim(so)[2] - 2, lower.tail = F) * sign(avg_logFC)) %>%
    as.tibble()
  #Get upregulated
  demarkers <- markers %>%
    dplyr::filter(p_val_adj < 0.05, avg_logFC > log(2))
  demarkers %>%
    group_by(cluster) %>%
    count()
  demarkers.mat <- demarkers %>%
    dplyr::mutate(value = 1) %>%
    dplyr::select(cluster, ENSG, value) %>%
    dplyr::mutate(cluster = paste0("Cluster", cluster)) %>%
    cast_sparse(ENSG, cluster, value)
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
    dplyr::select(cluster, ENSG, value) %>%
    dplyr::mutate(cluster = paste0("Cluster", cluster)) %>%
    cast_sparse(ENSG, cluster, value)
  notkeep <- keep %>%
    dplyr::filter(ENSG %ni% row.names(demarkers_down.mat))
  missing <- sparseMatrix(dims = c(dim(notkeep)[1], length(colnames(demarkers_down.mat))), i={}, j={})
  row.names(missing) <- notkeep$ENSG
  colnames(missing) <- colnames(demarkers_down.mat)
  demarkers_down.df <- rbind(demarkers_down.mat, missing) %>%
    data.frame()
  #Format continuous
  markers.mat <- markers %>%
    dplyr::select(cluster, ENSG, tstat) %>%
    dplyr::mutate(cluster = paste0("Cluster", cluster)) %>%
    cast_sparse(ENSG, cluster, tstat)
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
  
  # Compress if directed
  if (compress) {
    system(paste0("gzip ../features/", name, "/projected_pcaloadings_clusters", suffix, ".txt"))
    system(paste0("gzip ../features/", name, "/average_expression", suffix, ".txt"))
    system(paste0("gzip ../features/", name, "/diffexprs_genes_clusters", suffix, ".txt"))
    system(paste0("gzip ../features/", name, "/diffexprs_down_genes_clusters", suffix, ".txt"))
    system(paste0("gzip ../features/", name, "/diffexprs_tstat_clusters", suffix, ".txt"))
  }

  return(demarkers)
}

# Plot top DE genes on UMAP
PlotAndSaveDEGenesOnUMAP <- function(so, demarkers, name, suffix = "", height = 10, rank_by_tstat = FALSE, display = T, raster_dpi = 100) {
  if (rank_by_tstat) {
    topdegenes <- demarkers %>%
      group_by(cluster) %>% 
      top_n(n = 2, wt = tstat)
    topdegenes <- topdegenes %>%
      group_by(cluster) %>% 
      slice(c(1,2), with_ties = F)
    topdegenes.df <- bind_cols(data.frame(t(so@assays$RNA@scale.data[topdegenes$ENSG,])),
                               data.frame(so@reductions$umap@cell.embeddings)) %>%
      as.tibble()
  } else {
    topdegenes <- demarkers %>%
      group_by(cluster) %>% 
      top_n(n = 2, wt = avg_logFC)
    topdegenes <- topdegenes %>%
      group_by(cluster) %>% 
      slice(c(1,2), with_ties = F)
    topdegenes.df <- bind_cols(data.frame(t(so@assays$RNA@scale.data[topdegenes$ENSG,])),
                               data.frame(so@reductions$umap@cell.embeddings)) %>%
      as.tibble()
  }
  p <- topdegenes.df %>%
    reshape2::melt(id.vars = c("UMAP_1", "UMAP_2")) %>%
    merge(., keep, by.x = "variable", by.y = "ENSG") %>%
    as.tibble() %>%
    dplyr::mutate(value = case_when(value > 3 ~ 3,
                                    value < -3 ~ -3,
                                    T ~ value)) %>%
    dplyr::rename("exprs" = value) %>%
    ggplot(., aes(x = UMAP_1, y = UMAP_2)) + 
    geom_point_rast(aes(color = exprs), size = 1, raster.dpi = raster_dpi) +
    scale_color_gradientn(colors = jdb_palette("solar_extra")) +
    pretty_plot() +
    facet_wrap(~symbol, ncol = 4)
  if (display) {
    plot(p)
  }
  ggsave(p + theme(legend.position = "none"), filename = paste0("../plots/", name, "/umap_degenes", suffix, ".pdf"), device = cairo_pdf, width = 7, height = height, family = "Helvetica")
}

