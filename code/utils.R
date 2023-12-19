# Not %in% function
"%ni%" <- Negate("%in%")

# t-test
DiffTTest <- function (object, cells.1, cells.2, genes.use = NULL, print.bar = TRUE,
          assay.type = "RNA")
{
  data.test <- GetAssayData(object = object, assay.type = assay.type,
                            layer = "data")
  genes.use <- SetIfNull(x = genes.use, default = rownames(x = data.test))
  if (print.bar) {
    iterate.fxn <- pblapply
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
ConvertToENSGAndProcessMatrix <- function(mat, row_id_type, keep, conversiondir) {
  ### Convert to ENSG
  if (row_id_type == "mouse_symbol") {
    sym2ensmusg <- read.table(paste0(conversiondir, "symbol2ensmusg.txt"), header = F, stringsAsFactors = F, col.names = c("symbol", "ENSMUSG"))
    ### Drop rows that don't exist in column
    mat <- mat[row.names(mat) %in% sym2ensmusg$symbol,]
    ### Rename
    idx <- match(row.names(mat), sym2ensmusg$symbol)
    row.names(mat) <- sym2ensmusg$ENSMUSG[idx]
    ### Update row type variable
    row_id_type <- "ENSMUSG"
  }
  if (row_id_type == "human_symbol") {
    sym2ensg <- read.table(paste0(conversiondir, "ensg2symbol.txt"), header = F, stringsAsFactors = F, col.names = c("ENSG", "symbol"))
    ### Drop rows that don't exist in column
    mat <- mat[row.names(mat) %in% sym2ensg$symbol,]
    ### Rename
    idx <- match(row.names(mat), sym2ensg$symbol)
    row.names(mat) <- sym2ensg$ENSG[idx]
    ### Update row type variable
    row_id_type <- "ENSG"
  }
  if (row_id_type == "ENSMUSG") {
    ensmusg2ensg <- read.table(paste0(conversiondir, "ensmusg2ensg.txt"), header = F, stringsAsFactors = F, col.names = c("ENSMUSG", "ENSG"))
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
  ### Filter to genes of interest
  if (is.data.frame(keep)) {
    rowkeep <- row.names(mat) %in% keep$ENSG
    mat <- mat[rowkeep,]
    ### Add missing genes
    notkeep <- keep %>%
      dplyr::filter(ENSG %ni% row.names(mat))
    missing <- sparseMatrix(dims = c(dim(notkeep)[1], length(colnames(mat))), i={}, j={})
    row.names(missing) <- notkeep$ENSG
    colnames(missing) <- colnames(mat)
    mat <- rbind(mat, missing)
  }
  ### Return
  return(mat)
}

# Plot and save variable genes
PlotAndSaveHVG <- function(so, outputdir, name, display = T) {
  p <- so@assays$RNA@meta.data %>%
    dplyr::rename("variable_gene" = vf_vst_counts_variable) %>%
    ggplot(.) +
    geom_histogram(aes(x = log10(vf_vst_counts_variance.standardized), fill = variable_gene), bins = 50) +
    pretty_plot() +
    scale_fill_manual(values = jdb_palette("FantasticFox")[c(3,5)]) +
    xlab("log10 mean-standardized variance")
  if (display) {
    plot(p)
  }
  ggsave(p, filename = paste0(outputdir, "plots/", name, "/variablegenes.pdf"), device = cairo_pdf, width = 6, height = 4, family = "Helvetica")
}

# Plot and save PCA Elbow
PlotAndSavePCAElbow <- function(so, num_pcs, outputdir, name, display = T, hline=-1) {
  p <- cbind("PC" = seq(1:num_pcs), "stdev" = so@reductions$pca@stdev) %>%
    as_tibble() %>%
    ggplot(., aes(x =  PC, y = stdev)) +
    geom_point() +
    pretty_plot() +
    ylab("standard deviation") +
    xlab("principal components")
  if (hline != -1) {
    filename <- paste0(outputdir, "plots/", name, "/maxpcaelbow.pdf")
    p <- p + geom_hline(yintercept=hline, colour = "red")
  } else {
    filename <- paste0(outputdir, "plots/", name, "/pcaelbow.pdf")
  }
  if(display) {
    plot(p)
  }
  ggsave(p, filename = filename, device = cairo_pdf, width = 6, height = 4, family = "Helvetica")
}

# Plot and save UMAP Clusters
PlotAndSaveUMAPClusters <- function(so, clust_col, outputdir, name, suffix = "", display = T, raster_dpi = 300, width = 6, height = 5) {
  clusters.df <- bind_cols(data.frame("cluster" = clust_col),
                           data.frame(so@reductions$umap@cell.embeddings)) %>%
    as_tibble()
  p <- clusters.df %>%
    ggplot(., aes(x = UMAP_1, y = UMAP_2)) +
    geom_point_rast(aes(color = cluster), size = 0.5, raster.dpi = raster_dpi) +
    scale_color_manual(values = c(jdb_palette("corona"), jdb_palette("corona"), jdb_palette("corona"), jdb_palette("corona"), jdb_palette("corona"))) +
    pretty_plot() +
    guides(col = guide_legend(ncol = 1)) +
    theme(legend.key.size = unit(0.1, 'lines'))
  if(display) {
    plot(p)
  }
  ggsave(p + theme(legend.position = "right"), filename = paste0(outputdir, "plots/", name, "/umap_clusters", suffix, ".pdf"), device = cairo_pdf, width = width, height = height, family = "Helvetica")
}

# Plot and save PCs on UMAP
PlotAndSavePCsOnUMAP <- function(so, outputdir, name, display=T, raster_dpi=100, number_pcs=24) {
  if (number_pcs > 24) {
    number_pcs <- 24
  }
  pcs.df <- bind_cols(data.frame(so@reductions$pca@cell.embeddings[,1:number_pcs]),
                      data.frame(so@reductions$umap@cell.embeddings)) %>%
    as_tibble()
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
  ggsave(p + theme(legend.position = "none"), filename = paste0(outputdir, "plots/", name, "/umap_pcs.pdf"), device = cairo_pdf, width = 7, height = 10, family = "Helvetica")
}

# Plot and save ICs on UMAP
PlotAndSaveICsOnUMAP <- function(so, outputdir, name, display=T, raster_dpi=100) {
  ics.df <- bind_cols(data.frame(so@reductions$ica@cell.embeddings[,1:24]),
                      data.frame(so@reductions$umap@cell.embeddings)) %>%
    as_tibble()
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
  ggsave(p + theme(legend.position = "none"), filename = paste0(outputdir, "plots/", name, "/umap_ics.pdf"), device = cairo_pdf, width = 7, height = 10, family = "Helvetica")
}

# Plot and save known marker genes on UMAP
PlotAndSaveKnownMarkerGenesOnUMAP <- function(so, keep, marker_genes, outputdir, name, display=T, raster_dpi=100) {
  knownmarker_genes <- keep %>%
    dplyr::select(ENSG, symbol) %>%
    dplyr::filter(symbol %in% marker_genes)
  knownmarkers.df <- bind_cols(data.frame(t(so@assays$RNA$scale.data[knownmarker_genes$ENSG,])),
                               data.frame(so@reductions$umap@cell.embeddings)) %>%
    as_tibble()
  p <- knownmarkers.df %>%
    reshape2::melt(id.vars = c("UMAP_1", "UMAP_2")) %>%
    merge(., knownmarker_genes, by.x = "variable", by.y = "ENSG") %>%
    as_tibble() %>%
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
  ggsave(p + theme(legend.position = "none"), filename = paste0(outputdir, "plots/", name, "/umap_knownmarkers.pdf"), device = cairo_pdf, width = 7, height = 10, family = "Helvetica")
}

# Save global features
SaveGlobalFeatures <- function(so, keep, outputdir, name, compress=T, ica_skipped=FALSE) {
  # Write out projected gene loadings across all cells
  so@reductions$pca@feature.loadings.projected %>%
    data.frame() %>%
    rownames_to_column(., var = "ENSG") %>%
    as_tibble() %>%
    arrange(factor(ENSG, levels = keep$ENSG)) %>%
    fwrite(.,
           paste0(outputdir, "features/", name, "/projected_pcaloadings.txt"),
           quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
  if (compress)
    system(paste0("gzip -f ", outputdir, "features/", name, "/projected_pcaloadings.txt"))
  if (!ica_skipped) {
    # Write out ICA across all cells
    so@reductions$ica@feature.loadings.projected %>%
      data.frame() %>%
      rownames_to_column(., var = "ENSG") %>%
      as_tibble() %>%
      arrange(factor(ENSG, levels = keep$ENSG)) %>%
      fwrite(.,
             paste0(outputdir, "features/", name, "/projected_icaloadings.txt"),
             quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
    if (compress)
      system(paste0("gzip -f ", outputdir, "features/", name, "/projected_icaloadings.txt"))
  }
}

# Define fast t-test function for sparse matrices
row_t_welch2 <- function (x, y, alternative = "two.sided", mu = 0, conf.level = 0.95) {
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
  rnames <- dimnames(x)[[1]]
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
    X = seq_along(cluster_levels),
    #X = 1:3,
    FUN = function(x) {
      ind1 <- so@meta.data[,clusters] == cluster_levels[x]
      ind2 <- so@meta.data[,clusters] != cluster_levels[x]
      ind1[is.na(ind1)] <- FALSE
      ind2[is.na(ind2)] <- FALSE
      set1 <- so@assays$RNA@layers$data[,ind1]
      set2 <- so@assays$RNA@layers$data[,ind2]
      if(type == "sparse") {
        out <- row_t_welch2(set1, set2)
      } else {
        out <- row_t_welch(set1 %>% data.matrix(), set2 %>% data.matrix())
      }
      if (nrow(out) == nrow(so@assays$RNA@layers$scale.data)) {
        rownames(out) <- rownames(so@assays$RNA@layers$scale.data)
      } else {
        # This happens when variable genes are manually specified (i.e. with multiple inputfiles)
        rownames(out) <- rownames(so@assays$RNA@features)
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
WithinClusterFeatures <- function(so, keep, clusters, clus, outputdir, name, suffix="", compress=T, type = "dense") {

  # Calculate differentially expressed genes
  #markers <- FindAllMarkers(so, test.use = "t", logfc.threshold = 0)
  markers <- FindAllMarkers2(so, clusters, type)
  # For some reason FindAllMarkers2 makes ENSG dissappear
  if (all(c("ENSG", "symbol") %in% colnames(markers)))
    markers <- markers %>%
      merge(., keep  %>% dplyr::select(ENSG, symbol), by = "ENSG") %>%
      #dplyr::mutate(p_val = ifelse(p_val < 10^-200, 10^-200, p_val),
      #              tstat = qt(p_val * 2, dim(so)[2] - 2, lower.tail = F) * sign(avg_logFC)) %>%
      as_tibble()
  else
    markers <- markers %>%
      merge(., keep  %>% dplyr::select(ENSG), by = "ENSG") %>%
      #dplyr::mutate(p_val = ifelse(p_val < 10^-200, 10^-200, p_val),
      #              tstat = qt(p_val * 2, dim(so)[2] - 2, lower.tail = F) * sign(avg_logFC)) %>%
      as_tibble()
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

  so.clus.pcs <- lapply(seq_along(clus), function(x) {PC_within_cluster(so, clus[x])})
  so.clus.pcs[sapply(so.clus.pcs, is.null)] <- NULL
  so.clus.pcs <- so.clus.pcs %>%
    do.call(cbind, .)

  # Write out projected gene loadings within clusters
  so.clus.pcs %>%
    data.frame() %>%
    rownames_to_column(., var = "ENSG") %>%
    as_tibble() %>%
    arrange(factor(ENSG, levels = keep$ENSG)) %>%
    fwrite(.,
           paste0(outputdir, "features/", name, "/projected_pcaloadings_clusters", suffix, ".txt"),
           quote = F, row.names = F, col.names = T, sep = "\t")

  # Write out normalized expression within clusters and across all cells
  so.ae <- AverageExpression(so, layer = "scale.data")$RNA
  colnames(so.ae) <- paste0("Cluster", colnames(so.ae))
  # so.ae$Allcells <- apply(so@assays$RNA@layers$scale.data, 1, mean)
  so.ae %>%
    data.frame() %>%
    # apply(., 1, mean) %>%
    # data.frame() %>%
    rownames_to_column(., var = "ENSG") %>%
    as_tibble() %>%
    arrange(factor(ENSG, levels = keep$ENSG)) %>%
    fwrite(.,
           paste0(outputdir, "features/", name, "/average_expression", suffix, ".txt"),
           quote = F, row.names = F, col.names = T, sep = "\t")
  so.ae %>%
    data.frame() %>%
    apply(., 1, mean) %>%
    data.frame() %>%
    rownames_to_column(., var = "ENSG") %>%
    as_tibble() %>%
    arrange(factor(ENSG, levels = keep$ENSG)) %>%
    fwrite(.,
           paste0(outputdir, "features/", name, "/average_expression_rowmeans", suffix, ".txt"),
           quote = F, row.names = F, col.names = T, sep = "\t")


  # Write differential expression (DE genes) between clusters
  demarkers.df %>%
    rownames_to_column(., var = "ENSG") %>%
    as_tibble() %>%
    arrange(factor(ENSG, levels = keep$ENSG)) %>%
    fwrite(.,
           paste0(outputdir, "features/", name, "/diffexprs_genes_clusters", suffix, ".txt"),
           quote = F, row.names = F, col.names = T, sep = "\t")

  # Write differential expression (DE genes) between clusters (downregulated)
  demarkers_down.df %>%
    rownames_to_column(., var = "ENSG") %>%
    as_tibble() %>%
    arrange(factor(ENSG, levels = keep$ENSG)) %>%
    fwrite(.,
           paste0(outputdir, "features/", name, "/diffexprs_down_genes_clusters", suffix, ".txt"),
           quote = F, row.names = F, col.names = T, sep = "\t")

  # Write differential expression (t-stat) between clusters
  markers.df %>%
    rownames_to_column(., var = "ENSG") %>%
    as_tibble() %>%
    arrange(factor(ENSG, levels = keep$ENSG)) %>%
    fwrite(.,
           paste0(outputdir, "features/", name, "/diffexprs_tstat_clusters", suffix, ".txt"),
           quote = F, row.names = F, col.names = T, sep = "\t")

  # Compress if directed
  if (compress) {
    system(paste0("gzip -f ", outputdir, "features/", name, "/projected_pcaloadings_clusters", suffix, ".txt"))
    system(paste0("gzip -f ", outputdir, "features/", name, "/average_expression", suffix, ".txt"))
    system(paste0("gzip -f ", outputdir, "features/", name, "/diffexprs_genes_clusters", suffix, ".txt"))
    system(paste0("gzip -f ", outputdir, "features/", name, "/diffexprs_down_genes_clusters", suffix, ".txt"))
    system(paste0("gzip -f ", outputdir, "features/", name, "/diffexprs_tstat_clusters", suffix, ".txt"))
  }

  return(demarkers)
}

# Plot top DE genes on UMAP
PlotAndSaveDEGenesOnUMAP <- function(so, keep, demarkers, outputdir, name, suffix = "", height = 10, rank_by_tstat = FALSE, display = T, raster_dpi = 100) {
  if (rank_by_tstat) {
    topdegenes <- demarkers %>%
      group_by(cluster) %>%
      top_n(n = 2, wt = tstat)
    topdegenes <- topdegenes %>%
      group_by(cluster) %>%
      slice(c(1,2))  # removed with_ties=F
    topdegenes.df <- bind_cols(data.frame(t(so@assays$RNA@layers$scale.data[topdegenes$ENSG,])),
                               data.frame(so@reductions$umap@cell.embeddings)) %>%
      as_tibble()
  } else {
    topdegenes <- demarkers %>%
      group_by(cluster) %>%
      top_n(n = 2, wt = avg_logFC)
    topdegenes <- topdegenes %>%
      group_by(cluster) %>%
      slice(c(1,2))
    topdegenes.df <- bind_cols(data.frame(t(so@assays$RNA@layers$scale.data[topdegenes$ENSG,])),
                               data.frame(so@reductions$umap@cell.embeddings)) %>%
      as_tibble()
  }

  p <- topdegenes.df %>%
    reshape2::melt(id.vars = c("UMAP_1", "UMAP_2")) %>%
    merge(., keep, by.x = "variable", by.y = "ENSG") %>%
    as_tibble() %>%
    dplyr::mutate(value = case_when(value > 3 ~ 3,
                                    value < -3 ~ -3,
                                    T ~ value)) %>%
    dplyr::rename("exprs" = value) %>%
    ggplot(., aes(x = UMAP_1, y = UMAP_2)) +
    geom_point_rast(aes(color = exprs), size = 1, raster.dpi = raster_dpi) +
    scale_color_gradientn(colors = jdb_palette("solar_extra")) +
    pretty_plot() +
  if ("symbol" %in% colnames(topdegenes.df))
    p <- p + facet_wrap(~symbol, ncol = 4)
  if (display)
    plot(p)

  ggsave(p + theme(legend.position = "none"), filename = paste0(outputdir, "plots/", name, "/umap_degenes", suffix, ".pdf"), device = cairo_pdf, width = 7, height = height, family = "Helvetica", limitsize = FALSE)
}

read_sparse_mat <- function(x, rowannot, colannot, keep, converttoensg, conversiondir, rowIdType="ENSG") {
  if (endsWith(toLower(x), ".rds")) {
    mat <- get_from_rds(x)
  } else if (rowannot  == "none") {
    mat <- data.frame(fread(x), row.names=1) %>%
        data.matrix() %>%
        Matrix(sparse = TRUE)
  } else {
    mat <- readMM(x) %>%
      data.matrix() %>%
      Matrix(sparse = TRUE)
    rows <- data.frame(fread(rowannot, header=FALSE), row.names=1)
    cols <- data.frame(fread(colannot, header=FALSE), row.names=1)
    rownames(mat) <- rownames(rows)
    colnames(mat) <- rownames(cols)
  }
  if (converttoensg)
    mat <- ConvertToENSGAndProcessMatrix(mat, rowIdType, keep, conversiondir)
  ### Reorder according to keep so that we can cbind later
  if (is.data.frame(keep)) mat[match(keep$ENSG, rownames(mat)),] else mat
}

# Simple logger
simpleLogger <- setRefClass("simpleLogger",
                            fields=list(logfile="character", logprefix="character", to_log="vector", to_print="vector"))
simpleLogger$methods(
  initialize=function(logfile_path="log.txt", loglevel="TRACE", printlevel="INFO", logprefix="") {
    # Log prefix can be used if you want to make 1 log file for multiple runs for instance, it's just prefixed on each print
    if (grepl("/", logfile_path)) {
      logdir <- paste0(strsplit(logfile_path, "/")[[1]][1:(str_count(logfile_path, "/"))], collapse="/")
      if (!(dir.exists(logdir)))
        cat("Directory", logdir, "does not exist, creating it instead\n")
        dir.create(logdir, showWarnings = FALSE)
    }
    .self$logfile <- logfile_path
    cat("", file=.self$logfile)
    lvls <- c("ERROR", "WARN", "INFO", "DEBUG", "TRACE")
    if (!(printlevel %in% lvls)) {
      warning("printlevel should be one of c(\"ERROR\", \"WARN\", \"INFO\", \"DEBUG\"), resetting to INFO")
      printlevel <- "INFO"
    }
    .self$to_print <- lvls[1:which(lvls == printlevel)]

    if (!(loglevel %in% lvls)) {
      warning("loglevel should be one of c(\"ERROR\", \"WARN\", \"INFO\", \"DEBUG\", \"TRACE\"), resetting to TRACE")
      loglevel <- "TRACE"
    }
    .self$to_log <- lvls[1:which(lvls == loglevel)]
    .self$logprefix <- logprefix
  },
  format_message=function(level, message) {
    if (logprefix == "") {
      paste(strptime(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(level, ":"), message)
    } else {
      paste(paste(.self$logprefix, "-"), strptime(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(level, ":"), message)
    }
  },
  log = function(msg, lvl, terminate=FALSE) {
    msg_fmt <- .self$format_message(lvl, msg)
    if (lvl %in% .self$to_log) cat(msg_fmt, "\n", file=.self$logfile, append=TRUE)
    if (lvl %in% .self$to_print) {
      if (lvl == "ERROR") {
        if (terminate) stop(msg) else cat(msg_fmt, "\n")
      } else if (lvl == "WARN") {
        if (lvl %in% .self$to_print) warning(msg)
      } else {
        if (lvl %in% .self$to_print) cat(msg_fmt, "\n")
      }
    }
  },
  error=function(msg, terminate=TRUE) .self$log(msg, "ERROR", terminate),
  warn=function(msg)  .self$log(msg, "WARN"),
  info=function(msg)  .self$log(msg, "INFO"),
  debug=function(msg) .self$log(msg, "DEBUG"),
  trace=function(msg) .self$log(msg, "TRACE"),
  change_prefix=function(newprefix) .self$logprefix <- newprefix,
  trace_expr=function(expr) {
    .self$trace(paste0("\n  expression called:\n", paste0("-- ", deparse(substitute(expr)), collapse="\n")))
    expr_stderr <- capture.output(expr_stdout <- capture.output(
      invisible(eval(expr = expr, envir = parent.env(environment())))
      , type="output"), type="message")
    on.exit(
      if ((length(expr_stdout) > 0) || (length(expr_stderr) > 0))
        .self$trace(paste0("\n", paste0("+stdout-> ", expr_stdout, collapse="\n"),
                           "\n", paste0("+stderr-> ", expr_stderr, collapse="\n")))
    )
  }
)
as.character.simpleLogger <- function(x)
  paste0("simpleLogger class logging to ", x$logfile, " with level ", x$to_log,"and printing level ", x$to_print)

verify_input <- function(args, logger) {
  # This might feel extra, but it massively helps in (a) capturing errors early and (b) producing more easy to understand errors
  verify_numeric <- function(logger, x, label, positive=FALSE) {
    if (as.numeric(x) != x)
      logger$error(paste0(label, "(", x, ") invalid, should be numeric"))
    if (positive)
      if (x <= 0)
        logger$error(paste0(label, "(", x, ") invalid, should be greater than 0"))
  }
  verify_int <- function(logger, x, label, positive=FALSE) {
    verify_numeric(logger, x, label, positive)
    if (x != round(x))
      logger$error(paste0(label, "(", x, ") invalid, should be integer"))
  }
  verify_str_of_numeric <- function(logger, x, label, positive=TRUE) {
    if (is.na(suppressWarnings(as.numeric(x))))
      logger$error(paste0(label, "(", x, ") invalid, should be interpretable as numeric"))
    verify_numeric(logger, as.numeric(x), label, positive)
  }
  verify_str_of_int <- function(logger, x, label, positive=TRUE) {
    verify_str_of_numeric(logger, x, label, positive)
    verify_int(logger, as.numeric(x), label, positive)

  }
  verify_bool <- function(logger, x, label, strict=TRUE)
    if (strict) {
      if (!(x %in% c(TRUE, FALSE)))
        logger$error(paste0(label, "(", x, ") invalid, should be boolean"))
    } else {
      if (x || !x)
        logger$error(paste0(label, "(", x, ") invalid, should be interptetable as boolean"))
    }
  verify_file <- function(logger, x, label, or_dir=FALSE) {
    if (!(file.exists(x) || (dir.exists(x) && or_dir)))
      logger$error(paste0(label, ": ", x, ": file or directory not found"))
  }
  verify_int(logger, args$cores, "Number of cores", positive = TRUE)
  if (args$cores > parallel::detectCores())
    logger$warn(paste0("Number of specified cores (",args$cores,") greater than nummber of detected cores (",parallel::detectCores(), "), you've been warned"))
  if (args$numberPcs != "elbow")
    verify_str_of_int(logger, args$numberPcs, "Number of PCs", positive = TRUE)
  verify_int(logger, args$varGenes, "Number of variable genes", positive = TRUE)
  verify_numeric(logger, args$clusRes, "Cluster resolution", positive = TRUE)
  if (length(args$inputData) > 1) {
    for (i in args$inputData)
      verify_file(logger, i, "inputData", or_dir=TRUE)
  } else {
    verify_file(logger, args$inputData, "inputData", or_dir = TRUE)
  }
  if (length(args$inputAnnot) > 1) {
    for (i in args$inputAnnot)
      verify_file(logger, i, "inputAnnot")
  } else {
    if (!(args$inputAnnot %in% c("none", "header")))
      if (!(file.exists(args$inputAnnot)))
        logger$error("Input annotation should be an existing file, 'header' or 'none'")
  }
  if (length(args$rowAnnot) > 1) {
    for (i in args$rowAnnot)
      verify_file(logger, i, "rowAnnot")
  } else {
    if (all(args$rowAnnot != "none"))
      if (!(file.exists(args$rowAnnot)))
        logger$error("Input row annotation should be an existing file or 'none'")
  }
  if (length(args$colAnnot) > 1) {
    for (i in args$colAnnot)
      verify_file(logger, i, "colAnnot")
  } else {
  if (all(args$colAnnot != "none"))
    if (!(file.exists(args$colAnnot)))
      logger$error("Input col annotation should be an existing file or 'none'")
  }
  if (all((args$rowAnnot != "none") | (args$rowAnnot != "none")))
    if (all((args$colAnnot == "none") | (args$rowAnnot == "none")))
      logger$error("row and col annotations should always be used together")
  verify_bool(logger, args$isProcessed, "is processed")
  if (!(is.character(args$rowIdType)))
    logger$error("row ID type should be character")
  verify_int(logger, args$minFeatures, "Minimum number of features", positive = FALSE)
  verify_bool(logger, args$plotClustree, "plot clustree")
  verify_int(logger, args$DEGenesPlotHeight, "DEGenes plot height", positive = TRUE)
  verify_bool(logger, args$dontCompress, "dontCompress")
  verify_bool(logger, args$display, "display")
  verify_bool(logger, args$generateBatchId, "generateBatchId")
  if (args$generateBatchId && (length(args$inputData) == 1))
    logger$error("generateBatchId can only be used with multiple input files")
  verify_int(logger, args$integrationAnchorDims, "integrationAnchorDims", positive = TRUE)
  verify_bool(logger, args$useIntegrationSampleSizeReference, "useIntegrationSampleSizeReference")
  if (args$useIntegrationSampleSizeReference && (length(args$inputData) == 1))
    logger$error("useIntegrationSampleSizeReference can only be used with multiple input files")
  if (!dir.exists(args$outputDir))
    logger$error("outputDir does not exist")
  if (args$geneAnnot != "")
    verify_file(logger, args$geneAnnot, "geneAnnot")
  if (args$conversionDir != "")
    if (!(dir.exists(args$conversionDir)))
      logger$error("conversionDir: directory not found")
}

file_ref_to_vec <- function(logger, x, alternative_strings=c()) {
  # Convert reference to file(s) (either filename, directoryname, or file pattern) to vector of files
  if (length(x) > 1) {
    return(x)
  } else if (file.exists(x)) {
    return(x)
  } else if (dir.exists(x)) {
    return(paste0(gsub("/$", "", x), "/", dir(x)))
  } else if (length(Sys.glob(x)) > 0) {
    return(Sys.glob(x))
  } else if (x %in% alternative_strings) {
    return(x)
  } else {
    logger$error(paste0(deparse(substitute(x)), ": File(s) or directory not found."))
  }
}

get_from_rds <- function(x) {
  pre_load <- ls()
  load(x)
  obj_names <- ls()[!(ls() %in% pre_load)]
  obj_names <- obj_names[obj_names != "pre_load"]
  if (length(obj_names) != 1)  # 2 because pre_load is also new
    logger$error("Cannot handle multiple saved objects currently")
  mat <- get(obj_names) # make copy
  rm(list=obj_names)    # Remove the original
  mat
}