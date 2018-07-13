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

# UMAP 
RunUMAP <- function(
  object,
  cells.use = NULL,
  dims.use = 1:5,
  reduction.use = 'pca',
  genes.use = NULL,
  assay.use = 'RNA',
  max.dim = 2L,
  reduction.name = "umap",
  reduction.key = "UMAP",
  n_neighbors = 30L,
  min_dist = 0.3,
  metric = "correlation",
  seed.use = 42,
  ...
) {
  if (!py_module_available(module = 'umap')) {
    stop("Cannot find UMAP, please install through pip (e.g. pip install umap-learn).")
  }
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
    py_set_seed(seed = seed.use)
  }
  cells.use <- SetIfNull(x = cells.use, default = colnames(x = object@data))
  if (is.null(x = genes.use)) {
    dim.code <- GetDimReduction(
      object = object,
      reduction.type = reduction.use,
      slot = 'key'
    )
    dim.codes <- paste0(dim.code, dims.use)
    data.use <- GetDimReduction(
      object = object,
      reduction.type = reduction.use,
      slot = 'cell.embeddings'
    )
    data.use <- data.use[cells.use, dim.codes, drop = FALSE]
  } else {
    data.use <- GetAssayData(object = object, assay.type = assay.use, slot = 'scale.data')
    genes.use <- intersect(x = genes.use, y = rownames(x = data.use))
    if (!length(x = genes.use)) {
      stop("No genes found in the scale.data slot of assay ", assay.use)
    }
    data.use <- data.use[genes.use, cells.use, drop = FALSE]
    data.use <- t(x = data.use)
  }
  parameters.to.store <- as.list(x = environment(), all = TRUE)[names(formals("RunUMAP"))]
  object <- SetCalcParams(
    object = object,
    calculation = "RunUMAP",
    ... = parameters.to.store
  )
  umap_import <- import(module = "umap", delay_load = TRUE)
  umap <- umap_import$UMAP(
    n_neighbors = as.integer(x = n_neighbors),
    n_components = as.integer(x = max.dim),
    metric = metric,
    min_dist = min_dist
  )
  umap_output <- umap$fit_transform(as.matrix(x = data.use))
  colnames(x = umap_output) <- paste0(reduction.key, 1:ncol(x = umap_output))
  rownames(x = umap_output) <- cells.use
  object <- SetDimReduction(
    object = object,
    reduction.type = reduction.name,
    slot = "cell.embeddings",
    new.data = as.matrix(x = umap_output)
  )
  object <- SetDimReduction(
    object = object,
    reduction.type = reduction.name,
    slot = "key",
    new.data = reduction.key
  )
  return(object)
}