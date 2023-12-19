# Load libraries
suppressPackageStartupMessages(library(tidyverse, verbose=FALSE))
suppressPackageStartupMessages(library(data.table, verbose=FALSE))
suppressPackageStartupMessages(library(BuenColors, verbose=FALSE))
suppressPackageStartupMessages(library(Seurat, verbose=FALSE))
suppressPackageStartupMessages(library(irlba, verbose=FALSE))
suppressPackageStartupMessages(library(Matrix, verbose=FALSE))
suppressPackageStartupMessages(library(future, verbose=FALSE))
suppressPackageStartupMessages(library(reticulate, verbose=FALSE))
suppressPackageStartupMessages(library(ggrastr, verbose=FALSE))
suppressPackageStartupMessages(library(tidytext, verbose=FALSE))
suppressPackageStartupMessages(library(matrixTests, verbose=FALSE))
suppressPackageStartupMessages(library(pathviewr, verbose=FALSE))
suppressPackageStartupMessages(library(clustree, verbose=FALSE))
source("utils.R")
source("make_features_argparser.R")

# Possible names indicating predefined clusters
alternate_PreDefinedClusters_names <- c(
  "Main_cell_type", "ctype", "cell.type", "cell.types", "celltype", "CellType", "Cell_type_refined",
  "cell_ontology_class", "cell_annotation", "Cell.Labels", "subtype", "structure_name", "ClusterID", "Sub_Cluster",
  "Clusters_full", "Cluster", "Manuscript_Identity", "group", "Labels", "class", "cluster_id", "TISSUE_ID", "Organ",
  "x"
)

# Looks a bit hacky maybe, but this allows defaults to be defined in 1 place to prevent conflicts
#  when running interactive (or sourced) vs via command line
make_features <- function(name, inputData,
                          logger            = simpleLogger(paste0(name, ".log")),
                          cores             = parser$defaults[[which(parser$args == "--cores")]],
                          numberPcs         = parser$defaults[[which(parser$args == "--numberPcs")]],
                          varGenes          = parser$defaults[[which(parser$args == "--varGenes")]],
                          clusRes           = parser$defaults[[which(parser$args == "--clusRes")]],
                          inputAnnot        = parser$defaults[[which(parser$args == "--inputAnnot")]],
                          rowAnnot          = parser$defaults[[which(parser$args == "--rowAnnot")]],
                          colAnnot          = parser$defaults[[which(parser$args == "--colAnnot")]],
                          isProcessed       = FALSE, # flag default is always FALSE
                          rowIdType         = parser$defaults[[which(parser$args == "--rowIdType")]],
                          minFeatures       = parser$defaults[[which(parser$args == "--minFeatures")]],
                          plotClustree      = FALSE,
                          dontCompress      = FALSE,
                          DEGenesPlotHeight = parser$defaults[[which(parser$args == "--DEGenesPlotHeight")]],
                          generateBatchId   = FALSE,
                          display           = FALSE,
                          integrationAnchorDims = parser$defaults[[which(parser$args == "--integrationAnchorDims")]],
                          useIntegrationSampleSizeReference = FALSE,
                          markerGenes       = NULL,
                          outputDir         = parser$defaults[[which(parser$args == "--outputDir")]],
                          geneAnnot         = parser$defaults[[which(parser$args == "--geneAnnot")]],
                          conversionDir     = parser$defaults[[which(parser$args == "--conversionDir")]]) {

  inputData <- file_ref_to_vec(logger, inputData)
  rowAnnot <- file_ref_to_vec(logger, rowAnnot, alternative_strings = c("none", "header"))
  colAnnot <- file_ref_to_vec(logger, colAnnot, alternative_strings = c("none", "header"))
  # log input arguments
  fun_args <- within(as.list(environment()), rm(logger))

  logger$debug(paste0("make_features function called with the following arguments:\n  ",
                      paste0(lapply(seq_along(fun_args), function(y, n, i) { paste(n[[i]], "=", as.character(y[[i]])) }, y=fun_args, n=names(fun_args)), collapse="\n  ")
  ))
  verify_input(fun_args, logger = logger)

  if (numberPcs != "elbow")
    numberPcs <- as.numeric(numberPcs)
  compress <- !dontCompress
  if (!(endsWith(outputDir, "/")))
    outputDir <- paste0(outputDir, "/")
  convertToEnsg <- conversionDir != ""
  if ((conversionDir != "") && !endsWith(conversionDir, "/"))
    conversionDir <- paste0(conversionDir, "/")

  # Set up parallelization
  # Remember to use htop to delete forgotten forks
  Sys.setenv(R_FUTURE_FORK_ENABLE = T)
  options(future.globals.maxSize = 2048 * 1024^2)
  # This should (emphasis on should) remove the need to manually kill workers...
  #  never hurts to check of course
  oplan <- plan(strategy = "multicore", workers = cores)
  on.exit(plan(oplan), add=TRUE)

  # Setup
  dir.create(paste0(outputDir, "plots/", name), recursive = TRUE, showWarnings = FALSE)
  dir.create(paste0(outputDir, "features/", name), recursive = TRUE, showWarnings = FALSE)
  dir.create(paste0(outputDir, "data/", name), recursive = TRUE, showWarnings = FALSE)

  #------------------------------------------------LOAD AND FORMAT DATA-----------------------------------------------#
  if (geneAnnot != "") {
    logger$trace_expr(keep <- read.table(geneAnnot, sep = "\t", header = T, stringsAsFactors = F, col.names = c("ENSG", "symbol", "chr", "start", "end", "TSS")))
  } else {
    keep <- NULL  # Added later
  }

  # Read in data and annotations
  logger$info("Loading data and annotations")
  if (length(inputData) > 1) {# Read all matrices into a list
    # assumes rownames to be present as first column in each file
    if (length(rowAnnot) == 1)
      rowAnnot <- rep(rowAnnot, length(inputData))
    if (length(colAnnot) == 1)
      colAnnot <- rep(colAnnot, length(inputData))

    logger$trace_expr(datalist <- lapply(seq_along(inputData),
                                         function(i)
                                           read_sparse_mat(inputData[i], keep = keep, rowannot = rowAnnot[i], colannot = colAnnot[i], converttoensg = convertToEnsg, conversiondir = conversionDir, rowIdType = rowIdType))
    )
    # Bind assuming same row order
    logger$trace_expr(mat <- do.call("cbind", datalist))
    colnames(mat) <- make.names(colnames(mat), unique=T)
  } else {
    if (endsWith(tolower(inputData), ".rds")) {
      logger$info("Reading RDS file")
      logger$trace_expr(mat <- get_from_rds(inputData))
    } else if (rowAnnot != "none") {
      # MM format:
      logger$trace_expr(mat <- readMM(inputData) %>%
        data.matrix() %>%
        Matrix(sparse = TRUE))
      logger$trace_expr(col.annot <- data.frame(fread(colAnnot, header=F), row.names=1))
      logger$trace_expr(row.annot <- data.frame(fread(rowAnnot, header=F), row.names=1))
      colnames(mat) <- rownames(col.annot)
      rownames(mat) <- rownames(row.annot)
    } else {
      logger$trace_expr(mat <- data.frame(fread(inputData, sep = "\t"))[,-1] %>%
        data.matrix() %>%
        Matrix(sparse = TRUE))
      logger$trace_expr(rownames(mat) <- data.frame(fread(inputData), select = 1, skip = 1, sep = "\t")[,1])
    }
  }

  # read annotations, if any
  if (all(inputAnnot != "none")) {
    if (all(rowAnnot != "none")) {
      # MM format:
      logger$trace_expr(mat.annot <- data.frame(fread(inputAnnot, header=T)))
      logger$trace_expr(mat.annot <- mat.annot[!duplicated(mat.annot$V1), ])
      rownames(mat.annot) <- mat.annot[,1]
      logger$trace_expr(mat <- mat[,colnames(mat) %in% rownames(mat.annot)])
    } else if (file.exists(inputAnnot)) {
      logger$trace_expr(mat.annot <- data.frame(fread(inputAnnot), row.names=1, header=T))
    } else if (inputAnnot == "header") {
      logger$trace_expr(mat.annot <- data.frame(Id=dimnames(mat)[[2]]))
      logger$trace_expr(mat.annot$PreDefinedClusters <- sapply(mat.annot$Id, function(i) last(strsplit(i, "\\.")[[1]])))
    }
    if (!("PreDefinedClusters" %in% colnames(mat.annot))) {
      for (alt_name in alternate_PreDefinedClusters_names)
        if (alt_name %in% colnames(mat.annot)) {
          colnames(mat.annot)[colnames(mat.annot) == alt_name] <- "PreDefinedClusters"
          break
        }
    }
  }

  # Convert to ENSG, drop duplicates, and fill in missing genes
  if (convertToEnsg)
    logger$trace_expr(mat <- ConvertToENSGAndProcessMatrix(mat, row_id_type = rowIdType, keep=keep, conversiondir = conversionDir))

  #--------------------------------------------------COMPUTE FEATURES-------------------------------------------------#
  logger$info("Converting to Seurat object")
  if (all(inputAnnot != "none")) {
    # Replace invalid dimname characters like R would to ensure column names and rownames match
    logger$trace_expr(rownames(mat.annot) <- gsub("-", ".", rownames(mat.annot)))
    logger$trace_expr(colnames(mat) <- gsub("-", ".", colnames(mat)))
    logger$trace_expr(so <- CreateSeuratObject(counts = mat, project = name, min.features = minFeatures, meta.data = mat.annot))
  } else if (generateBatchId) {
    mat.annot <- matrix(0, nrow = dim(mat)[2], ncol = 1)
    curr_batch_start_ind <- 1
    for (i in seq_along(datalist)) {
      num_cells <- dim(datalist[[i]])[2]
      mat.annot[curr_batch_start_ind:(curr_batch_start_ind + num_cells - 1), 1] <- toString(i)
      curr_batch_start_ind <- curr_batch_start_ind + num_cells
    }
    mat.annot <- data.frame(mat.annot)
    rownames(mat.annot) <- colnames(mat)
    colnames(mat.annot) <- c("BATCH_ID")
    logger$trace_expr(so <- CreateSeuratObject(counts = mat, project = name, min.features = minFeatures, meta.data = mat.annot))

  } else {
    logger$trace_expr(so <- CreateSeuratObject(counts = mat, project = name, min.features = minFeatures))
  }
  # Clean up
  rm(mat)
  gc(verbose=FALSE, full=TRUE) # May take a few seconds, but gest to always GC after rm of big things

  # QC
  if (!isProcessed){
    logger$info("Running QC")
    logger$trace_expr(so <- subset(so,
                                   subset = nFeature_RNA > quantile(so$nFeature_RNA, 0.05) &
                                     nFeature_RNA < quantile(so$nFeature_RNA, 0.95)))
  }

  if (length(inputData) > 1) {
    # Regress out batch effects
    splitbys <- c("BATCH_ID", "Platform", "donor", "Seq_batch_ID")
    splitby <- splitbys[splitbys %in% colnames(so@meta.data)]
    if (length(splitby) == 0)  {
      logger$error(paste0("Tried splitting Seurat by any of ", paste0(splitbys, collapse=","), " but failed"))
    } else if (length(splitby) > 1) {
      logger$error(paste0("Tried splitting Seurat by any of ", paste0(splitbys, collapse=","), " but found multiple"))
    }
    logger$trace_expr(so.list <- SplitObject(so, split.by = splitby))
    for (i in seq_along(so.list)) {
      logger$trace_expr(so.list[[i]] <- NormalizeData(so.list[[i]]))
      logger$trace_expr(so.list[[i]] <- FindVariableFeatures(so.list[[i]], nfeatures = varGenes))
    }
    if (useIntegrationSampleSizeReference) {
      logger$trace_expr(so.anchors <- FindIntegrationAnchors(object.list = so.list, dims = 1:integrationAnchorDims, k.filter=100, anchor.features = vargenes, reference = which.max(lapply(so.list, FUN = function(x) dim(x)[2]))))
    } else {
      logger$trace_expr(so.anchors <- FindIntegrationAnchors(object.list = so.list, dims = 1:integrationAnchorDims, k.filter=100, anchor.features = varGenes))
    }
    logger$trace_expr(so.integrated <- IntegrateData(anchorset = so.anchors, dims = 1:integrationAnchorDims, features.to.integrate = rownames(so)))
    # This is a bit of a hack, which makes it so that we don't have to change our pipeline
    logger$trace_expr(so <- CreateSeuratObject(CreateAssay5Object(counts = GetAssayData(so, assay = "RNA",layer = "counts"), data = GetAssayData(so.integrated,assay = "integrated", layer = "data")), meta.data = so@meta.data))


    # This is never normalized, only scaled
  } else {
    logger$info("Normalizing data")
    logger$trace_expr(so <- NormalizeData(so, normalization.method = "LogNormalize", scale.factor = 1000000))
  }
  logger$info("Scaling data")
  logger$trace_expr(so <- ScaleData(so))

  if (length(inputData) == 1) {
    # Identify variable genes
    logger$info("Identifying variable genes")
    logger$trace_expr(so <- FindVariableFeatures(so, nfeatures = varGenes))
    logger$trace_expr(PlotAndSaveHVG(so, outputDir, name, display=display))
  } else {
    VariableFeatures(so) <- so.integrated@assays$integrated@var.features
  }

  if (!is.numeric(numberPcs)) { # if 'elbow':
    logger$info("Determining number of PCs to use")
    # Run maximal PCA
    logger$trace_expr(max_pcs <- min(c(length(row.names(so)), length(colnames(so))))-1)
    logger$trace_expr(so2 <- RunPCA(so, npcs = max_pcs))  # dropped: , features = row.names(so)
    # Find elbow and use row number as number of PCs
    logger$trace_expr(numberPcs <- find_curve_elbow(cbind(1:max_pcs, so2@reductions$pca@stdev)))
    logger$info(paste("Number of PCs to use:", numberPcs))
    logger$trace_expr(min_dev <- so2@reductions$pca@stdev[numberPcs])
    logger$trace_expr(PlotAndSavePCAElbow(so2, max_pcs, outputDir, name, hline=min_dev, display=display))
    # Redo PCA with chosen number of PCs:
  }
  logger$info("Running PCA")
  # Run PCA
  logger$trace_expr(so <- RunPCA(so, npcs = numberPcs))  # features can't be dropped here
  # Project PCA to all genes
  logger$trace_expr(so <- ProjectDim(so, do.center = TRUE))
  # Plot Elbow
  logger$trace_expr(PlotAndSavePCAElbow(so, numberPcs, outputDir, name, display=display))


  # Run ICA
  logger$info("Running ICA")
  ica_skipped <- FALSE
  try_out <- try({
    logger$trace_expr(so <- RunICA(so, nics = numberPcs))
    # Project ICA to all genes
    logger$trace_expr(so <- ProjectDim(so, reduction = "ica", do.center = T))
  }, silent=TRUE)
  if (class(try_out) == "try-error") {
    error_type <- attr(try_out,"condition")
    if (error_type$message == "unused argument (features = NULL)") {
      ica_skipped <- TRUE
      logger$info("Caught expected error, skipping ICA")
    } else {
      stop(error_type)
    }
  }

  logger$info("Clustering")
  # Cluster cells
  logger$trace_expr(so <- FindNeighbors(so, dims = 1:numberPcs, nn.eps = 0))

  # Make clustree plot
  if (plotClustree) {
    logger$info("Generating clustree plot")
    so2 <- so
    magnitude <- floor(log10(clusRes))  # Get order of magnitude
    # Resolutions to try for plotting purposes, 5 steps of magnitude lower and 5 higher, e.g. with plot_res 0.6 this will be:
    # 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1
    # only values above 0 are kept, so e.g. with clusRes 3 this will be
    # 1 2 3 4 5 6 7 8
    res_to_try <- c(clusRes-(5:1)*10^(magnitude), clusRes, clusRes+(1:5)*10^(magnitude))
    res_to_try <- res_to_try[res_to_try > 0]
    for (res in res_to_try) {
      logger$debug(paste0("Running res ",res," for clustree"))
      logger$trace_expr(so2 <- FindClusters(so2, resolution = res, n.start = 100))
    }
    logger$trace_expr(p <- clustree(so2, prefix="RNA_snn_res."))
    # add a red bar at chosen clusRes
    p <- p + geom_hline(yintercept = which(res_to_try[length(res_to_try):1] == clusRes)-1, colour="red")  # Add bar for chosen RES
    ggsave(p, filename = paste0(outputDir, "plots/", name, "/clustree.pdf"), device = cairo_pdf, width = 6, height = 4, family = "Helvetica")
    rm (so2)
    logger$info("Continuing base clustering")
  }

  logger$trace_expr(so <- FindClusters(so, resolution = clusRes, n.start = 100))

  # UMAP dim reduction
  logger$info("Running UMAP dimension reduction")
  logger$trace_expr(so <- RunUMAP(so, reduction.key = "UMAP_",dims = 1:numberPcs, min.dist = 0.4, n.epochs = 500,
                                  n.neighbors = 10, learning.rate = 0.1, spread = 2))
  # Plot UMAP clusters
  logger$trace_expr(PlotAndSaveUMAPClusters(so, so@meta.data$seurat_clusters, outputDir, name, display=display))

  # Plot known clusters on UMAP (if applicable)
  if (all(inputAnnot != "none")) {
    logger$trace_expr(
      PlotAndSaveUMAPClusters(so, so@meta.data$PreDefinedClusters, outputDir, name, display=display, suffix = "_pre_def")
    )
  }

  # Hack to keep all features
  if (geneAnnot == "")
    logger$trace_expr(keep <- data.frame(ENSG=rownames(so@assays$RNA)))

  # Plot PCs on UMAP
  logger$trace_expr(PlotAndSavePCsOnUMAP(so, outputDir, name, number_pcs = numberPcs, display=display))
  # Plot ICs on UMAP
  if (!ica_skipped) logger$trace_expr(PlotAndSaveICsOnUMAP(so, outputDir, name, display=display))
  # Plot known marker genes on UMAP
  if (length(markerGenes[markerGenes != ""]) > 0)
    logger$trace_expr(PlotAndSaveKnownMarkerGenesOnUMAP(so, keep, markerGenes, outputDir, name, display=display))

  # Save global features
  logger$info("Saving features")
  logger$trace_expr(SaveGlobalFeatures(so, keep, outputDir, name, compress=compress, ica_skipped))

  if (length(inputData) > 1) {
    rownames(so@assays$RNA@layers$scale.data) <- so.integrated@assays$integrated@var.features
  } else {
    rownames(so@assays$RNA@layers$scale.data) <- rownames(so@assays$RNA@features)
  }

  # Compute any cluster dependent features (DE genes, within-cluster PCs, etc.) and save them
  logger$info("Computing dependent features")
  # Seurat clusters
  Idents(object=so) <- "seurat_clusters"
  clus <- levels(so@meta.data$seurat_clusters)
  logger$trace_expr(demarkers <- WithinClusterFeatures(so, keep, clusters="seurat_clusters", clus=clus, outputdir=outputDir, name=name, compress=compress))
  # Pre-defined cluster dependent features (if applicable)
  if (inputAnnot != "none") {
    Idents(object=so) <- "PreDefinedClusters"
    clus <- unique(so@meta.data$PreDefinedClusters)
    logger$trace_expr(demarkers_pre_def <- WithinClusterFeatures(so, keep, clusters="PreDefinedClusters", clus=clus, outputdir=outputDir, name=name, suffix = "_pre_def", compress=compress))
  }

  logger$info("Plotting dependent features")
  # Plot DE genes on UMAP
  logger$trace_expr(PlotAndSaveDEGenesOnUMAP(so, keep, demarkers, outputDir, name, display=display, height = DEGenesPlotHeight, rank_by_tstat = TRUE))
  # Plot DE genes from pre-defined clusters on UMAP (if applicable)
  if (inputAnnot != "none") {
    logger$trace_expr(PlotAndSaveDEGenesOnUMAP(so, keep, demarkers_pre_def, outputDir, name, display=display, suffix = "_pre_def", height = DEGenesPlotHeight, rank_by_tstat = TRUE))
  }

  # Save Seurat object
  logger$info("Saving Seurat object")
  logger$trace_expr(saveRDS(so, paste0(outputDir, "data/", name, "/so.rds")))
  logger$info("Finished")
}