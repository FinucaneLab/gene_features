library(rjson)
source("make_features_argparser.R")

extract <- function(l, indicator, default) {
  if (startsWith(l, indicator)) {
    gsub(indicator, "", l)
  } else {
    default
  }
}
# New defaults
newConfig <- list(overrideDefaults=list(
  outputDir="../",
  geneAnnot="../resources/gene_annot_jun10.txt",
  conversionDir="../resources"
))
to_skip <- c(
  "_old_code_to_new_config.R",
  "utils.R",
  "make_features_fun.R",
  "make_features_argparser.R",
  "make_features_batch.R",
  "make_features.R"
)
for (file in dir()) {
  if (endsWith(file, ".R") && !(file %in% to_skip)) {
    scriptLines <- readLines(file)
    name <- gsub("\\.R", "", file)
    isProcessed <- FALSE
    generateBatchId <- FALSE
    useIntegrationSampleSizeReference <- FALSE
    markerGenes <- c()
    for (line in scriptLines) {
      lineTrim <- gsub(" ", "", line)
      if (startsWith(lineTrim, "marker_genes<-c"))
        markerGenes <- strsplit(gsub("\"", "", gsub(")", "", gsub("marker_genes<-c\\(", "", lineTrim))), ",")[[1]]
      numberPcs <- extract(lineTrim, "number_pcs<-", 20)
      varGenes <- extract(lineTrim, "vargenes<-", 3000)
      clusRes <- extract(lineTrim, "clus_res<-", 0.8)
      inputAnnot <- "TODO"  # Not easily parseable this way, do manually
      rowAnnot <- "TODO" # Not easily parseable this way, do manually
      colAnnot <- "TODO" # Not easily parseable this way, do manually
      isProcessed <- (isProcessed || grepl("subset=nFeature_RNA>quantile\\(", lineTrim))
      rowIdType <- gsub("\")", "", extract(lineTrim, "mat<-ConvertToENSGAndProcessMatrix\\(mat, \"", "human_symbol\")"))
      minFeatures <- extract(lineTrim, "so<-CreateSeuratObject\\(counts=mat,project=name,min.features=", "200)")
      minFeatures <- if (endsWith(minFeatures, ")")) gsub(")", "", minFeatures) else strsplit(minFeatures, ",")[[1]][1]
      DEGenesPlotHeight <- strsplit(extract(lineTrim, "PlotAndSaveDEGenesOnUMAP\\(so,demarkers,name,height=", "30,"), ",")[[1]][1]
      generateBatchId <- (generateBatchId || grepl("colnames\\(mat.annot)=c\\(\"BATCH_ID", lineTrim))
      integrationAnchorDims <- strsplit(extract(lineTrim, "so.anchors<-FindIntegrationAnchors\\(object.list=so.list,dims=1:", "30,"), ",")[[1]][1]
      useIntegrationSampleSizeReference <- (useIntegrationSampleSizeReference || grepl(",reference=reference_dataset", lineTrim))
    }
    newConfig[[name]] <- list(
      inputData="TODO",
      numberPcs=as.numeric(numberPcs),
      varGenes=as.numeric(varGenes),
      clusres=as.numeric(clusRes),
      inputAnnot=inputAnnot,
      rowAnnot=rowAnnot,
      colAnnot=colAnnot,
      isProcessed=isProcessed,
      rowIdType=rowIdType,
      minFeatures=as.numeric(minFeatures),
      DEGenesPlotHeight=as.numeric(DEGenesPlotHeight),
      generateBatchId=generateBatchId,
      integrationAnchorDims=as.numeric(integrationAnchorDims),
      useIntegrationSampleSizeReference=useIntegrationSampleSizeReference,
      markerGenes=markerGenes
    )
  }
}

write(toJSON(newConfig, indent=4),"../INPROGRESS_public_data_config.json")