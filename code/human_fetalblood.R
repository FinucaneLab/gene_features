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
name <- "human_fetalblood"
number_pcs <- 60
vargenes <- 2500
clus_res <- 0.2

# Setup
dir.create(paste0("../plots/", name))
dir.create(paste0("../features/", name))

# Notes on data:
# 

#------------------------------------------------LOAD AND FORMAT DATA-----------------------------------------------#

# Load this in in case we need it later
keep <- read.table("../resources/gene_annot_jun10.txt", sep = "\t", header = T, stringsAsFactors = F, col.names = c("ENSG", "symbol", "chr", "start", "end", "TSS"))

# Read in data and annotations
# We'll have to read in each batch and cbind
read_sparse_mat <- function(x) {
  mat <- readMM(x) %>%
    Matrix(sparse = TRUE) %>%
    as("dgCMatrix")
  file <- sub("matrix", "genes", x)
  file <- sub("mtx", "tsv", file)
  rows <- data.frame(fread(file, header=F), row.names=1)
  file <- sub("matrix", "barcodes", x)
  file <- sub("mtx", "tsv", file)
  cols <- data.frame(fread(file, header=F), row.names=1)
  rownames(mat) <- rownames(rows)
  colnames(mat) <- rownames(cols)
  mat <- ConvertToENSGAndProcessMatrix(mat, "ENSG")
  ### Reorder according to keep so that we can cbind later
  mat <- mat[match(keep$ENSG, rownames(mat)),]
  return(mat)
}

filelist <- Sys.glob(paste0("../data/", name, "/*.mtx.gz"))
# Read all matrices into a list
datalist <- lapply(filelist, read_sparse_mat)

# Bind assuming same row order
mat <- do.call("cbind", datalist)
colnames(mat) <- make.names(colnames(mat), unique=T)

# Make batch indicators
mat.annot <- matrix(0, nrow = dim(mat)[2], ncol = 1)
curr_batch_start_ind = 1
for (i in 1:length(datalist)) {
  num_cells <- dim(datalist[[i]])[2]
  mat.annot[curr_batch_start_ind:(curr_batch_start_ind + num_cells - 1), 1] = toString(i)
  curr_batch_start_ind <- curr_batch_start_ind + num_cells
}
mat.annot <- data.frame(mat.annot)
rownames(mat.annot) <- colnames(mat)
colnames(mat.annot) = c("BATCH_ID")

