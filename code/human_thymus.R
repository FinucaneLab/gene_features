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
options(future.globals.maxSize = 4096 * 1024^2)
plan(strategy = "multicore", workers = 32)

# Parameters
name <- "human_thymus"
number_pcs <- 60
vargenes <- 2500
clus_res <- 0.4

# Setup
dir.create(paste0("../plots/", name))
dir.create(paste0("../features/", name))

# Notes on data:
# Stuck: can't read in h5ad file (duplicate row names error)

#------------------------------------------------LOAD AND FORMAT DATA-----------------------------------------------#

# Read in data and annotations
file <- paste0("../data/", name, "/HTA08.v01.A06.Science_human_tcells.h5ad")
tmp_so <- ReadH5AD(file, assay = "RNA", layers = "data")

