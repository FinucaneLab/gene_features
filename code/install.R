# Install all the things
if (!require("devtools")) {install.packages("devtools")}
library(devtools)
if (!require("tidyverse")) {install.packages("tidyverse")}
if (!require("data.table")) {install.packages("data.table")}
if (!require("Matrix")) {install.packages("Matrix")}
if (!require("irlba")) {install.packages("irlba")}
if (!require("Seurat")) {install.packages("Seurat")}
if (!require("argparser")) {install.packages("argparser")}
if (!require("pathviewr")) {install.packages("pathviewr")}
if (!require("clustree")) {install.packages("clustree")}
if (!require("rjson")) {install.packages("rjson")}
if (!require("MUDAN")) {devtools::install_github("JEFworks/MUDAN")}
if (!require("BuenColors")) {devtools::install_github("caleblareau/BuenColors")}


