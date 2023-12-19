#!/usr/bin/env Rscript
source("make_features_argparser.R")
source("make_features_fun.R")

# Parse arguments and start logger
config <- parse_args(parser)

# Parse possible vectors:
for (argName in c("inputData", "inputAnnot", "rowAnnot", "colAnnot"))
  if (grepl(",", config[[argName]]))
    config[[argName]] <- strsplit(config[[argName]], ",")[[1]]

config$logger <- simpleLogger(paste0(config$outputDir, "/logs/", config$name, ".log"),
                              loglevel="TRACE",
                              printlevel="INFO")

do.call(make_features, within(config, rm(help, opts))[-1])
