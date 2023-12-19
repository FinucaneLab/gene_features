#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(rjson, verbose=FALSE))
source("make_features_argparser.R")
source("make_features_fun.R")

# For interactive_testing:
#  Will setup the environment for the indicated name (logging to interactive.log)
INTERACTIVE_TEST <- FALSE
interactive_name <- "human_kidney_test"

batch_parser <- arg_parser("Generate and format batch of PoPS feature matrices. see make_features.R for options and defaults. Please turn comma-separated arguments into lists (see batch_config.json)")
batch_parser <- add_argument(batch_parser, "--config", help="config JSON file", type="character")
batch_parser <- add_argument(batch_parser, "--name", help="batch name", type="character", default="batch")
batch_parser <- add_argument(batch_parser, "--log", help="logfile", type="character", default="batch")
# Layout of the json should be as follows:
# {
#  name : { optionName: Value, optionName2: Value2 },
#  name2 : { optionName: Value, optionName2: Value2 }
# }
# Defaults will be used for options not included
# custom defaults can be set with the 'overrideDefaults' name, e.g.:
# {
#  "overrideDefaults": {"cores": 6}
#  name : { optionName: Value, optionName2: Value2 },
#  name2 : { optionName: Value, optionName2: Value2 }
# }
# See ../batch_config.json for an example

if (INTERACTIVE_TEST) {
  config <- list(config="../batch_config.json", name="interactive")
} else {
  config <- parse_args(batch_parser)
}

logger <- simpleLogger(paste0(config$log),
                       loglevel="TRACE",
                       printlevel="INFO")

defaults <- lapply(parser$args[-1], function(i) {
  default <- parser$defaults[[which(parser$args == i)]]
  if (is.na(default) && parser$is.flag[[which(parser$args==i)]])
    FALSE
  else
    default
})
names(defaults) <- gsub("--", "", parser$args[-1])
defaults <- within(defaults, rm(help, opts, name, inputData)) # Remove argparser arguments and required ones from defaults

# Parse config JSON and prepare arguments list
batch_config <- fromJSON(file=config$config)
if ("overrideDefaults" %in% names(batch_config)) {
  if ("inputData" %in% batch_config[["overrideDefaults"]])
    logger$error("inputData cannot be set as default")
  for (key in names(batch_config[["overrideDefaults"]])) {
    if (!(key %in% names(defaults)))
      logger$error(paste0("key \"", key,"\" in overrideDefaults not recognized"))
    logger$debug(paste0("Setting default: ", key," = ", batch_config[["overrideDefaults"]][[key]]))
    defaults[[key]] <- batch_config[["overrideDefaults"]][[key]]
  }
}

fun_args_list <- list()
for (runName in names(batch_config)[names(batch_config) != "overrideDefaults"]) {
  fun_args_list[[runName]] <- defaults
  fun_args_list[[runName]]$name <- runName
  fun_args_list[[runName]]$logger <- logger
  if (!("inputData" %in% names(batch_config[[runName]])))
    logger$error(paste0("Error in ", runName, ": inputData is a required argument"))
  for (key in names(batch_config[[runName]])) {
    if (!((key %in% names(defaults)) || (key == "inputData")))
      logger$error(paste0("Error in ", runName, ": ", key," is a not a recognized argument"))
    fun_args_list[[runName]][[key]] <- batch_config[[runName]][[key]]
  }
}

logger$info("Finished parsing JSON, starting batch")

if (INTERACTIVE_TEST) {
  lapply(seq_along(fun_args_list[[interactive_name]]), function(i) assign(names(fun_args_list[[interactive_name]])[i], fun_args_list[[interactive_name]][[i]], envir = .GlobalEnv))
} else {
  for (runName in names(fun_args_list)) {
    logger$change_prefix(runName)
    logger$info(paste0("Starting generation of features for ", runName))
    tryCatch(
    {do.call(make_features, fun_args_list[[runName]])},
    error = function(cond) logger$error(paste0("FAILED: ", conditionMessage(cond)), terminate=FALSE)
    )

  }
}
