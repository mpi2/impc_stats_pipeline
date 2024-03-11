args <- commandArgs(trailingOnly = TRUE)
library(foreach)
library(DRrequiredAgeing)
library(SmoothWin)
library(nlme)
library(base64enc)
library(jsonlite)

generate_data <- function(args, thresh = 4) {
  message("Reading the solr query for the procedures (procedure_group)...",
          args[2])
  library(DRrequiredAgeing)
  trash <- mainAgeing(
    file = args[1],
    subdir = args[2],
    seed = 123456,
    BatchProducer = TRUE,
    cpu = 1,
    memory = 8000,
    controlSize = 1500,
    extraBatchParameters = NULL,
    combineEAandLA = FALSE,
    solrBaseURL = NULL
  )
  trash <- NULL
  gc()
  return(NULL)
}
generate_data(args)
