args = commandArgs(trailingOnly = TRUE)
library(foreach)
#library(DRrequired)
library(DRrequiredAgeing)
library(doParallel)
library(parallel)
library(foreach)
library(SmoothWin)
library(nlme)
library(base64enc)
library(jsonlite)
#library(rjson)


jobCreator = function(path = getwd(),
                      pattern = '.Rdata',
                      JobListFile = 'DataGenerationJobList.bch') {
  if (!dir.exists(file.path('DataGeneratingLog')))
    dir.create(
      path = file.path('DataGeneratingLog'),
      recursive = TRUE,
      showWarnings = FALSE
    )
  if (file.exists(JobListFile))
    unlink(JobListFile)
  files = list.files(
    path = path,
    pattern = pattern,
    full.names = TRUE,
    recursive = FALSE,
    include.dirs = FALSE,
    all.files = TRUE
  )
  proc = tools::file_path_sans_ext(basename(files))
  write(
    paste0(
      'bsub -J IMPC_stats_pipeline_lsf_jobs ',
      ' -e ',
      file.path('DataGeneratingLog', paste0(proc, '_errorlog.log')),
      ' -o ',
      file.path('DataGeneratingLog', paste0(proc, '_outputlog.log')),
      ' -n 1 -q bigmem -M 80000 Rscript InputDataGenerator.R "',
      files,
      '" "',
      proc,
      '"'
    ),
    file = JobListFile
  )
}

#jobCreator()

GenerateData = function(args, thresh = 4) {
  message('Reading the solr query for the procedures (procedure_group) ...',
          args[2])
  # Sys.sleep(20)
  # while (system(command = 'bjobs -r | wc -l',
  # 							wait = TRUE,
  # 							intern = TRUE) > thresh) {
  # 	message(Sys.time() + '. waiting for 20s ...',
  # 					args[2])
  # 	Sys.sleep(20)
  # }
  library(DRrequiredAgeing)
  trash = mainAgeing(
    file = args[1],
    subdir = args[2],
    coreRatio = 0 / 14,
    seed = 123456,
    BatchProducer = TRUE,
    activateMulticore = FALSE,
    cpu = 1,
    memory = 8000,
    controlSize = 1500,
    extraBatchParameters = NULL,
    combineEAandLA = FALSE,
    solrBaseURL = NULL
  )
  trash = NULL
  gc()
  return(NULL)
}
GenerateData(args)
