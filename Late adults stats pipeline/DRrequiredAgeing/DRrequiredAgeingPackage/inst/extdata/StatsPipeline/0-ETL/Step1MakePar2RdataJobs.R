# Jobs creator for Par2Rdata
f <- function(path = getwd(), mem = 25000, pattern = ".parquet") {
  files <- list.files(
    path,
    pattern = pattern,
    full.names = TRUE,
    recursive = TRUE,
    include.dirs = FALSE
  )
  write(
    paste0(
      "bsub -J IMPC_stats_pipeline_lsf_jobs -M ", mem,
      " -e ", files, ".err",
      " -o ", files, ".log Rscript Step2Parquet2Rdata.R ",
      files
    ),
    file = "jobs_step2_Parquet2Rdata.bch"
  )
  write("", file = "step1_completed.log")
}
