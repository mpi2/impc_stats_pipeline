# Jobs creator for Par2Rdata
f <- function(path = getwd(),
              pattern = ".parquet",
              mem = "5G",
              time = "00:05:00") {
  files <- list.files(
    path,
    pattern = pattern,
    full.names = TRUE,
    recursive = TRUE,
    include.dirs = FALSE
  )
  write(
    paste0(
      "sbatch --job-name=impc_stats_pipeline_job --mem=", mem,
      " --time=", time,
      " -e ", files, ".err",
      " -o ", files, ".log --wrap='Rscript Step2Parquet2Rdata.R ",
      files,
      "'"
    ),
    file = "jobs_step2_Parquet2Rdata.bch"
  )
  write("", file = "../step1_completed.log")
}
