# Jobs creator for MergeRdatas
f <- function(path = file.path(getwd(), "ProcedureScatterRdata"),
              mem = "50G",
              time = "01:30:00") {
  dirs <- list.dirs(path = path,
                    full.names = TRUE,
                    recursive  = FALSE)
  unique_dirs <- unique(na.omit(dirs))
  write(
    paste0(
      "sbatch --job-name=impc_stats_pipeline_job --mem=", mem,
      " --time=", time,
      " -e ", unique_dirs, "/step4_merge_rdatas.err",
      " -o ", unique_dirs, "/step4_merge_rdatas.log",
      " --wrap Rscript Step4MergingRdataFiles.R ", unique_dirs
    ),
    file = "jobs_step4_MergeRdatas.bch"
  )
  write("", file = "step3_completed.log")
}
