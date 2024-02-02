# Jobs creator for MergeRdatas
f <- function(path = file.path(getwd(), "ProcedureScatterRdata"), mem = 80000) {
  dirs <- list.dirs(path = path,
                    full.names = TRUE,
                    recursive  = FALSE)
  unique_dirs <- unique(na.omit(dirs))
  write(
    paste0(
      "bsub -q bigmem -J IMPC_stats_pipeline_lsf_jobs -M ", mem,
      " -e ", unique_dirs, "/step4_merge_rdatas.err",
      " -o ", unique_dirs, "/step4_merge_rdatas.log ",
      "Rscript Step4MergingRdataFiles.R ", unique_dirs
    ),
    file = "jobs_step4_MergeRdatas.bch"
  )
  write("", file = "step3_completed.log")
}
