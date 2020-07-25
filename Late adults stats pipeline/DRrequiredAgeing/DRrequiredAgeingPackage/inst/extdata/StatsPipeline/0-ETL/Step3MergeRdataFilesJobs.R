f = function(path = file.path(getwd(), 'ProcedureScatterRdata'),
             mem = 40000) {
  dirs = list.dirs(path = path,
                   full.names = TRUE ,
                   recursive  = FALSE)
  ############## jobs creator#
  write(
    paste0(
      'bsub -M ',
      mem,
      ' -e "step4_MergeRdatas_error.log" -o "step4_MergeRdatas_output.log" Rscript Step4MergingRdataFiles.R "',
      unique(na.omit(dirs)),
      '"'
    ),
    file = 'jobs_step4_MergeRdatas.bch'
  )
  write('',file = 'Step3 completed.log')
}

