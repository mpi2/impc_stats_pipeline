############# jobs creator for Par2Rdata
f = function(path = getwd(), mem = 25000,pattern = '.parquet') {
  files = list.files(
    path = paste0(path, '/'),
    pattern = pattern,
    full.names = TRUE,
    recursive = TRUE,
    include.dirs = FALSE
  )
  write(
    paste0(
      'bsub -M ',
      mem,
      ' -e "step2_Par2Rdata_error.log" -o "Step2_Par2Rdata_output.log" Rscript Step2Parquet2Rdata.R "',
      files,
      '" ""'
    ),
    file = 'jobs_step2_Parquet2Rdata.bch'
  )
  write('',file = 'Step1 completed.log')
}

