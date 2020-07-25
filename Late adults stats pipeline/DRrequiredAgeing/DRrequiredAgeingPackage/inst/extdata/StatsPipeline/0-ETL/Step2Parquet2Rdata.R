args = commandArgs(trailingOnly = TRUE)
####################################################
##################### STEP 1 #######################
####################################################
f = function(files) {
  library(miniparquet)
  df = lapply(seq_along(files),
              function(i) {
                message(i, '|', length(files), ' ~> ', files[i])
                r = miniparquet::parquet_read(files[i])
                return(r)
              })
  ###############
  procedure_list = na.omit(unique(unlist(lapply(df, function(x) {
    unique(x$procedure_group)
  }))))

  for (proc in procedure_list) {
    message(which(procedure_list == proc),
            '/',
            length(procedure_list),
            ' Checking for ',
            proc)
    rdata = data.table::rbindlist(lapply(df, function(x) {
      r = subset(x, x$procedure_group == proc)
      return(r)
    }))

    if (!is.null(rdata) &&
        nrow    (rdata) > 0) {
      rdata = rdata[!duplicated(rdata),]
    }

    outDir = file.path('ProcedureScatterRdata', proc)
    if (!dir.exists(outDir))
      dir.create(outDir, recursive = TRUE)
    save(rdata,
         file = file.path(
           outDir,
           paste(
             round(runif(1) * 10 ^ 6)  ,
             proc                  ,
             '.Rdata'              ,
             sep = '_'             ,
             collapse = '_'
           )
         ),
         compress = FALSE)
    rm(rdata)
    gc()
  }
}

##############
if(!is.null(args[1])   &&
   length(args[1]) > 0 &&
   args[1] != '') {
  message('Parquet2Rdata for each procedure in progress ...')
  f(args[1])
}
###############




