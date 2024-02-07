args = commandArgs(trailingOnly = TRUE)
####################################################
##################### STEP 4 #######################
####################################################
f = function(RootDir) {
  for (dir in RootDir) {
    message('Merging data in ', dir)
    rdata0 = NULL
    fs = list.files(
      path = dir,
      pattern = "\\.Rdata",
      full.names = TRUE,
      recursive = FALSE,
      ignore.case = TRUE
    )
    if (length(fs) < 1)
      next
    for (f in fs) {
      if (!grepl(pattern = "\\.Rdata", x = f))
        next
      message('\t file: ', f)
      load(f)
      rdata0 = rbind(rdata0, rdata)
      rm(rdata)
    }
    outDir = file.path(getwd(), "Rdata")
    if (!dir.exists(outDir)) {
      dir.create(outDir, recursive = TRUE)
    }
    save(rdata0, file = file.path(outDir                              ,
                                  paste0(
                                    # Sys.Date(),
                                    #'_',
                                    unique(rdata0$procedure_group),
                                    '.Rdata',
                                    collapse = '_'
                                  )))
    rm(rdata0)
    gc()
  }
}

###########
filepath <- as.character(args[1])
if (!is.null(filepath) &&
      length(filepath) > 0 &&
      filepath != "") {
  message("Parquet2Rdata for each procedure in progress ...")
  f(filepath)
}
###########
# Do not forget to clean up scatter rdata files
###########
