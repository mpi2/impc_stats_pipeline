args = commandArgs(trailingOnly = TRUE)

runner = function(file) {
  library(DRrequiredAgeing)
  mainAgeing(
    file = suppressWarnings(tail(UnzipAndfilePath(file), 1)),
    subdir = 'Results_IMPC_SP_Windowed',
    concurrentControlSelect = FALSE,
    seed = 123456,
    # For windowing only,
    messages = FALSE,
    outdelim = '\t',
    encode = FALSE,
    BatchProducer = FALSE,
    # Windowing
    activeWindowing = TRUE,
    check = 1,
    storeplot = TRUE,
    plotWindowing = TRUE,
    ####
    debug = FALSE,
    noSpaceAllowed = TRUE,
    ####
    virtualDrive = FALSE,
    ####
    MMOptimise = c(1,1,1,1,1,1),
    FERRrep = 0, # decision on 19-8-2020 12.00 AM to turn this option off
    ### Just for the certain situations (must be removed in FUTURE)
    OverwriteExistingFiles = FALSE,
    storeRawData = TRUE,
    outlierDetection = FALSE,
    compressRawData = TRUE,
    writeOutputToDB = FALSE,
    solrBaseURL = '',
    measureExecutionTime = FALSE,
    onlyFillNotExisitingResults = FALSE,
    DRversion = 'DRversionNotSpecified'
  )
}

ignore.my.name = runner(file = args[1])
