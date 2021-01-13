args = commandArgs(trailingOnly = TRUE)

runner = function(file) {
  library(DRrequiredAgeing)
  mainAgeing(
    file = suppressWarnings(tail(UnzipAndfilePath(file), 1)),
    subdir = 'Results_IMPC_SP',
    concurrentControlSelect = FALSE,
    seed = 123456,
    # For windowing only,
    messages = FALSE,
    outdelim = '\t',
    encode = FALSE,
    BatchProducer = FALSE,
    # Windowing
    activeWindowing = FALSE,
    check = 1,
    storeplot = FALSE,
    plotWindowing = FALSE,
    ####
    debug = FALSE,
    noSpaceAllowed = TRUE,
    ####
    virtualDrive = FALSE,
    ####
    MMOptimise = c(1,1,1,1,1,1),
    FERRrep = 0, # decision on 19-8-2020 12.00 AM to turn this option off
    ####
    activateMulticore = FALSE,
    coreRatio = 1,
    MultiCoreErrorHandling = 'stop',
    inorder = FALSE,
    verbose = TRUE,
    ### Just for the certain situations (must be removed in FUTURE)
    OverwriteExistingFiles = FALSE,
    storeRawData = TRUE,
    outlierDetection = FALSE,
    compressRawData = TRUE,
    writeOutputToDB = FALSE,
    solrBaseURL = '',
    measureExecutionTime = FALSE,
    onlyFillNotExisitingResults = FALSE
  )
}

ignore.my.name = runner(file = args[1])
