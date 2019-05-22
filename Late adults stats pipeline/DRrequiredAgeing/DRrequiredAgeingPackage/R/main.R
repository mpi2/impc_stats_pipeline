## You must check the 'check' parameter
mainAgeing = function(file = 'http://ves-ebi-d0:8090/mi/impc/dev/solr/experiment/select?q=*%3A*&fq=procedure_stable_id%3AIMPC_ECG_002&rows=590000&wt=csv&indent=true'            ,
                      sep = ','                            ,
                      na.strings = 'NA'                    ,
                      normalisedPhenlist = FALSE           ,
                      subdir = 'Results'                   ,
                      seed = 123456                        ,
                      readCategoriesFromFile  = TRUE       ,
                      OverwriteExistingFiles  = FALSE      ,
                      onlyFillNotExisitingResults = FALSE  ,
                      WhiteListMethods  = NULL             ,
                      # Carefully use this option!
                      # It can remove the entire result file (some colonies in a single output file)
                      # Only for multicore
                      activateMulticore       = TRUE       ,
                      coreRatio = 5 / 14                   ,
                      concurrentControlSelect = FALSE      ,
                      MultiCoreErrorHandling  = 'pass'     ,
                      inorder                 = FALSE      ,
                      verbose                 = TRUE       ,
                      # Only for simulations
                      simulation              = FALSE      ,
                      Simulation.iteration    = 1          ,
                      # Only for windowing
                      activeWindowing = TRUE               ,
                      sensitivity = c(1, 1, 1, 0)          ,
                      pvalThreshold = c(0, 0, 0, 0)        ,
                      check = 2                            ,
                      direction = c(1, 1)                  ,
                      weightORthreshold = 'weight'         ,
                      predFunction = function(m) {
                        predict(m)
                      }                                    ,
                      residFunction = function(m) {
                        resid(m)
                      }                                    ,
                      messages = TRUE                      ,
                      threshold = sqrt(.Machine$double.eps) * 10,
                      outdelim = '\t'                      ,
                      debug = TRUE                         ,
                      encode = FALSE                       ,
                      noSpaceAllowed = TRUE                ,
                      plotWindowing = TRUE                 ,
                      storeplot = TRUE                     ,
                      virtualDrive   = FALSE               ,
                      checkNamesForMissingColNames = TRUE  ,
                      # Raw data
                      storeRawData           = TRUE        ,
                      compressRawData        = TRUE        ,
                      # Only for Batch generator
                      BatchProducer          =  FALSE      ,
                      cpu = 4                              ,
                      memory = 9000                        ,
                      nMax                   = 10000       ,
                      ChunkSize              = 24          ,
                      MinColoniesInChunks    = 32          ,
                      controlSize            = 1500        ,
                      ### Just for debuging
                      superDebug             = FALSE       ,
                      extraBatchParameters   = '-m "rh7-hosts-ebi5-12 rh7-hosts-ebi5-13 rh7-hosts-ebi5-14 rh7-hosts-ebi5-15 rh7-hosts-ebi5-16 rh7-hosts-ebi5-17 rh7-hosts-ebi5-18 rh7-hosts-ebi5-19 rh7-hosts-ebi5-20 rh7-hosts-ebi5-24 rh7-hosts-ebi5-25 rh7-hosts-ebi5-26 rh7-hosts-ebi5-27 rh7-hosts-ebi6-00 rh7-hosts-ebi6-01 rh7-hosts-ebi6-02 rh7-hosts-ebi6-03 rh7-hosts-ebi6-04 rh7-hosts-ebi6-05 rh7-hosts-ebi6-06 rh7-hosts-ebi6-07 rh7-hosts-ebi6-08 rh7-hosts-ebi6-09 rh7-hosts-ebi6-10 rh7-hosts-ebi6-11 rh7-hosts-ebi6-12 rh7-hosts-ebi6-13 rh7-hosts-ebi6-14 rh7-hosts-ebi6-15 rh7-hosts-ebi6-16 rh7-hosts-ebi6-17"',
                      ...) {
  message0('DRrequiredAgeing loaded')
  message0(
    Sys.time(),
    '  ############################################################\n',
    Sys.time(),
    '  % Make sure ExceptionList and CategoryList are updated!    %\n',
    Sys.time(),
    '  % Please check SENSITIVITY for the windowing algorithm     %\n',
    Sys.time(),
    '  % Sensitivy:',
    sensitivity,
    '                                                              \n',
    Sys.time(),
    '  ############################################################',
    ShowTimeTag = FALSE
  )
  message0('Process started ...')
  message0('Machine info:  ', paste(Sys.info(), collapse = ', '))
  message0('Loading dependent packages ...')
  requireNamespace('PhenStat')
  requireNamespace('PhenStatAgeing'     )
  requireNamespace('doParallel')
  requireNamespace('parallel')
  requireNamespace('foreach')
  requireNamespace('SmoothWin')
  requireNamespace('nlme')
  requireNamespace('base64enc')
  requireNamespace('RJSONIO'    )
  requireNamespace('jsonlite'   )
  # Config files
  message0('Loading configuration ...')
  methodmap                      = readConf('MethodMap.conf')
  equationmap                    = readConf('EquationMap.conf')
  CategoryMap                    = readConf('CategoryMap.conf')
  initial                        = readConf('Initialize.conf')
  exceptionList                  = readFile(file = 'ExceptionMap.list')
  #CategoricalCategoryBlackList   = readFile(file = 'CategoricalCategoryBlackList.list')
  # Main subdirectory/working directory
  message0('Preparing the working directory ...')
  cwd = getwd()
  wd  = file.path(cwd,
                  paste(subdir, sep = '_', collapse = '_'))
  dir.create0(wd, recursive = TRUE)
  if (virtualDrive) {
    message0('Creating a virtual drive ... ')
    system('subst U: /D', wait = TRUE)
    system(paste0('subst U: "', wd, '"'), wait = TRUE)
    wd = 'U:'
  }
  message0('Setting the working directory to: \n\t\t ===> ', wd)
  setwd(dir = wd)
  ##################
  set.seed(seed)
  # Read file
  message0('Reading the input file ...\n\t ~>', file)
  if (!file.exists(file))
    message0('File is not local or does not exist!')

  rdata = read.csv(
    file = file                                    ,
    check.names      = checkNamesForMissingColNames,
    sep              = sep                         ,
    na.strings       = na.strings                  ,
    stringsAsFactors = TRUE
  )
  #### Temporary for ageing pipeline only
  # rdataEarly = read.csv(
  #   file = gsub(
  #     pattern = 'LA_',
  #     replacement = '_',
  #     gsub(
  #       pattern     = 'http://ves-ebi-d1.ebi.ac.uk:8988',
  #       replacement = 'http://ves-ebi-d0.ebi.ac.uk:8986',
  #       x = file
  #     )
  #   ),
  #   check.names      = checkNamesForMissingColNames,
  #   sep              = sep                         ,
  #   na.strings       = na.strings                  ,
  #   stringsAsFactors = TRUE
  # )
  # com_cols  = intersect(colnames(rdata), colnames(rdataEarly))
  # rdata     = rbind(rdata[, com_cols], rdataEarly[, com_cols])

  message0('Input file dimentions: ',
           paste0(dim(rdata), collapse  = ', '))
  rdata = rdata[!is.na(rdata$phenotyping_center), ] # Just to remove NA centers
  new.data              = rdata
  new.data              = new.data[order(Date2Integer(new.data$date_of_experiment)), ]
  #########
  new.data$colony_id    = as.character(new.data$colony_id)
  #new.data$colony_id[new.data$biological_sample_group %in% "control"] = NA
  new.data$external_sample_id = as.factor(new.data$external_sample_id)
  ################
  # Start analysis
  ################
  # Initializing cores
  message0('Initialising cores ...')
  crs = cores0(coreRatio = coreRatio, activate  = activateMulticore)
  closeAllConnections()
  registerDoSEQ()
  message0('The detected OS: ', .Platform$OS.type)
  if (.Platform$OS.type == 'windows') {
    cl = makeCluster(crs,
                     outfile = outMCoreLog(wd))

  } else{
    cl = makeForkCluster(crs,
                         outfile = outMCoreLog(wd))
  }
  registerDoParallel(cl, cores = crs)
  # End of multicore initialization
  # Get possible categories for the categorical variables
  message0('Loading the list of possible categories for categorical variables ...')
  CatList = GetPossibleCategories (procedure = NULL, file = readCategoriesFromFile)
  message0('Filtering the dataset in progress ....')
  Strtime      = Sys.time()
  procedures   = as.character(unique(na.omit(new.data$procedure_group)))
  for (procedure in procedures) {
    ###
    n2.9 = base::subset(new.data,  new.data$procedure_group %in% procedures)
    parameters  = as.character(unique(na.omit(n2.9$parameter_stable_id)))
    for (parameter in parameters) {
      FactorLevels = ReadFactorLevelsFromSolr(parameter = parameter, CatList = CatList)
      ### counter starts here ....
      counter   = 1
      outP      = list()
      n3.0 = base::subset(n2.9,  n2.9$parameter_stable_id %in% parameter)
      centers   = as.character(unique(na.omit(n3.0$phenotyping_center)))
      for (center in centers) {
        n3.1     = base::subset(n3.0, n3.0$phenotyping_center %in% center)
        strains  = as.character(unique(na.omit(n3.1$strain_accession_id)))
        for (strain in strains) {
          n3.2  = base::subset(n3.1,  n3.1$strain_accession_id %in% strain)
          metas = as.character(unique(na.omit(n3.2$metadata_group)))
          for (meta in metas) {
            n3.3   = base::subset(n3.2,  n3.2$metadata_group %in% meta)
            n3.3.c = base::subset(n3.3,  n3.3$biological_sample_group %in% 'control')
            n3.3.m = base::subset(n3.3, !(n3.3$biological_sample_group %in% 'control'))
            zygositys = as.character(unique(na.omit(n3.3.m$zygosity)))
            for (zyg in zygositys) {
              n3.3.m_zyg = base::subset(n3.3.m, n3.3.m$zygosity %in% zyg)
              colonys    = as.character(unique(na.omit(n3.3.m_zyg$colony_id)))
              nColonies  = length(colonys)
              if (BatchProducer && nColonies > 0) {
                #nMax                 = 10000
                ChunkSizeFromNumbers = ((nrow(n3.3.c) < nMax) * max(1, round(nColonies /
                                                                               ChunkSize)) +
                                          (nrow(n3.3.c) >= nMax) * nColonies)
                minCol = ((nrow(n3.3.c) < nMax) * MinColoniesInChunks + (nrow(n3.3.c) >=
                                                                           nMax) * 1)
                ColonyChunks         = chunkVector(
                  x = colonys,
                  n = ChunkSizeFromNumbers,
                  min = minCol,
                  activate = (nColonies >= MinColoniesInChunks) &&
                    (nrow(n3.3.c) >= controlSize)
                )
                outpDir  = file.path0(
                  wd,
                  paste0(Sys.Date(), '_', subdir, '_RawData/'),
                  check = FALSE,
                  create = TRUE,
                  IncludedFileName = FALSE
                )
                SubSubDirOrdFileName = RemoveSpecialChars(paste(
                  Sys.Date()             ,
                  #RandomRegardSeed()    ,
                  #procedure             ,
                  parameter              ,
                  center                 ,
                  zyg                    ,
                  strain                 ,
                  meta                   ,
                  collapse = '_'
                ))
                outpfile = file.path0(
                  outpDir,
                  SubSubDirOrdFileName,
                  check = FALSE,
                  create = TRUE,
                  IncludedFileName = TRUE
                )
                mess = paste0(
                  Sys.time(),
                  '. Processed file: ',
                  SubSubDirOrdFileName,
                  '. #Colonies = ',
                  nColonies,
                  ', #Controls = ',
                  nrow(n3.3.c),
                  ', Chunks = ',
                  length(ColonyChunks)
                )
                message0(mess, ShowTimeTag = FALSE)
                write(
                  x = mess,
                  file = paste0(
                    Sys.Date(),
                    '_',
                    subdir,
                    '_DataGenerationLog.log'
                  ),
                  10 ^ 5,
                  append = TRUE
                )
                counter = 1
                for (ChunkedColonies in ColonyChunks) {
                  BatchData     = rbind (subset(n3.3.m_zyg, n3.3.m_zyg$colony_id %in% ChunkedColonies),
                                         n3.3.c)
                  BatchFileName = file.exists0(
                    paste0(
                      outpfile,
                      '_C',
                      length(ColonyChunks),
                      '_',
                      RandomRegardSeed(),
                      '_',
                      counter,
                      '.csv'
                    ),
                    overwrite = OverwriteExistingFiles
                  )
                  if (all(dim(BatchData) > 0)) {
                    write.csv(BatchData,
                              file = BatchFileName,
                              row.names = FALSE)
                    out = BatchGenerator(
                      file                 = BatchFileName    ,
                      dir                  = outpDir          ,
                      procedure            = procedure        ,
                      parameter            = parameter        ,
                      center               = center           ,
                      cpu                  = cpu              ,
                      memory               = memory           ,
                      extraBatchParameters =  extraBatchParameters
                    )
                    write(
                      x = out,
                      file = paste0(outpDir, '/', subdir, '_Batch.bch'),
                      ncolumns = 10 ^ 5,
                      append = TRUE
                    )
                    counter = counter + 1
                  }
                  rm0(c('BatchData', 'BatchFileName'), silent = TRUE)
                }
              } else{
                message0(
                  ' [',
                  paste(
                    procedure,
                    parameter,
                    center   ,
                    strain   ,
                    meta     ,
                    zyg      ,
                    length(colonys),
                    sep = ']~>['
                  ),
                  ']\n'
                )
                ### Single or multiple cores?
                `%activemulticore%` = ifelse (activateMulticore &&
                                                !BatchProducer,
                                              `%dopar%`,
                                              `%do%`)
                if (activateMulticore &&
                    !BatchProducer) {
                  message0('Multicore processing in progress ...')
                } else{
                  message0('Single core processing in progress ...')
                }
                i = 1
                MultiCoreRes = foreach::foreach (
                  i = 1:length(colonys),
                  .packages = c(
                    'PhenStat'    ,
                    'SmoothWin'   ,
                    'base64enc'   ,
                    'nlme'        ,
                    'RJSONIO'     ,
                    'jsonlite'    ,
                    'PhenStatAgeing',
                    'DRrequiredAgeing'
                  ),
                  .errorhandling = c(MultiCoreErrorHandling),
                  .verbose = verbose                        ,
                  .inorder = inorder
                ) %activemulticore% {
                #for (i in  1:length(colonys)){
                  message0('*~*~*~*~*~* ', i, '|', length(colonys), ' *~*~*~*~*~*')
                  for (sim.index in 1:ifelse(simulation, Simulation.iteration, 1)) {
                    # Removing the old objects if exist
                    ObjectsThatMustBeRemovedInEachIteration()
                    # Initialization before starting the analysis
                    note   = list()
                    colony = colonys[i]
                    message0('Current colony: ',colony)

                    n3.4 = base::subset(n3.3.m_zyg,	n3.3.m_zyg$colony_id %in% c(colony))
                    n3.5 = rbind (n3.4, n3.3.c)
                    note = c(note,
                             list(
                               bodyweight_included_in_data = CheckIfNameExistInDataFrame(obj = n3.5,
                                                                                         name = 'weight',
                                                                                         checkLevels = FALSE)
                             ))

                    # Imaginary URLs
                    note$gene_page_url        = GenePageURL       = GenePageURL(n3.5)
                    note$bodyweight_page_url = BodyWeightCurvURL = BodyWeightCurvURL(n3.5)

                    ReadMeTxt       = ReadMe (obj = n3.4, URL = GenePageURL)

                    # Define response column [do not move me!]
                    depVariable = getResponseColumn(n3.5$observation_type)
                    depVar      = depVariable$column
                    message0('Dependent variable: ', depVar)
                    note$response_type       = paste0(depVar,
                                                      '_of_type_',
                                                      paste(depVariable$lbl, sep = '.'))
                    note$observation_type    =
                      if (!is.null(unique(n3.5$observation_type))) {
                        paste(unique(n3.5$observation_type),
                              sep = '~',
                              collapse = '~')
                      } else{
                        NULL
                      }
                    note$data_type           =
                      if (!is.null(unique(n3.5$data_type))) {
                        paste(unique(n3.5$data_type),
                              sep = '~',
                              collapse = '~')
                      } else{
                        NULL
                      }

                    minSampRequired = ifelse(
                      is.ABR(x = parameter),
                      as.numeric(initial$min_ABR_mut_each_sex),
                      ifelse(
                        is.numeric(n3.5[, depVar]),
                        as.numeric(initial$min_num_mut_each_sex),
                        2
                      )
                    )

                    # add missing levels to categorical variables
                    if (!is.numeric(n3.5[, depVar])) {
                      AllLevels = mapLevelsToFactor(levels = levels(n3.5[, depVar]),
                                                    newlevels = FactorLevels$levels)
                      levels(n3.5[, depVar]) = AllLevels$levels
                      #####
                      note = c(note,
                               FactorLevels$note,
                               AllLevels$note)
                      #####
                      SexGenResLevels = min(2 * 2 * length(AllLevels$levels), 4)
                    } else{
                      SexGenResLevels = 4
                    }
                    if (!depVariable$accepted)
                      return('Not a proper dataset!')

                    if (simulation && is.numeric(n3.5[, depVar])) {
                      message0('Simulation in progress ... Round ',
                               sim.index)
                      n3.5_tmp = mimicControls(
                        df = n3.5,
                        removeMutants = (sim.index == 1) ,
                        ArtifLabel  = 'experimental'     ,
                        mutLabel = 'experimental'        ,
                        baselines = 'control'            ,
                        neutralise = TRUE                ,
                        resample = (sim.index != 1)      ,
                        depVariable = depVar             ,
                        sex  = 'sex'                     ,
                        minSampRequired = minSampRequired,
                        SexGenResLevels = SexGenResLevels,
                        indicator       = sim.index      ,
                        plot = superDebug
                      )
                      n3.5 = n3.5_tmp$df
                      note = list(note , simulation_details = n3.5_tmp$note)
                    }

                    # Summary statistics
                    n3.5_summary = SummaryStatisticsOriginal(x = n3.5, depVar = depVar)
                    note         = c(note, n3.5_summary)

                    # Remove zero frequency categories
                    n3.5.1_F_list = RemoveZeroFrequencyCategories(
                      x = n3.5,
                      minSampRequired = minSampRequired,
                      depVar = depVar,
                      totalLevels = SexGenResLevels
                    )
                    n3.5.1 = n3.5.1_F_list$x
                    note   =  c(note, n3.5.1_F_list$note)

                    # Remove var categories
                    n3.5.1_v_list = RemoveZerovarCategories(
                      x = n3.5.1,
                      depVar = depVar,
                      minvar = 0,
                      method = getMethodi(
                        var = parameter,
                        type = ifelse(
                          is.numeric(n3.5.1[, depVar]),
                          'numeric',
                          'charachter'
                        ),
                        methodMap = methodmap
                      )
                    )
                    n3.5.1 = n3.5.1_v_list$x
                    note   = c(note, n3.5.1_v_list$note)

                    OrgSpecIds      = OtherExtraColumns(
                      obj = n3.5,
                      ColNames = c(
                        'external_sample_id',
                        'sex',
                        'biological_sample_group',
                        depVar,
                        'date_of_experiment',
                        'weight'
                      ),
                      names = c(
                        # all lower case
                        'original_external_sample_id',
                        'original_sex',
                        'original_biological_sample_group',
                        'original_response',
                        'original_date_of_experiment',
                        'original_body_weight'
                      )
                    )
                    note = c(note, OrgSpecIds)
                    message0('Creating output directory and file name ...')
                    SubSubDirOrdFileName = file.path0(
                      RemoveSpecialChars(center)    ,
                      RemoveSpecialChars(procedure) ,
                      RemoveSpecialChars(parameter) ,
                      RemoveSpecialChars(colony)    ,
                      RemoveSpecialChars(zyg)       ,
                      RemoveSpecialChars(meta)      ,
                      create = FALSE,
                      check = noSpaceAllowed
                    )
                    FileName = 'output'
                    outDir   = file.path0(
                      wd,
                      SubSubDirOrdFileName,
                      create = TRUE,
                      check  = FALSE,
                      IncludedFileName = FALSE
                    )
                    outpfile = outpfile2 = paste0(outDir,
                                                  '/',
                                                  FileName,
                                                  collapse = '')
                    message0('Output directory: \n \t\t =>=>=> ',
                             outpfile)
                    if (onlyFillNotExisitingResults) {
                      if (any(file.exists(paste(
                        outpfile, c('NotProcessed.tsv', 'Successful.tsv'), sep = '_'
                      )))) {
                        message0('File already exists then skipped!')
                        return(NULL)
                      } else{
                        message0('Result does not exist! Adding in progress ...')
                        rmme = lapply(list.files(dirname(outpfile), full.names = TRUE), function(x) {
                          if (!is.null(x)    &&
                              file.exists(x) &&
                              (
                                grepl(
                                  pattern = '.Rdata',
                                  x = x,
                                  fixed = TRUE
                                ) ||
                                grepl(
                                  pattern = 'Failed_critical_error',
                                  x = x,
                                  fixed = TRUE
                                )
                              )
                          )
                          file.remove(x)
                        })
                        write(outpfile, file = 'DoesNotExists.log', append = TRUE)
                      }
                    }

                    ####
                    if (storeRawData) {
                      # There is a second snippet for the rawdata + weights
                      RawoutputFile = RawoutputFile0 = file.exists0(paste(outpfile2,
                                                                          'rawData.csv',
                                                                          sep = '_'),
                                                                    overwrite  = OverwriteExistingFiles)
                      ReadMeFile  = file.exists0(file.path(dirname(RawoutputFile), 'ReadMe.txt'))
                      message0(
                        'writting the raw data file to disk ... \n \t\t ===> ',
                        RawoutputFile
                      )
                      write.csv(
                        x = n3.5            ,
                        row.names = FALSE   ,
                        file = RawoutputFile
                      )
                      write(
                        x        = ReadMeTxt ,
                        file     = ReadMeFile,
                        ncolumns = 1
                      )
                      if (compressRawData) {
                        comRes = compressFiles(
                          fileList  = c(ReadMeFile , RawoutputFile),
                          dir       = dirname (RawoutputFile)      ,
                          filename  = basename(RawoutputFile)      ,
                          overwrite = OverwriteExistingFiles
                        )
                        if (comRes$status == 0)
                          message0('Compression successful')
                        else
                          message0('Compression failed')
                        RawoutputFile0 = comRes$file
                      }
                    }
                    note$input_file           = relativePath(path = file, reference  = wd)
                    note$output_raw_data_file = relativePath(path = if (storeRawData) {
                      RawoutputFile0
                    } else{
                      NULL
                    },
                    reference = wd)
                    note$read_me_file         = relativePath(path = if (storeRawData &&
                                                                        !compressRawData) {
                      ReadMeFile
                    } else{
                      NULL
                    },
                    reference = wd)
                    ###'
                    isException    = IsInList(
                      item = c(parameter, procedure),
                      list = exceptionList,
                      message = 'Value found in the skip list'
                    )
                    n3.5.2  = droplevels0(n3.5.1)
                    MergLev = MergeLevels(x = n3.5.2[, depVar],
                                          listOfLevelMaps = CategoryMap)
                    n3.5.2[, depVar] = MergLev$x
                    n3.5.2           = droplevels0(n3.5.2[!is.na(n3.5.2[, depVar]),])
                    n3.5.2OnlyKO     = subset(n3.5.2,n3.5.2$biological_sample_group %in% 'experimental')
                    note$relabeled_levels_categorical_variables_only  = MergLev$note
                    if (!is.null(n3.5.2) &&
                        # data.frame is not zero
                        min0(dim(n3.5.2)) > 0 &&
                        # is it  really exist!
                        length(unique(n3.5.2$biological_sample_group)) > 1 &&
                        # include mut and cont
                        min0(table(n3.5.2$biological_sample_group)) >= minSampRequired &&
                        max0(table(n3.5.2OnlyKO$biological_sample_group, n3.5.2OnlyKO$sex)) > 1 &&
                        # include at least 4/2 of each genotype
                        #length(unique(n3.5.2$colony_id)) > 1  &&
                        length(RepBlank(
                          unique(n3.5.2$colony_id),
                          match = c('', NA, 'NA')
                        )) > 1 &&
                        # include 2 colonies (cont & mut)
                        checkGenInCol(n3.5.2) &&
                        # each sex and genotype
                        depVariable$accepted  &&
                        length(na.omit(n3.5.2[, depVar])) > 0 &&
                        # response is not empty!
                        # there must be variation in data
                        NonZeroVariation(n3.5.2[, depVar]) &&
                        !isException &&
                        columnLevelsVariationRadio(dataset = n3.5.2, columnName = depVar) > 0.005 &&
                        RR_thresholdCheck(data = n3.5.2,depVar = depVar,parameter = parameter,methodmap = methodmap)$criteria_result
                    ) {
                      message0('Analysing the dataset in progress ...')
                      message0('Creating PhenList object ...')
                      a = PhenStat::PhenList(
                        n3.5.2,
                        testGenotype             = 'experimental',
                        refGenotype              = 'control',
                        dataset.colname.genotype = 'biological_sample_group',
                        dataset.colname.sex      = 'sex',
                        dataset.colname.weight   = 'weight',
                        dataset.colname.batch    = 'date_of_experiment'
                      )
                      a_summary_before_concurrent = SummaryStatisticsOriginal(
                        x = a@datasetPL,
                        depVar = depVar,
                        sex = 'Sex',
                        genotype = 'Genotype',
                        label = 'phenlist_data_summary_statistics'
                      )
                      note = c(note, a_summary_before_concurrent)
                      #
                      PhenListSpecIds = OtherExtraColumns (
                        obj = a@datasetPL,
                        ColNames = 'external_sample_id',
                        names = 'phenlist_data_spec_ids'
                      )
                      note = c(note, PhenListSpecIds)
                      ### Get method of analysis
                      method = getMethodi(
                        var = parameter,
                        type = ifelse(
                          is.numeric(a@datasetPL[, depVar]),
                          'numeric',
                          'charachter'
                        ),
                        methodMap = methodmap
                      )
                      # WhiteList methods
                      if (!is.null(WhiteListMethods) &&
                          !(method %in% WhiteListMethods)) {
                        message0('Black list applied. Method = ', method)
                        return(FALSE)
                      }
                      #### Check for concurrent control selection
                      aTmp = concurrentContSelect(
                        activate = concurrentControlSelect &&
                          (method %in% 'MM') && !activeWindowing,
                        PhenListObj = a,
                        depVar = depVar,
                        minSampRequired = minSampRequired
                      )
                      a    = aTmp$obj
                      note = c(note, aTmp$note)
                      #### END OF concurrent

                      # Just before analysing data (because cuncurrent sampling)
                      if (concurrentControlSelect &&
                          (method %in% 'MM') && !activeWindowing) {
                        a_phenlist_concurrent_summary = SummaryStatisticsOriginal(
                          x = a@datasetPL,
                          depVar = depVar,
                          sex = 'Sex',
                          genotype = 'Genotype',
                          label = 'phenlist_and_cuncurrent_data_summary_statistics'
                        )
                        note  =  c(note, a_phenlist_concurrent_summary)
                      }
                      # check the Weight column
                      message0('Checking whether Weight column exists in the raw data ...')
                      if (!CheckIfNameExistInDataFrame(a@datasetPL, 'Weight')) {
                        note$existence_of_weight_column =
                          'Weight column does not exist in the raw data'
                      }
                      # Equation type
                      equationType = ifelse(
                        CheckIfNameExistInDataFrame(a@datasetPL, 'Weight'),
                        getEquation(var = parameter,
                                    equationMap = equationmap),
                        'withoutWeight'
                      )
                      # This is the main engine!
                      note = c(note, list(
                        bodyweight_initially_included_in_model = ifelse(method %in% 'MM', equationType, FALSE)
                      ))
                      if (normalisedPhenlist){
                        a = normalisePhenList(phenlist = a, colnames = c(depVar, 'Weight'))
                      }
                      message0('Fitting the model ...')
                      message0('Method: ', method, '\n\t Equation:', equationType)
                      c.ww0 =	PhenStatWindow(
                        phenlistObject = a,
                        parameter = parameter,
                        minObs = minSampRequired,
                        method = method,
                        depVariable = depVar,
                        equation = equationType,
                        threshold = threshold,
                        pvalThreshold = pvalThreshold,
                        sensitivity = sensitivity,
                        messages = messages,
                        main = paste(
                          unique(n3.5.2$procedure_name)[1],
                          '\n',
                          unique(n3.5.2$parameter_name)[1],
                          '\n',
                          unique(center)[1],
                          unique(colony)[1],
                          unique(zyg)[1],
                          sep = '-',
                          collapse = ','
                        ),
                        seed = seed,
                        check = check,
                        storeplot = storeplot,
                        windowing = activeWindowing,
                        plot = plotWindowing,
                        PicDir =  dirname(outpfile),
                        filename = basename(outpfile),
                        OverwriteExistingFiles = OverwriteExistingFiles,
                        superDebug = superDebug,
                        predFunction = predFunction,
                        residFunction = residFunction,
                        weightORthreshold = weightORthreshold,
                        direction  = direction
                      )
                      note = c(
                        note                               ,
                        c.ww0$note                         ,
                        applied_method = c.ww0$method      ,
                        image_url      = relativePath(
                          path = c.ww0$graphFileName ,
                          reference = wd
                        )
                      )
                      ExtraCols = c('external_sample_id')
                      ####
                      message0('Preparing the output from VectorOutput function ...')
                      c.ww.vec       = VectorOutput0(
                        c.ww0     = c.ww0,
                        ExtraCols    = ExtraCols,
                        activeWindowing = activeWindowing
                      )
                      ##
                      outP    = SuccessfulOutput(
                        args = list(
                          c.ww0 = c.ww0                                  ,
                          depVariable = depVariable                      ,
                          depVar      = depVar                           ,
                          c.ww.vec = c.ww.vec                            ,
                          procedure = procedure                          ,
                          parameter = parameter                          ,
                          center  = center                               ,
                          n3.5   = n3.5                                  ,
                          strain = strain                                ,
                          meta =   meta                                  ,
                          zyg = zyg                                      ,
                          colony = colony                                ,
                          note = note                                    ,
                          PhenListSpecIds = PhenListSpecIds              ,
                          OrgSpecIds = OrgSpecIds                        ,
                          BodyWeightCurvURL = BodyWeightCurvURL          ,
                          GenePageURL = GenePageURL                      ,
                          encode = encode                                ,
                          wd     = wd
                        )
                      )# c(as.list(environment()), ls()))
                      if ((
                        NullOrError(c.ww0$NormalObj)          ||
                        !NullOrError(c.ww0$NormalObj$messages) ||
                        (
                          NullOrError(c.ww0$WindowedObj) &&
                          activeWindowing && (c.ww0$method %in% 'MM')
                        )
                      ) &&
                      debug)
                        save(
                          n3.5,
                          file = paste0(
                            outpfile,
                            RemoveSpecialChars(colony),
                            '_',
                            RandomRegardSeed(),
                            '.Rdata'
                          )
                        )
                      ##
                      StoreRawDataAndWindowingWeights(
                        storeRawData = storeRawData,
                        activeWindowing = activeWindowing,
                        c.ww0 = c.ww0,
                        RawoutputFile = RawoutputFile,
                        orgData = n3.5.2,
                        compressRawData = compressRawData,
                        files  = c(RawoutputFile, ReadMeFile),
                        dir        = dirname (RawoutputFile) ,
                        filename   = basename(RawoutputFile) ,
                        ReadMeTxt  = ReadMeTxt,
                        ReadMeFile = ReadMeFile,
                        methodmap  = methodmap

                      )
                      SucFaiFile = paste(
                        outpfile2,
                        ifelse(
                          !NullOrError(c.ww0$NormalObj) && NullOrError(c.ww0$NormalObj$messages),
                          'Successful.tsv',
                          'Failed_critical_error.tsv'
                        ),
                        sep =
                          '_'
                      )

                      write.table(
                        x = paste(outP,
                                  collapse =   outdelim),
                        file = file.exists0(SucFaiFile, overwrite = OverwriteExistingFiles),
                        append = TRUE,
                        row.names = FALSE,
                        col.names = FALSE,
                        quote = FALSE
                      )
                    } else {
                      message0('Dataset not processed ...')
                      optFail     = NotProcessedOutput(
                        args = list(
                          procedure         = procedure                           ,
                          parameter         = parameter                           ,
                          center            = center                              ,
                          n3.5.2            = n3.5.2                              ,
                          n3.5              = n3.5                                ,
                          strain            = strain                              ,
                          meta              = meta                                ,
                          zyg               = zyg                                 ,
                          colony            = colony                              ,
                          depVariable       = depVariable                         ,
                          #c.ww.vec         = c.ww.vec                            ,
                          note              = note                                ,
                          isException       = isException                         ,
                          minSampRequired   = minSampRequired                     ,
                          depVar            = depVar                              ,
                          GenePageURL       = GenePageURL                         ,
                          BodyWeightCurvURL = BodyWeightCurvURL                   ,
                          OrgSpecIds        = OrgSpecIds                          ,
                          encode 			      = encode                              ,
                          methodmap         = methodmap
                        )
                      )
                      #optFail     = NotProcessedOutput(args = c(as.list(environment()), ls()))
                      NotProcFile = paste(outpfile2,
                                          'NotProcessed.tsv',
                                          sep =
                                            '_')
                      write.table(
                        paste(optFail,
                              collapse =   outdelim),
                        file =   file.exists0(NotProcFile, overwrite = OverwriteExistingFiles),
                        append = TRUE,
                        row.names = FALSE,
                        col.names = FALSE,
                        quote = FALSE
                      )
                      if(simulation)
                        break
                    }
                  }

                  message0(
                    'Finished in ',
                    round(difftime(Sys.time() , Strtime, units = 'sec'), 2),
                    '(s).\n
                  -----------------------------------
                  \n\n '
                  )
                  counter  = counter  + 1
                }
              }
            }
          }
        }
      }
    }
  }
  message0('Closing Connections ...')
  stopCluster(cl)
  registerDoSEQ()
  closeAllConnections()
  stopImplicitCluster()
  message0('Finished.')
  setwd(cwd)
}
