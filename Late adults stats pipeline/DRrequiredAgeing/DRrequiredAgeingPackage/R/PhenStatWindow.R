PhenStatWindow = function (phenlistObject                                ,
                           method = 'MM'                                 ,
                           depVariable = 'Value'                         ,
                           equation = 'withWeight'                       ,
                           parameter = NULL                              ,
                           minObs = 4                                    ,
                           sensitivity = c(1, 1, 1)                        ,
                           pvalThreshold = 0                             ,
                           #### for windowing only
                           windowing = TRUE                              ,
                           seed = 123456,
                           weightFUN = function(x) {
                             # nlme::varComb(nlme::varIdent(form =  ~ 1 |
                             #                                Genotype),
                             #               nlme::varFixed(~ 1 /
                             #                                ModelWeight))
                             nlme::varFixed(~ 1 /
                                              (ModelWeight))
                           },
                           check = 2                                     ,
                           messages = FALSE                              ,
                           main = ''                                     ,
                           threshold = 10 ^ -18                          ,
                           plot   = TRUE                                 ,
                           storeplot = TRUE                              ,
                           PicDir = NULL                                 ,
                           OverwriteExistingFiles  = FALSE               ,
                           filename = RandomRegardSeed(1)                ,
                           superDebug = FALSE                            ,
                           predFunction = predFunction                   ,
                           residFunction = residFunction                 ,
                           weightORthreshold = 'weight'                  ,
                           maxPeaks = 15                                 ,
                           direction = direction                         ,
                           ########
                           outlierDetection = FALSE                      ,
                           #######
                           FERROptimise                                  ,
                           MMOptimise                                    ,
                           FERRrep                                       ,
                           ...)
{
  requireNamespace('PhenStat')
  requireNamespace('SmoothWin')
  requireNamespace('OpenStats')
  requireNamespace('nlme')
  set.seed(seed)
  ### Outliwe detection
  if (outlierDetection && method %in% 'MM') {
    message0('Outlier detection in progress ...')
    phenlistObject = PhenListOutlierDetection(phenlistObject     ,
                                              plot   = !storeplot,
                                              active = outlierDetection)
  } else{
    phenlistObject@datasetPL$outlierWeight = 1
  }
  ###########################
  # Do not remove lines below (necessary for windowing)
  phenlistObject@datasetPL = SimplesortDataset(x    = phenlistObject@datasetPL,
                                               cols = c('external_sample_id', 'discrete_point'))
  if (CheckIfNameExistInDataFrame(phenlistObject@datasetPL, 'Batch')) {
    message0('Sorting the dataset on Batch')
    phenlistObject@datasetPL = sortDataset(x        = phenlistObject@datasetPL,
                                           BatchCol = 'Batch')
  }
  ###########################
  if (method %in% 'MM' &&
      unique(phenlistObject@datasetPL$observation_type) %in% 'time_series') {
    RandEffTerm = if (CheckIfNameExistInDataFrame(phenlistObject@datasetPL, 'LifeStage'))
      as.formula('~ 1 | external_sample_id')
    else
      as.formula('~ 1 | external_sample_id/Batch')
    CorrEffect  = nlme::corAR1(form = RandEffTerm) #corSymm()
    # windowing does not apply to time series
    message0('Time series detected! Windowing does not apply to the time series ...')
    windowing = FALSE
  } else{
    RandEffTerm = if (CheckIfNameExistInDataFrame(phenlistObject@datasetPL, 'LifeStage'))
      as.formula('~ 1 | Batch')
    else
      as.formula('~ 1 | Batch')
    CorrEffect  = NULL
  }
  ###########################
  # Run normal (not windowed) models
  ###########################
  note = windowingNote = graphFileName = object0 = NULL
  message0(method, ' in progress ....')
  object0 =   OpenStats::OpenStatsAnalysis(
    OpenStatsListObject   = phenlistObject,
    method                = method,
    MM_random             = RandEffTerm,
    correlation           = CorrEffect,
    MM_optimise           = MMOptimise,
    FERR_FullComparisions = FERROptimise,
    FERR_rep              = FERRrep,
    MM_BodyWeightIncluded = ifelse(equation %in% 'withWeight', TRUE, FALSE),
    debug                 = TRUE
  )
  note$'Normal analysis'$'Step 1-1 messages' = object0$messages
  # If not possible use MM (only for ABR)
  if (is.ABR  (x            = parameter) &&
      !is.null (object0$messages))
  {
    method = 'MM'
    message0('Running the MM (only ABR) ... ')
    object0 =   OpenStats::OpenStatsAnalysis(
      OpenStatsListObject = phenlistObject,
      method              = method,
      MM_BodyWeightIncluded = ifelse(equation %in% 'withWeight', TRUE, FALSE),
      MM_random   = RandEffTerm,
      correlation = CorrEffect,
      debug       = TRUE
    )
    note$'Normal analysis'$'Step 1-2 Mixed Model messages' = object0$messages
  }
  ############## If not possible add jiter
  if (is.ABR(x = parameter) &&
      !is.null (object0$messages)) {
    message0('Running the MM+Jitter (only ABR) ... ')
    method = 'MM'
    phenlistObject@datasetPL[, depVariable] = jitter(phenlistObject@datasetPL[, depVariable], 0.1)
    object0 =   OpenStats::OpenStatsAnalysis(
      OpenStatsListObject   = phenlistObject,
      method                = method,
      MM_BodyWeightIncluded = ifelse(equation %in% 'withWeight', TRUE, FALSE),
      MM_random   = RandEffTerm,
      correlation = CorrEffect,
      debug       = TRUE
    )
    note$'Normal analysis'$'Step 1-3 Mixed Model with added Jitter messages' = object0$messages
  }

  if(!is.null(object0$messages))
    object0 = NULL

  ###########################
  #### This is for windowing only
  ###########################
  if (windowing && method %in% c('MM') &&
      is.numeric(phenlistObject@datasetPL[, depVariable]))
  {
    # Full model
    message0('Running the full model before applying windowing ... ')
    # objectNorm = OpenStatsAnalysis(
    #   OpenStatsListObject = phenlistObject,
    #   method              = method,
    #   MM_random   = RandEffTerm,
    #   correlation = CorrEffect,
    #   MM_BodyWeightIncluded = ifelse(equation %in% 'withWeight', TRUE, FALSE),
    #   MM_optimise = c(0, 1, 1),
    #   debug       = TRUE
    # )
    ############################### Be careful about this section!
    objectNorm = object0
    objectNorm$output$Final.Model = objectNorm$output$Initial.Model
    ###############################
    windowingNote$'Windowing analysis'$'Fully loaded model' = objectNorm$messages
    #######################################################
    message0('Start windowing ... ')
    #######################################################
    if (!is.null(objectNorm$output$Initial.Model)) {
      obj = objectNorm$output$Initial.Model
      phenlistObject@datasetPL = nlme::getData(obj) # !important
      ###
      tt = Date2Integer(phenlistObject@datasetPL$Batch)
      mm = mm.bck = which(!phenlistObject@datasetPL$Genotype %in% phenlistObject@refGenotype)
      message0('The number of modes: ', length(unique(mm)))
      # Windowing & DOE
      if (length(unique(tt[mm])) > maxPeaks) {
        message0(
          'More than ',
          maxPeaks,
          ' modes. A random sample of ',
          maxPeaks,
          ' will be used. total modes: ',
          length(unique(mm))
        )
        sa            = sample(unique(tt[mm]), maxPeaks)
        mm            = mm[which(tt[mm] %in% sa)]
      }
      ####
      message0('Windowing algorithm in progress ...')
      r = SmoothWin(
        object = obj                            ,
        data = phenlistObject@datasetPL         ,
        t = tt                                  ,
        m = mm                                  ,
        weightFUN = weightFUN                   ,
        messages = messages                     ,
        check = check                           ,
        seed = seed                             ,
        threshold = threshold                   ,
        simple.output = TRUE                    ,
        sensitivity = sensitivity               ,
        pvalThreshold = pvalThreshold           ,
        debug = superDebug                      ,
        residFun = residFunction                ,
        predictFun = predFunction               ,
        weightORthreshold = weightORthreshold   ,
        direction = direction                   ,
        min.obs = function(ignore.me.in.default) {
          message0('Total number of sex: ',
                   noSexes0(phenlistObject))
          lutm = length(unique(tt[mm]))
          r = ifelse(lutm > 1,
                     noSexes0(phenlistObject) * 35,
                     max(pi * sqrt(length(tt)), 35))
          r = max(r * lutm, length(mm), na.rm = TRUE)
          r = min(r       , length(tt), na.rm = TRUE)
          message0('min.obs =  ', r)
          return(r)
        },
        zeroCompensation = threshold*10^-3,
        externalWeight   = phenlistObject@datasetPL$outlierWeight
      )
      ##############################
      phenlistObject@datasetPL$AllModelWeights = we = we2 =  r$finalModel$FullWeight
      if (is.null(we) ||
          length(unique(we)) < 2 ||
          var(we, na.rm = TRUE) < threshold) {
        we2  = NULL
        windowingNote$'Windowing extra' = 'There is no variation in the weights then the standard model is applied.'
      } else{
        ####################
        vMutants  = MeanVarOverTime(
          mm = (1:length(tt))[mm.bck],
          tt = tt,
          data = phenlistObject@datasetPL[, depVariable]
        )
        VControls = MeanVarOverTime(
          mm = (1:length(tt))[-mm.bck],
          tt = tt,
          data = phenlistObject@datasetPL[, depVariable]
        )
        message0('Disabled but: Mutant sd = ',
                 vMutants,
                 ', Control sd = ',
                 VControls)
        #we2[mm]  = we2[mm] * vMutants
        #we2[-mm] = we2[-mm] * VControls
      }

      MM_optimise2 = MMOptimise
      MM_optimise2[2] = 0

      if (!is.null(we2)) {
        message0('Reclaculating an optimised model using the windowing weights ...')

        # in this step we remove off threshold weights (to solve the computation error of small values)
       phenlistObject2 = removeOffweightsFromPLobject(phenlistObject,threshold)

        objectf = OpenStatsAnalysis(
          OpenStatsListObject = phenlistObject2,
          method              = method,
          MM_random   = RandEffTerm   ,
          correlation = CorrEffect    ,
          MM_BodyWeightIncluded = ifelse(equation %in% 'withWeight', TRUE, FALSE),
          MM_weight   = nlme::varComb(nlme::varFixed(~ 1 / AllModelWeights)),
          debug       = TRUE,
          MM_optimise = MM_optimise2
        )
      } else{
        message0('Ops no, the soft windowing failed! Just a normal full model will be estimated ...')
        objectf = OpenStatsAnalysis(
          OpenStatsListObject = phenlistObject,
          method              = method,
          MM_random   = RandEffTerm   ,
          correlation = CorrEffect    ,
          MM_BodyWeightIncluded = ifelse(equation %in% 'withWeight', TRUE, FALSE),
          debug       = TRUE,
          MM_optimise = c(0, 1, 1, 1, 1, 1)
        )
      }
      windowingNote$'Windowing analysis'$'Final model' = objectf$messages
      # Full model windowing
      message0('Fitting the windowing weights into the full PhenStat model ...')
      # objectfulw = OpenStatsAnalysis(
      #   OpenStatsListObject = phenlistObject,
      #   method              = method,
      #   MM_random   = RandEffTerm,
      #   correlation = CorrEffect,
      #   MM_BodyWeightIncluded = ifelse(equation %in% 'withWeight', TRUE, FALSE),
      #   MM_weight = nlme::varComb(
      #     nlme:: varFixed( ~ 1 / AllModelWeights)
      #   ),
      #   MM_optimise = c(0,1,1),
      #   debug       = TRUE
      # )
      ############################### Be careful about this section!
      objectfulw = objectf
      objectfulw$output$Final.Model = objectfulw$output$Initial.Model
      ###############################
      windowingNote$'Windowing analysis'$'Full model windowed result' = objectfulw$messages
      # Plotting
      args = list(
        # for output
        r = r           ,
        tt = tt         ,
        mm = mm.bck     ,
        we = we         ,
        threshold         = threshold         ,
        maxPeaks          = maxPeaks          ,
        ## for plot
        plot              = plot              ,
        PicDir            = PicDir            ,
        filename          = filename          ,
        storeplot         = storeplot         ,
        phenlistObject    = phenlistObject    ,
        check             = check             ,
        main              = main              ,
        depVariable       = depVariable       ,
        threshold         = threshold         ,
        external_sample_ids = nlme::getData(obj)[,'external_sample_id'],
        observation_id      = nlme::getData(obj)[,'observation_id']
      )
      #args    = c(as.list(environment()), list())
      windowingNote$'Window parameters' = WindowingDetails(args)
      graphFileName  = PlotWindowingResult(args = args, overwrite = OverwriteExistingFiles)
      ObjectsThatMustBeRemovedInEachIteration(c('args', 'objectNorm','objectfulw'))
    } else{
      message0('An error in forming the full model (windowing only) ... ')
      objectf = objectNorm = objectfulw = NULL# object0
      r       = we         = NULL
    }
  } else{
    objectf   = objectNorm = objectfulw = NULL# object0
    r         = we         = NULL
  }

  outputList = list(
    InputObject     = phenlistObject,
    # Full model (not windowed)
    FullObj         = objectNorm,
    # Optimised model (not windowed)
    NormalObj       = object0   ,
    # Optimized windowed model
    WindowedObj     = objectf   ,
    # Full model windowed
    FullWindowedObj = objectfulw,
    WinDetails      = r             ,
    weight          = we            ,
    method          = method        ,
    graphFileName   = ifelse(is.null(graphFileName), 'NA', graphFileName),
    note = c(note, windowingNote)
  )
  if (superDebug) {
    message0('SuperDebug is activated. Writting the Rdata file ... ')
    SupDebFile = file.exists0(file.path(PicDir,
                                        'superDebugAnalysis.Rdata'))
    message0('Superdebug file: \n\t\t ~> ', SupDebFile)
    agg = c(as.list(environment()), list())
    save(agg,
         file = SupDebFile)
  }
  message0('Running GC() ....')
  gc()
  message0('End of the core process ....')
  return(
    outputList
  )
}
