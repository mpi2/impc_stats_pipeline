# Startup message
.onAttach <- function(lib, pkg) {
  packageStartupMessage(
    paste0(
      '\n >===============================================================================<',
      '\n "DRrequired" is developed by International Mouse Phenotyping Consortium (IMPC) ',
      '\n More details https://www.mousephenotype.org/                                            ',
      '\n Contact us hamedhm@ebi.ac.uk                                                  ',
      '\n >===============================================================================<'
    ),
    domain = NULL,
    appendLF = TRUE
  )
}

complete.cases0 = function(x, ...) {
  if (!is.null(x) && length(x) > 0) {
    r = complete.cases(x)
    if (is.null(r))
      r = NULL
  } else{
    message0('Null input in the "complete.case". Null retured!')
    r = NULL
  }
  return(r)
}

# Objects that must be removed at each iterations
ObjectsThatMustBeRemovedInEachIteration = function(x = NULL, ...) {
  if (is.null(x)) {
    rm0(
      c(
        'note',
        'minSampRequired',
        'n3.4',
        'n3.5',
        'lvls',
        'counts',
        'n3.5.1',
        'note',
        'SubSubDirOrdFileName',
        'outDir',
        'outpfile',
        'AllLevels',
        'GenePageURL',
        'BodyWeightCurvURL',
        'depVariable',
        'depVar',
        'n3.5_summary',
        'n3.5.1_F_list',
        'n3.5.1_v_list',
        'isException',
        'ExtraCols',
        'a',
        'PhenListSpecIds',
        'OrgSpecIds',
        'a_summary_before_concurrent',
        'a_phenlist_concurrent_summary',
        'method',
        'args',
        'aTmp',
        'c.ww0',
        'c.ww.vec',
        'outP',
        'optFail'
      )
    )
  } else{
    rm0(x)
  }
}

getNotEmptyValue = function(x, head = Inf) {
  ux    = unique(x)
  neux  = which(ux != '' | !is.null(ux))
  if (length(neux) > 0) {
    r = head(ux [neux], head)
  } else{
    r = 'NotExist'
  }
  return(as.character(r))
}

# Only works for the data from Solr
# Automatic select the corresponding column for the type of data (only for category and data_point)
getResponseColumn = function(x, activate = TRUE) {
  if (activate) {
    if (length(x) > 0) {
      lbl      = unique(x)
      column   = paste(lbl, collapse = '-', sep = '-')
      accepted = TRUE

      if (length(lbl) > 1)
        accepted = FALSE

      if (lbl == 'categorical')
        column = 'category'
      else if (lbl == 'unidimensional')
        column = 'data_point'
      else
        accepted = FALSE
    } else{
      accepted = FALSE
      column = lbl = 'Unknown'
    }
  } else{
    lbl    = NA
    column = lbl = NA
    accepted = FALSE
  }
  return(list(
    column = column,
    accepted = accepted,
    lbl = lbl
  ))
}


IsInBlackListCategories = function(x, len = 1, blackList = NULL) {
  note = NULL
  if (!is.null(blackList) && !is.numeric(x) && nlevels(x) == len) {
    r    = levels(x) %in% blackList
    note = 'The categorical variable has only one level that is found in the skipt list'
  } else{
    r    = FALSE
  }
  return(list(result = r, note = note))
}



#
MergeLevels = function(x, listOfLevelMaps) {
  note = NULL
  if (!is.null(x) && !is.numeric(x) &&  length(x) > 0 &&
      nlevels(as.factor(x)) > 0) {
    x  = as.factor(x)
    nl = nlevels(x)
    l  = levels(x)
    for (i in 1:nl) {
      if (l[i] %in% names(listOfLevelMaps)) {
        newCat = as.character(listOfLevelMaps[l[i] == names(listOfLevelMaps)])
        levels(x)[which(levels(x) == l[i])] = ifelse(newCat != 'NA', newCat, NA)
        note = c(note, paste0(as.character(l[i]), '-->', newCat))
      }
    }
  }
  return(list(x = x, note = if (is.null(note)) {
    note
  } else{
    paste(note, collapse = ', ')
  }))
}


# Gene page URL
GenePageURL = function(obj) {
  if (any(dim(obj) == 0) || is.null(obj)) {
    return ('Empty dataset has no URL!')
  }
  activation  = c(
    'gene_accession_id',
    'allele_accession_id',
    'zygosity',
    'parameter_stable_id',
    'pipeline_stable_id',
    'phenotyping_center',
    'metadata_group'
  ) %in% names(obj)
  ##########
  if (prod(activation) > 0) {
    r = paste(
      'https://www.mousephenotype.org/data/charts?'  ,
      'accession='                                   ,
      getNotEmptyValue(obj$gene_accession_id, 1)     ,
      '&allele_accession_id='                        ,
      getNotEmptyValue(obj$allele_accession_id, 1)   ,
      '&zygosity='                                   ,
      getNotEmptyValue(obj$zygosity, 1)        ,
      '&parameter_stable_id='                        ,
      getNotEmptyValue(obj$parameter_stable_id, 1)   ,
      '&pipeline_stable_id='                         ,
      getNotEmptyValue(obj$pipeline_stable_id, 1)    ,
      '&phenotyping_center='                         ,
      getNotEmptyValue(obj$phenotyping_center, 1)    ,
      '&metadata_group='                             ,
      getNotEmptyValue(obj$metadata_group, 1)        ,
      sep = ''                                       ,
      collapse = ''
    )
    r = URLencode(r, repeated = TRUE)
  } else{
    r = NA
  }
  return(r)
}

# BWT URL
BodyWeightCurvURL = function(obj) {
  if (any(dim(obj) == 0) || is.null(obj)) {
    return ('Empty dataset has no BWT URL!')
  }
  # allele_accession_id to make sure that we are safe for external application of the R package
  activation = c('allele_accession_id', 'gene_accession_id') %in% names(obj)
  if (prod(activation) > 0) {
    r = paste(
      'https://www.mousephenotype.org/data/charts?'    ,
      'accession='                                     ,
      getNotEmptyValue(obj$gene_accession_id, 1)       ,
      '&parameter_stable_id='                          ,
      'IMPC_BWT_008_001'                               ,
      '&chart_type=TIME_SERIES_LINE'                   ,
      sep = ''                                         ,
      collapse = ''
    )
    r = URLencode(r, repeated = TRUE)
  } else{
    r = NA
  }
  return(r)
}

#### create dir
dir.create0 = function(x, ...) {
  if (!dir.exists(x))
    dir.create(x, ...)
  return(x)
}
###
LogicalToNote = function(x) {
  if (x) {
    r = 'Enabled'
  } else{
    r = 'Disabled'
  }
  return(r)
}
#### Concurrent sample selection
concurrentContSelect = function(activate = TRUE,
                                PhenListObj,
                                minSmp  = 6,
                                depVar  = 'data_point',
                                control = 'control'   ,
                                sexColName      = 'Sex'     ,
                                GenotypeColName = 'Genotype',
                                BatchColName    = 'Batch'   ,
                                minSampRequired,

                                plot = FALSE) {
  note = NULL
  PhenListObjBck = PhenListObj
  # Old keys
  # cun_con_sel | cun_con_selection | concurrent_control_selection
  # Part 1 messages
  if (is.numeric(PhenListObj@datasetPL[, depVar])) {
    m  =  (PhenListObj@datasetPL[PhenListObj@datasetPL[, GenotypeColName] != control, ])
    m1 =  (PhenListObj@datasetPL[PhenListObj@datasetPL[, BatchColName] %in% m[, BatchColName], ])
    if (length(unique(m[, BatchColName])) == 1) {
      tb = table(m1[, sexColName], m1[, GenotypeColName])
      if (control %in% colnames(tb)) {
        if (min0(tb[, which(colnames(tb) == control)]) >= minSmp) {
          if (activate) {
            message0 ('Concurrent control selection in progress ...')
            PhenListObj@datasetPL = m1
            note$concurrent_control_selection = paste0(LogicalToNote(activate),
                                                       ". Single Batch in use")
          } else{
            note$concurrent_control_selection = paste0(LogicalToNote(activate),
                                                       '. Single Batch possible')
          }
        } else{
          note$concurrent_control_selection = paste0(
            LogicalToNote(activate),
            '. There are not enough controls in the dataset at the mutant dates  (min = ',
            minSmp,
            ' but there are ',
            min0(tb[, which(colnames(tb) == control)]),
            ')'
          )
        }
      } else{
        note$concurrent_control_selection = paste0(
          LogicalToNote(activate),
          ". The is no control in the dataset at the same day as mutants DOE"
        )
      }
    } else{
      note$concurrent_control_selection = paste0(
        LogicalToNote(activate),
        ". Concurrent control selection not possible. Mutants are scattered on different dates"
      )
    }
  } else {
    tb  = table(PhenListObj@datasetPL[, BatchColName])
    if (length(tb) > 1) {
      note$concurrent_control_selection = paste0(LogicalToNote(activate),
                                                 '. Multi batch found in data (categorical data)')
    } else if (length(tb) == 1) {
      note$concurrent_control_selection = paste0(LogicalToNote(activate),
                                                 '. Original mutants in single batch (categorical data)')
    } else{
      note$concurrent_control_selection = paste0(LogicalToNote(activate),
                                                 '. But batch not clear. Please check data')
    }
  }
  # Part 2 extra checks
  if (activate && is.numeric(PhenListObj@datasetPL[, depVar])) {
    ### Check whether there is any variation in data
	### No worries it only applies to MM framework
    checkNoZeroVar = RemoveZerovarCategories(
      x = PhenListObj@datasetPL,
      depVar = depVar,
      sex = sexColName,
      genotype = GenotypeColName
    )
    #PhenListObj@datasetPL = checkNoZeroVar$x
    PhenListObj = PhenStat::PhenList(
      dataset = checkNoZeroVar$x,
      testGenotype = PhenListObjBck@testGenotype,
      refGenotype = PhenListObjBck@refGenotype
    )
    note =  c(note, checkNoZeroVar$note)


    conditions =  (
      !is.null(PhenListObj@datasetPL) &&
        # data.frame is not zero
        min0(dim(PhenListObj@datasetPL)) > 0 &&
        # is it  really exist!
        length(unique(PhenListObj@datasetPL[, GenotypeColName])) > 1 &&
        # include mut and cont
        min0(table(PhenListObj@datasetPL[, GenotypeColName])) >= minSampRequired &&
        # include at least 4/2 of each genotype
        checkGenInCol(
          PhenListObj@datasetPL,
          sex = sexColName,
          genotype = GenotypeColName
        ) &&
        # each sex and genotype
        # there must be variation in data
        NonZeroVariation(PhenListObj@datasetPL[, depVar])
    )

    if (!conditions) {
      PhenListObj = PhenListObjBck
      note$Concurrent_control_selection_error = paste0(
        LogicalToNote(activate),
        'Cuncurrent control selection would cause failure in model then DISABLED'
      )
    }
    ####################
  }
  if (plot)
    plot(PhenListObj, depVariable = depVar)

  return (list(obj = PhenListObj, note = note))
}



# Variation in data
NonZeroVariation = function(x) {
  x = na.omit(x)
  result = ifelse(length(x) > 0 &&
                    length(unique(x)) > 1, TRUE, FALSE)
  return(result)
}

# Is this an ABR screen
is.ABR = function(x) {
  length(grep("_ABR", x, fixed = TRUE)) > 0
}

# Is this a 3I screen
is.3i = function(x) {
  length(
    grep(
      "MGP_BCI|MGP_PBI|MGP_ANA|MGP_CTL|MGP_EEI|MGP_BMI|MGP_MLN|MGP_IMM|MGP_BHP|MGP_MIC",
      x,
      fixed = FALSE
    )
  ) > 0

}

# Exception map
getException =  function(parameter,
                         exceptionMap = NULL)
{
  matched001 = sapply(names(exceptionMap), grepl, parameter)
  if (sum(matched001) > 0) {
    result = FALSE
  } else{
    result = TRUE
  }
  return(result)
}

# Method map
getMethodi =  function(var,
                       type = 'numeric',
                       methodMap = NULL)
{
  matched001 = sapply(names(methodMap), grepl, var)
  if (sum(matched001) > 0) {
    if (sum(matched001) > 1) {
      methodMapReduced = methodMap[matched001]
      mappedPatternLengths = nchar(names(methodMapReduced))
      method =
        unlist(methodMapReduced[which(mappedPatternLengths ==  max(mappedPatternLengths))])
    } else {
      method = unlist(methodMap[matched001])
    }
  } else{
    if (type == 'numeric') {
      method = 'MM'
    } else{
      method = 'FE'
    }
  }
  return(method)
}

# Equation map
getEquation =  function(var,
                        equationMap = NULL)
{
  matched001 = sapply(names(equationMap), grepl, var)
  if (sum(matched001) > 0) {
    if (sum(matched001) > 1) {
      equationMapReduced = equationMap[matched001]
      mappedPatternLengths = nchar(names(equationMapReduced))
      equation =
        unlist(equationMapReduced[which(mappedPatternLengths ==  max(mappedPatternLengths))])
    } else {
      equation = unlist(equationMap[matched001])
    }
  } else{
    equation = 'withWeight'
  }
  return(equation)
}


# Read config files
readConf = function(file, path = NULL, ...) {
  if (is.null(path))
    path = system.file("extdata", package = "DRrequired")

  r   = read.dcf(file.path(path, file), all = TRUE,  ...)
  return(r)
}
# Read config files
readFile = function(file,
                    path = NULL,
                    FUN = function(x, ...) {
                      as.vector(read.table(file = x)[, 1])
                    },
                    ...) {
  if (is.null(path))
    path = system.file("extdata", package = "DRrequired")

  r   = FUN(file.path(path, file),  ...)
  return(r)
}

# Replace blank elements of a vector with user define charachter
RepBlank = function(x, match = '', replace = 'control') {
  x[x %in% match] = replace
  return(x)
}

# To prevent:
#                   female male
# control           4       0
# experimental      0       5
checkGenInCol = function(x, sex = 'sex', genotype = 'biological_sample_group') {
  x  = droplevels0(x)
  r = FALSE
  # not null
  if (min0(dim(x)) > 0) {
    tb = table(x[, genotype], x[, sex])
    if (min0(dim(tb)) > 0) {
      if (dim(tb)[2] > 1)
        r = (prod(tb[, 1]) + prod(tb[, 2])) > 0
      else
        r = (prod(tb[, 1])) > 0
    }
  }
  return(r)
}


##' Catch and save both errors and warnings, and in the case of
##' a warning, also keep the computed result.
##'
##' title tryCatch both warnings (with value) and errors
##' param expr an R expression to evaluate
##' return a list with 'value' and 'warning', where
##'   'value' may be an error caught.
##' author Martin Maechler;
##' Copyright (C) 2010-2012  The R Core Team
tryCatch.W.E0 <- function(expr)
{
  W = c('No warning found')
  E = c('No error found')
  w.handler <- function(w) {
    W <<- w
  }
  e.handler <- function(e) {
    E <<- e
    value = NULL
  }
  list(
    value = tryCatch(
      expr,
      error = e.handler,
      warning = w.handler
    ),
    warning = W,
    error = E,
    both = NULL
  )
}




tryCatch.W.E <- function(expr,
                         name,
                         sep = '_',
                         col = '_',
                         ...) {
  r = tryCatch.W.E0(expr = expr, ...)
  if (NullOrError(r$value)) {
    if ('call' %in% names(r$error))
      r$error$call = NULL
    if ('call' %in% names(r$warning))
      r$warning$call = NULL

    ### Assembling warnings and errors in one object
    name1 =  paste(name,
                   'functional_warnings',
                   sep = sep,
                   collapse = col)
    v1 =   if (all(r$warning == '')) {
      'No functional warning'
    } else{
      unlist(r$warning)
    }
    name2 = paste(name,
                  'functional_errors',
                  sep = sep,
                  collapse = col)
    v2 =  if (all(r$error == '')) {
      'No functional warning'
    } else{
      unlist(r$error)
    }

    r$both = list(v1, v2)
    names(r$both) = c(name1, name2)
  }
  return(r)
}

### Jeremy original script to read the errors/messages from the function
capture.stderr <- function (..., file = NULL, append = FALSE) {
  args <- substitute(list(...))[-1L]
  rval <- NULL
  closeit <- TRUE
  if (is.null(file)) {
    file <-
      textConnection(paste("rval",
                           #runif(1, min = 1, max = 10 ^ 5),
                           collapse = '_',
                           sep = '_'), "w", local = TRUE)
  } else if (is.character(file)) {
    file <- file(file,
                 ifelse (append, 'a', 'w'))
  } else if (inherits(file, "connection")) {
    if (!isOpen(file)) {
      open(file,
           ifelse (append, 'a', 'w'))
    } else{
      closeit <- FALSE
    }
  } else{
    stop("'file' must be NULL, a character string or a connection")
  }

  sink(file, type = "message")
  on.exit({
    sink(type = "message")
    if (closeit)
      close(file)
  })
  pf <- parent.frame()
  evalVis <- function(expr)
    withVisible(eval(expr, pf))

  result01 <- c()
  for (i in seq_along(args)) {
    expr <- args[[i]]
    tmp <- switch(
      mode(expr),
      expression = lapply(expr, evalVis),
      call = ,
      name = list(evalVis(expr)),
      stop("bad argument")
    )
    result01 <- c(result01, tmp)
  }
  on.exit()
  sink(type = "message")
  if (closeit)
    close(file)
  if (is.null(rval))
    invisible(list("output" = NULL, "result" = result01))
  else
    list("output" = trimws(na.omit(rval)), "result" = result01)
}


### check for null or errorneous object
NullOrError = function(obj) {
  if (is.null(obj) ||
      class(obj) %in% c("simpleWarning", "warning", "condition")) {
    return(TRUE)
  } else{
    return(FALSE)
  }
}


# Generate Batch file only for EBI
BatchGenerator = function(file                       ,
                          dir = NULL                 ,
                          procedure = NULL           ,
                          parameter = NULL           ,
                          center = NULL              ,
                          cpu = 2                    ,
                          memory = 5000              ,
                          extraBatchParameters = NULL) {
  dirOut = file.path(dir, 'ClusterOut')
  dirErr = file.path(dir, 'ClusterErr')
  dir.create0(dirOut)
  dir.create0(dirErr)

  oname = file.path(dirOut, paste(procedure, '_', parameter, sep = ''))
  ename = file.path(dirErr, paste(procedure, '_', parameter, sep = ''))

  ro = paste(' -o ', paste0('"', oname, '.ClusterOut', '"'), sep = '')
  re = paste(' -e ', paste0('"', ename, '.ClusterErr', '"'), sep = '')
  rf = paste(
    'bsub '               ,
    extraBatchParameters  ,
    ' -n '                ,
    cpu                   ,
    ' -M '                ,
    memory                ,
    ' '                   ,
    '-R "rusage[mem='        ,
    memory                ,
    ']"'                  ,
    ' '                   ,
    re                    ,
    ' '                   ,
    ro                    ,
    ' Rscript function.R ',
    paste('"', file, '"', sep = ''),
    sep = ''
  )
  return(rf)
}


# remove special characters from a string
RemoveSpecialChars = function(x,
                              what = '[^0-9A-Za-z]',
                              replaceBy = '_',
                              message = FALSE) {
  what = gsub(
    pattern = ']',
    x = what,
    replacement = paste0(replaceBy, ']')
  )
  if (message)
    message0('pattern: ',
             what ,
             '; replaced by ',
             replaceBy)
  r = gsub(what, replaceBy , x , ignore.case = TRUE)
  ###
  if (any(nchar(r) < 1)) {
    RN = RandomRegardSeed(1, stringOutput = TRUE, round = 8)
    r[nchar(r) < 1] = paste('no_name',
                            RN,
                            sep = '_',
                            collapse = '_')
  }
  return(r)
}


# Create sub directory from a url if does not exist
CreateSubDirIfNotExist = function(x) {
  base = dirname(x)
  dir.create0(base)
  return(x)
}

# Unique and not NULL
UniqueAndNNull = function(x,
                          removeSpecials = TRUE,
                          removeNewLine  = TRUE,
                          newlineRepChar = ' ' ,
                          collapse = '~') {
  if (length(na.omit(x)) > 0)
    x = na.omit(x)
  if (length(x[!is.null(x)]) > 0)
    x = x[!is.null(x)]
  if (length(x[!x %in% '']) > 0)
    x = x[!x %in% '']
  if (length(x) > 0) {
    x = as.character(unique(x))
    x = paste0(x, collapse = collapse)
    if (removeNewLine)
      x = gsub(pattern = '\n',
               x = x,
               replacement = newlineRepChar)
    if (removeSpecials)
      x = RemoveSpecialChars(x = x, replaceBy = ' ')
  } else{
    x = 'NotSpecified'
  }
  return(x)
}

# remove useless characters
removeNAPipe = function(x,
                        pattern = 'NA',
                        replacement = '') {
  x = na.omit(x)
  if (length(x) > 0)
    x = gsub(
      pattern = pattern,
      replacement = replacement,
      x = x,
      fixed = FALSE
    )
  return(x)
}


# as.numeric with quote!
as.numeric0 = function(x) {
  r = suppressWarnings(as.numeric(x))
  return(r)
}

# is.numeric with quote!
is.numeric0 = function(x) {
  r0  = is.na(as.numeric0(x))
  r1  = is.na(x) | (x %in% c('na', 'NA', 'TRUE', 'FALSE'))
  r   = r0 | r1
  return(!r)
}

## solve "}"]  and ["
RemoveTwoSpecialFromOutput  =  function(x) {
  ind = grepl(pattern = '{',
              x = x,
              fixed = TRUE)
  for (i in 1:length(x)) {
    if (!ind[i]) {
      x[i] = paste(
        ifelse(length(x[i]) > 1, '[', ''),
        ifelse(is.numeric0(x[i]), as.numeric0(x[i]), paste0('"', x[i], '"')),
        ifelse(length(x[i]) > 1, ']', ''),
        sep = ''
      )
    } else{
      x[i] = paste('{', x[i], '}', sep = '')
    }
  }
  return(x)
}

### removeObject if exist
rm0 = function(x, silent = FALSE) {
  x      = unique(x)
  obList = ls(all.names = TRUE, pos = ".GlobalEnv")
  InObj  = obList %in% x
  if (sum(InObj) > 0) {
    rm(list = as.character(obList[InObj]), pos = ".GlobalEnv")
    if (!silent)
      message0(
        'Cleaning the environment. Removed objects from memory:\n\t ~~~~>>>> ',
        paste(obList[InObj], collapse = ', ')
      )
  } else{
    if (!silent)
      message0('Cleaning the environment. No object removed from the memory')
  }
}



# Generate random numbers even when seed is set (max 10^6 numbers)
RandomRegardSeed = function(n = 1,
                            max = 10 ^ 9,
                            decimal = TRUE,
                            round = 5,
                            stringOutput = FALSE,
                            what = '[^0-9A-Za-z]',
                            replaceBy = '_') {
  r = runif(n, 1, as.numeric(Sys.time()) * 10000) %% max
  if (decimal)
    r = r %% 1
  r = round(r, digits = round)
  if (stringOutput) {
    r = RemoveSpecialChars(as.character(r), what = what, replaceBy = replaceBy)
  }
  return(r)
}


# remove zero count categories (filter on sex and genotype)
RemoveZeroFrequencyCategories = function(x,
                                         minSampRequired,
                                         depVar = 'data_point',
                                         drop = TRUE,
                                         totalLevels = 4,
                                         sexCol = 'sex',
                                         genotypeCol = 'biological_sample_group',
                                         sep = '_',
                                         activated = TRUE) {
  note   = NULL
  if (all(dim(x) > 0) && !is.null(complete.cases0(x[, depVar]))) {
    x = x[complete.cases0(x[, depVar]), ,drop=FALSE]
    if (is.numeric(x[, depVar])) {
      lvls   = interaction(x[, sexCol], x[, genotypeCol], sep = sep, drop = drop)
    } else{
      lvls   = interaction(x[, sexCol], x[, genotypeCol], sep = sep, drop = drop)
      #lvls   = interaction(x[, sexCol], x[, genotypeCol], x[, depVar], sep = sep, drop = drop)
    }
    lvls   = lvls[!is.na(x[, depVar])]
    counts = table(lvls)
    newx   = (x[lvls %in% names(counts)[counts >= minSampRequired], ])
    #####
    if ((length(names(counts[counts >= minSampRequired])) != totalLevels ||
         !identical(droplevels0(newx), droplevels0(x))) &&
        length(lvls) > 0) {
      note$sex_genotype_data_included_in_analysis = c(
        paste(
          'Not all Sex*Genotype interactions exist in data. Threshold: ',
          paste0(minSampRequired),
          ' sample per group. INCLUDED categories: ',
          paste(names(counts)[counts >= minSampRequired],
                sep = ',',
                collapse = ', '),
          '',
          sep = '',
          collapse = ''
        )
      )
    }
  } else{
    newx = NULL
  }
  if (any(dim(newx) == 0) || any(dim(x) == 0))
    newx = NULL

  return(list(x = newx, note = note))
}

# Droplevels does not work with null/empty objects
droplevels0 = function(x, ...) {
  if (!is.null(x) && all(dim(x) > 0)) {
    r = droplevels(x = x, ...)
    message0('droplevel applied.')
  } else{
    r = x
    message0('droplevel did not applied.')
  }
  return(r)
}

# remove zero count categories (filter on sex and genotype)
RemoveZerovarCategories = function(x,
                                   depVar,
                                   minvar = 0,
                                   sex = 'sex',
                                   genotype = 'biological_sample_group',
                                   sep = '_',
                                   drop = TRUE,
								   method = 'MM') {
  note   = NULL
  # do not move me
  if (any(dim(x) == 0) || is.null(complete.cases0(x[, depVar])))
    return(list(x = NULL, note = note))

  #x      = droplevels0(x)
  x      = x[complete.cases0(x[, depVar]), , drop = FALSE]
  newx   = x
  if (is.numeric(x[, depVar]) && (method %in% 'MM')) {
    lvls   = interaction(x[, sex], x[, genotype], sep = sep, drop = drop)
    vars   = tapply(x[, depVar], INDEX = lvls, function(xx) {
      if (is.numeric(xx) && length(na.omit(xx)) > 1) {
        r = var(xx, na.rm = TRUE)
      } else{
        r = NULL
      }
      return(r)
    })

    if (!is.null(vars)) {
      Zvars  = names(vars)[vars <= minvar]
      newx   = x[!lvls %in% Zvars,]
      if (!identical(droplevels0(newx), droplevels0(x))) {
        note$sex_genotype_zero_variation = paste0(
          'Some categories have zero variation and then removed from the raw data. List of zero variation categories: ',
          paste0(Zvars, collapse = ', '),
          ''
        )
      }
    }
  }else{
    message0('Checking variation only applied to the MM method. The input mmethod: ',
             method)
			 }
  if (is.null(x) ||
      is.null(newx) || any(dim(newx) == 0) || any(dim(x) == 0))
    newx = NULL

  return(list(x = newx, note = note))
}




SummaryStatisticsOriginal = function(x,
                                     depVar,
                                     sex = 'sex',
                                     genotype = 'biological_sample_group',
                                     label = 'raw_data_summary_statistics',
                                     lower = FALSE,
                                     drop = TRUE,
                                     sep = '_',
                                     removeSpecialChars = FALSE,
                                     replace = '_') {
  r   = NULL
  # do not move me
  if (any(dim(x) == 0)|| is.null(complete.cases0(x[, depVar])))
    return('empty dataset')

  x = x[complete.cases0(x[, depVar]), , drop = FALSE]
  #x      = droplevels0(x)
  if (is.numeric(x[, depVar])) {
    lvls   = interaction(x[, sex], x[, genotype], sep = sep, drop = drop)
  } else{
    lvls   = interaction(x[, sex], x[, genotype], x[, depVar], sep = sep, drop = drop)
  }

  isNumeric = is.numeric(x[, depVar])
  summaryT   = as.list(tapply(x[, depVar], INDEX = lvls, function(xx) {
    if (isNumeric) {
      c  = ifelse(length(na.omit(xx)) > 0, length(na.omit(xx)), 0)
      m  = ifelse(length(na.omit(xx)) > 0, mean(xx, na.rm = TRUE), NA)
      sd = ifelse(length(na.omit(xx)) > 1, sd(xx, na.rm = TRUE)  , NA)
      r = list(
        count = c                       ,
        mean = m                        ,
        sd = sd                         ,
        normality_test = ifelse(
          length(xx)           > 3    &&
            length(unique(xx)) > 3    &&
            length(xx)         < 5000 &&
            var(xx)            != 0,
          shapiro.test(xx)$p.value,
          'Not possible(Possible causes: <3 or >5000 unique data points'
        )
      )
    } else{
      c = ifelse(length(na.omit(xx)) > 0, length(na.omit(xx)), 0)
      r = list(count = c)
    }
    return(r)
  }, default = -999.991233210123))
  ##
  fTmp = function(isNum) {
    if (isNum) {
      r = list(count = 0,
               mean = NA,
               sd = NA)
    } else{
      r = list(count = 0)
    }
    return(r)
  }
  summaryT[summaryT %in% c(-999.991233210123)] = fTmp(isNum = isNumeric)

  if (lower)
    nnames = tolower(names(summaryT))
  else
    nnames =  names(summaryT)

  if (removeSpecialChars) {
    nnames = RemoveSpecialChars(nnames, replaceBy = replace)
  }

  r = list(lbl = summaryT)
  names(r) = label
  return(r)
}


# Other columns from the data
OtherExtraColumns = function(obj,
                             ColNames,
                             names,
                             lower = TRUE,
                             what = '[^0-9A-Za-z]',
                             replaceBy = '_') {
  p4 = NULL
  ColNameInd = ColNames %in% colnames(obj)
  newNames = ColNames[ColNameInd]
  if (all(dim(obj) > 0) && length(newNames) > 0) {
    if (lower)
      newNames = tolower(newNames)
    p0 = RemoveSpecialChars(names, what = what, replaceBy = replaceBy)
    p4 = obj[, newNames, drop = FALSE]
    p4 = as.list(p4)
    names(p4) = p0[ColNameInd]
  } else{
    p4$Error_in_extra_columns = "Column names do not exist or the dataset is empty"
  }
  return(p4)
}


# min including null/empty
min0 = function(x, ...) {
  if (!is.null(x) && length(na.omit(x)) > 0)
    r = min(x, ...)
  else
    r = 0

  return(r)
}

max0 = function(x, ...) {
  if (!is.null(x) && length(na.omit(x)) > 0)
    r = max(x, ...)
  else
    r = 0
  return(r)
}


# Get possible categories for the categorical variables
GetPossibleCategories = function(procedure = NULL, file = TRUE) {
  requireNamespace('pingr')
  if (file) {
    path = file.path(system.file("extdata", package = "DRrequired"),
                     'AllCts.csv')
    message0('Reading the Category levels from the file : \n\t\t ====> ',
             path)
    CatList = read.csv(file = path)[, -1]
  } else{
    if (pingr::is_online()) {
      url = paste0(
        #'http://ves-ebi-d0:8986/solr/pipeline/select?q=procedure_stable_id:',
        'https://wwwdev.ebi.ac.uk/mi/impc/dev/solr/pipeline/select?q=procedure_stable_id:',
        ifelse(is.null(procedure), '*', paste0('"', procedure, '"')),
        '&fq=categories%3A*&fl=parameter_stable_id%2Ccategories&wt=csv&indent=true&rows=10000000'
      )
      message0('Reading the Category levels from : \n\t\t ====> ',
               url)
      CatList = read.csv(file = url)
      CatList = CatList[!duplicated(CatList), ]
      if (any(dim(CatList) == 0)) {
        CatList[1,] = NA
      } else{
        message0(nrow(CatList),
                 ' parameter(s) found...')
      }
    } else{
      message0('Warning. Please check the internet connection ....')
      CatList = data.frame(parameter_stable_id = NA, categories = NA)
    }
  }
  return(CatList)
}

# Find levels for the parameter
ReadFactorLevelsFromSolr = function(parameter,
                                    CatList,
                                    initReplacePattern = '~~~~',
                                    skip = '\\,',
                                    sep = ',',
                                    fixed = TRUE) {
  note = NULL
  ind = which(CatList$parameter_stable_id == parameter)
  if (is.null(CatList) || any(dim(CatList) < 1) || length(ind) < 1)
    return(list(levels = NULL, note = note))

  if (length(ind) > 1) {
    message0(
      'It seems there are two the same procedure with the same levels. Please check the dataset'
    )
  }

  catCorr = gsub(
    pattern =  skip,
    x = as.character(CatList$categories[ind]),
    replacement = initReplacePattern,
    fixed = fixed
  )

  Levels = unlist(strsplit(x = catCorr, split = sep, fixed = fixed))
  fLevels = gsub(
    pattern = initReplacePattern,
    replacement = skip,
    x = Levels,
    fixed = fixed
  )
  return(list(levels = fLevels, note = note))
}

# Fast replacement of nulls
replaceNull = function(x, replace = '-') {
  if (length(x) > 1) {
    x[is.null(x)] == replace
  } else{
    x = replace
  }
  return(x)
}

mapLevelsToFactor = function(levels, newlevels, name = 'response') {
  #############
  res  = NULL
  n1   = paste0('all_levels_in_',
                name,
                collapse = '')
  res$n1 = newlevels
  ind    = which(newlevels %in% levels)
  if (length(ind) < 1 ||
      length(newlevels) < 0 ||
      length(levels) < 0 ||
      is.null(newlevels) ||
      is.null(levels) ||
      length(ind) == length(levels)) {
    if (!is.null(res$n1))
      names(res) = n1
    return(list(levels = levels, note = res))
  }
  ############
  n0 = paste0('missing_levels_in_',
              name,
              collapse = '')
  res$n0 = newlevels[-ind]
  names(res) = c(n1, n0)

  flevels = c(levels, newlevels[-ind])
  return(list(levels = flevels, note = res))
}



RemoveNonExsistingFilesFromList = function(fileList) {
  logical = file.exists(fileList)
  if (sum(!logical) > 0) {
    message0('Non-exsisting file(s): \n\t\t ~> ',
             paste(fileList[!logical], collapse  = '\n\t\t ~> '))
  }
  return(fileList[logical])
}
# Compress a file
compressFiles = function(fileList,
                         dir,
                         filename,
                         overwrite = TRUE,
                         flags = '-D -j',
                         rmSource = TRUE,
                         quiet = TRUE) {
  if (length(fileList) < 1 ||
      length(RemoveNonExsistingFilesFromList(fileList)) < 1) {
    message0('No file to compress or files do not exists')
    return(NULL)
  }
  if (quiet)
    flags = paste0(flags, ' -q')
  fileList      = RemoveNonExsistingFilesFromList(fileList)
  fileListPaste = paste('\n\t\t ~~> ', fileList, sep = '')
  zipFile =  file.exists0(file.path(dir, paste0(filename,  '.zip')), overwrite = overwrite)
  MB      = 1024
  message0('Compressing file ...')
  result = zip(zipfile = zipFile,
               files = fileList,
               flags = flags)
  if (file.exists(zipFile) &&
      (file.size(zipFile) / MB) > 0 &&
      result == 0) {
    message0('Successfully compressed. file(s): ', fileListPaste)
    if (rmSource)
      unlink(fileList)
  } else{
    message0('Failed in compressing file. file(s): ', fileListPaste)
  }
  return(list(status = result, file = zipFile))
}

# Base64 encoder
base64 = function(x, active) {
  if (active) {
    message0('Based64 encryption as well as (gzip) compression in progress ...')
    r = base64enc::base64encode(memCompress(x, "g"))
  } else{
    r = x
  }
  return(r)
}



# rename if file exists (Danger! if overwrite==TRUE then it removes the files!)
file.exists0 = function(file, overwrite = FALSE, ...) {
  counter = 1
  ###
  fname = basename(file)
  dotLocation = which(strsplit(fname, "")[[1]] == ".")
  if (length(dotLocation) > 0) {
    name = substr(x = fname,
                  start = 0,
                  stop = dotLocation - 1)

    ext = substr(x = fname,
                 start =  dotLocation,
                 stop = nchar(fname))

  } else{
    name = fname
  }
  ###
  if (overwrite && file.exists(file, ...)) {
    unlink(file, recursive = FALSE, force = FALSE)
  }
  while (!overwrite && file.exists(file, ...)) {
    file = file.path(dirname(path = file),
                     paste0(name, '_', counter, ext))
    counter = counter + 1
  }
  if (counter > 1) {
    write(
      counter - 1,
      file = file.path0(
        dirname(file),
        'duplicates_detected.txt',
        IncludedFileName = TRUE,
        create = TRUE,
        check = FALSE
      )
    )
  }
  return(file)
}

# add bits to file name
addBitToFileName = function(file,
                            bit = '',
                            sep = '_',
                            fullpath = TRUE) {
  if (fullpath)
    file = paste(dirname(file), bit, basename(file), sep = sep)
  else
    file = paste(bit, basename(file), sep = sep)
  return(file)
}


# plot analised data with window
plot_win = function(phenlistObject, r, depVariable, check, ...) {
  col = pch = as.factor(interaction(
    phenlistObject@datasetPL$Genotype,
    phenlistObject@datasetPL$Sex
  ))
  par(mar = c(5.1, 4.1, 4.1, 8.1))
  requireNamespace('SmoothWin')
  plot(r,
       col = as.integer(col),
       pch = as.integer(pch),
       cex.main = .7,
       ...)
  legend(
    'topright',
    legend = unique(col),
    col = as.integer(unique(col)),
    pch = as.integer(unique(pch)),
    inset = c(-.24, 0),
    xpd = TRUE,
    cex = .7
  )

  if (check == 1) {
    rt = table(phenlistObject@datasetPL$Batch)
    rt2 = rt[which(rt < 2)]
    rdata2 = phenlistObject@datasetPL[phenlistObject@datasetPL$Batch %in% names(rt2), ]
    points(as.numeric(as.Date(rdata2$Batch)),
           rdata2[, depVariable],
           pch = 8,
           col = 'gray')
  }
}

# boxplot for window
boxplot_win = function(phenlistObject, we, threshold, ...) {
  boxplot(
    data_point ~ Genotype + Sex,
    data = phenlistObject@datasetPL[we > threshold, ],
    col = rep(2:3, 2),
    cex.axis = .7,
    cex.main = .7,
    ...
  )
}

# make JSON vector (with [ or not)
SingleOrMore = function(name, value, comma = FALSE) {
  NewVal = ifelse(
    is.numeric(value),
    paste0(value, collapse = ', '),
    paste0('"', value, '"', collapse = ', ')
  )
  if (length(value) > 1) {
    outp = paste0('"', name, '":[', NewVal, ']', collapse = '')
  } else{
    outp = paste0('"', name, '":', NewVal, '', collapse = '')
  }
  if (comma)
    outp = paste(outp, ', ', sep = '', collapse = '')
  return(outp)
}


# cut the filename is it is too long
CutFileNameIfTooLong = function(fullPath,
                                max = 255,
                                checkRoot = TRUE) {
  file = basename(fullPath)
  dir  = dirname(fullPath)
  nch  = nchar(file)
  message0('The number of characters in file name: ', nch)
  k = 7
  if (nch >= max) {
    message0('File name too long. cut 255 characters from the RIGHT')
    file = paste0(
      RandomRegardSeed(
        n = 1,
        decimal = TRUE,
        round = 6
      ) * 10 ^ 6,
      '_',
      substr(x = file, start =  k + 2, nch),
      collapse = ''
    )
  }
  return(ifelse(checkRoot, file.path(dir, file), file.path(file)))
}




# Model and test for errors/messages
ModeWithErrorsAndMessages = function(m2ethod,
                                     key = 'final_method',
                                     phenList,
                                     depVariable,
                                     equation,
                                     method,
                                     modelWeight = NULL,
                                     keepList = NULL,
                                     threshold,
                                     check,
                                     sep = '_',
                                     col = '_',
                                     name = NULL,
                                     ...) {
  note = NULL
  object0 = tryCatch.W.E(
    expr =
      print(
        testDataset(
          phenList = phenList,
          depVariable = depVariable,
          equation = equation,
          method = method,
          modelWeight = modelWeight,
          keepList = keepList,
          threshold = threshold,
          check = check
        )
      ),
    name = name,
    sep = sep,
    col = col
  )
  mess = capture.stderr(try(expr = testDataset(
    phenList = phenList,
    depVariable = depVariable,
    equation = equation,
    method = method,
    modelWeight = modelWeight,
    keepList = keepList,
    threshold = threshold,
    check = check
  ),
  silent = TRUE))
  if (!NullOrError(object0$value))
  {
    note$k1 = paste0(m2ethod, ' is applied')
    names(note)[1] = key
  } else{
    n1 = 'internal_errors'
    n2 = paste(name,
               'phenstat_errors',
               sep = sep,
               collapse = col)

    note$internal_errors = if (length(object0$both) < 1) {
      'No internal error found'
    } else{
      object0$both
    }
    note$phenstat_errors = if (length(mess$output) < 1 ||
                               mess$output == '') {
      'No PhenStat error found'
    } else{
      mess$output
    }
    names(note) = c(n1, n2)
  }
  return(list(obj  =  object0, mess  =  mess, note = note))
}



# plot windowing
PlotWindowingResult = function(args, overwrite = FALSE, ...) {
  if (args$plot) {
    message0('ploting the results ...')
    file = ifelse(
      is.null(args$PicDir),
      file.path(getwd(), args$filename),
      file.path(dir.create0(args$PicDir, recursive = TRUE), args$filename)
    )
    file = file.exists0(paste(file, '.png', sep = ''), overwrite = overwrite)
    file = CutFileNameIfTooLong(file)
    ##########
    if (args$storeplot)
      png(
        filename = file,
        width = 1000,
        height = 1000,
        res = 120
      )
    plot_win(
      phenlistObject = args$phenlistObject,
      main = args$main,
      r = args$r,
      check = args$check,
      depVariable = args$depVariable,
      xaxt = 'n'
    )
    tick =  round(seq(min(args$tt), max(args$tt), length.out = 5))
    axis(1, at = tick, labels = as.Date(tick, '1970-01-01'))
    if (args$storeplot)
      graphics.off()
    # if (args$storeplot && 1 == 0) {
    #   png(filename  = addBitToFileName(file = file, bit = 'boxplot'))
    #   boxplot_win(
    #     phenlistObject = args$phenlistObject,
    #     we = args$we,
    #     threshold = args$threshold
    #   )
    #   if (args$storeplot)
    #     graphics.off()
    # }
    return(file)
  } else{
    return(NULL)
  }
}


## date to integer
Date2Integer = function(x) {
  dates = as.Date(x)
  dates = as.numeric(dates)
  return(dates)
}


# check if variable name exists in the data frame
CheckIfNameExistInDataFrame = function(obj, name, checkLevels = TRUE) {
  if (!NullOrError(obj) && !is.null(name) && length(name) > 0) {
    r = name %in% names(obj)
    if (r &&
        checkLevels &&
        !is.numeric(obj[, name]) &&
        length(unique(obj[, name])) == 1) {
      r = FALSE
    }
  } else{
    r = FALSE
  }
  return(r)
}

# System cores set to zero for 1 core
cores0 = function(coreRatio = .7,
                  activate = TRUE) {
  requireNamespace('parallel')
  if (activate) {
    if (coreRatio < 1) {
      crs = max(ceiling(parallel::detectCores() * min(coreRatio, 1)), 1)
    } else{
      crs = max(min(parallel::detectCores() , coreRatio)  , 1)
    }
  } else{
    crs = 1
  }
  message0(
    'The total number of cores on this machine: ',
    detectCores(),
    '; and ',
    crs,
    ' would be used.'
  )
  return(crs)
}


# file.path to check for space and special characters in the path
file.path0 = function(...,
                      check = TRUE,
                      what = '[^0-9A-Za-z/:]',
                      replace = '_',
                      create = TRUE,
                      IncludedFileName = FALSE) {
  path = file.path(...)
  if (check) {
    path = RemoveSpecialChars (path, what = what, replaceBy = replace)
  }
  if (create) {
    if (!IncludedFileName) {
      dir.create0(path, recursive = TRUE)
    } else{
      dir.create0(dirname(path), recursive = TRUE)
    }
  }
  return(path)
}




# exist in list
IsInList = function(item = NULL,
                    list = NULL ,
                    message = '') {
  r = FALSE
  if (length(list) < 1 || length(item) < 1)
    return(r)
  if (any(list %in% item)) {
    r = TRUE
    message0(message)
  }
  return(r)
}


### Multicore log
outMCoreLog = function(wd, dir = 'Multicore_logs', fname = '_MulCoreLog.txt') {
  path = file.path0(
    wd,
    dir,
    paste0(
      Sys.Date(),
      '_',
      RandomRegardSeed(1),
      '_',
      RemoveSpecialChars(paste(head(Sys.info(

      ), 3), collapse = '-')),
      fname
    ),
    check = FALSE,
    create = TRUE,
    IncludedFileName = TRUE
  )
  return(path)
}

sortList = function(x,...){
  x[order(names(x),...)]
}

# A new vector output for this package only
VectorOutput0 = function(c.ww0,
                         ExtraCols,
                         activeWindowing,
                         na = 'string',
                         null = 'null') {
  ####
  # The first piece (Normal results)
  funAux1 = function(x) {
    r =  paste0(paste(
      '"',
      VectorOutPutNames(),
      '"',
      paste(':',
            RemoveTwoSpecialFromOutput(x),
            '',
            sep = ''),
      sep = ''
    ),
    collapse = ',')
    return(r)
  }
  p1 = ifelse (!NullOrError(c.ww0$NormalObj$value),
               funAux1(
                 PhenStat::vectorOutput(c.ww0$NormalObj$value,
                                        othercolumns = ExtraCols)
               ),
               '')
  # The second piece (Windowing results)
  if (activeWindowing && !NullOrError(c.ww0$WindowedObj$value)) {
    p2 =  ifelse(
      !is(c.ww0$WindowedObj$value, 'simpleWarning'),
      funAux1 (
        PhenStat::vectorOutput(c.ww0$WindowedObj$value  ,
                               othercolumns = ExtraCols)
      ),
      ''
    )
    p3 =  ifelse(!is(c.ww0$FullObj$value, 'simpleWarning'),
                 funAux1 (
                   PhenStat::vectorOutput(c.ww0$FullObj$value  ,
                                          othercolumns = ExtraCols)
                 ),
                 '')
    p4 =  ifelse(
      !is(c.ww0$FullWindowedObj$value, 'simpleWarning'),
      funAux1 (
        PhenStat::vectorOutput(c.ww0$FullWindowedObj$value  ,
                               othercolumns = ExtraCols)
      ),
      ''
    )
  } else{
    p2 = p3 = p4 = ''
  }
  ###
  output = paste(
    '"normal_result":{',
    p1
    ,
    '}, "windowed_result":{',
    p2,
    '}, "full_model_result":{',
    p3,
    '}, "full_model_windowed":{',
    p4,
    '}',
    collapse = '',
    sep      = ''
  )

  list = RJSONIO::fromJSON(
    paste0('{', output, '}'),
    stringFun = function(x) {
      if (x %in% c('NA', 'TRUE', 'FALSE')) {
        checkQouteNAandNaN(
          pattern     = c('NA', 'TRUE', 'FALSE', '"TRUE"', '"FALSE"'),
          replacement = c(NA  ,  TRUE ,  FALSE,    TRUE,     FALSE),
          x = x
        )
        if (x %in% c('TRUE', 'FALSE'))
          x = as.logical(x)
      } else{
        x
      }
    },
    simplify = TRUE,
    simplifyWithNames = FALSE
  )
  jsonF = RJSONIO::toJSON(x = list,
                          #pretty = TRUE,
                          na = na,
                          null = null)
  return(list(
    json = jsonF,
    Initialjson = output,
    list = list
  ))
}


# Area under the curve
AUC = function(x, y) {
  lx = length(x)
  s = 0
  if (lx > 1) {
    for (i in 1:(lx - 1)) {
      s = s + abs(y[i] + y[i + 1]) * abs(x[i + 1] - x[i]) / 2
    }
  }
  return(s)
}

# Chunk a vector of in particular colonies into chunks (used mainly in Batch production)
chunkVector =  function(x,
                        n = 10,
                        min = 25,
                        activate = TRUE) {
  r = list()
  lx = length(x)
  if (lx > min && n > 1 && activate) {
    r = split(x, cut(seq_along(x), n, labels = FALSE))
  } else{
    r$'1' = x
  }
  return(r)
}

# Plan B for JSON not quouted keys!
checkQouteNAandNaN = function(pattern, replacement, x, ignoreCase = FALSE) {
  for (i in 1:length(pattern)) {
    pt  = pattern[i]
    rep = replacement[i]
    x = gsub(
      pattern = pt,
      replacement = rep,
      x = x,
      fixed = TRUE,
      ignore.case = ignoreCase
    )
  }
  return(x)
}







# Add windowing weights and store the data again
StoreRawDataAndWindowingWeights = function(storeRawData,
                                           activeWindowing        ,
                                           c.ww0                  ,
                                           ## 2
                                           orgData                ,
                                           RawoutputFile          ,
                                           files                  ,
                                           dir                    ,
                                           filename               ,
                                           ### 3
                                           ReadMeTxt              ,
                                           ReadMeFile             ,
                                           ### 4
                                           colnames = c('external_sample_id', 'AllModelWeights'),
                                           byid     = 'external_sample_id'                      ,
                                           compressRawData = TRUE) {
  if (storeRawData &&
      activeWindowing && !NullOrError(c.ww0$WindowedObj$value)
      && (c.ww0$method %in% 'MM')) {
    n3.5.w = merge(
      x = orgData,
      y = c.ww0$InputObject@datasetPL[, colnames],
      by = byid,
      all = TRUE
    )
    message0('writting the raw data + windowing weights to disk ... \n \t\t===> ',
             RawoutputFile)
    write.csv(x = n3.5.w,
              file = RawoutputFile,
              row.names = FALSE)
    write(x        = ReadMeTxt ,
          file     = ReadMeFile,
          ncolumns = 1)

    if (compressRawData) {
      comRes = compressFiles(
        fileList = files    ,
        dir = dir           ,
        filename = filename ,
        overwrite = TRUE
      )
      if (comRes$status == 0)
        message0('[data after adding the windowing weight column] Compression successful')
      else
        message0('[data after adding the windowing weight column] Compression failed')
    }
  }
}

# new message with timestamp in the begining
message0 = function(...,
                    ShowTimeTag = TRUE,
                    active = TRUE) {
  if (active) {
    if (ShowTimeTag)
      message(Sys.time(), '. ', ...)
    else
      message(...)
  }
}

## ReadMe
ReadMe = function(obj, URL = NULL, skip = NULL) {
  if (any(dim(obj)) > 0 && !is.null(obj)) {
    ReadMeTxt = paste(
      c(
        'Gene_symbol',
        'Procedure',
        'Parameter',
        'Centre',
        'Strain',
        #'Metadata',
        'Zygosity',
        'Colony_id',
        'Metadata_hash',
        'URL'
      ),
      c(
        UniqueAndNNull(obj$gene_symbol,removeSpecials = FALSE)         ,
        UniqueAndNNull(obj$procedure_group,removeSpecials = FALSE)     ,
        UniqueAndNNull(obj$parameter_stable_id,removeSpecials = FALSE) ,
        UniqueAndNNull(obj$phenotyping_center,removeSpecials = FALSE)  ,
        UniqueAndNNull(obj$strain_accession_id,removeSpecials = FALSE) ,
        #UniqueAndNNull(obj$metadata,removeSpecials = FALSE)           ,
        UniqueAndNNull(obj$zygosity,removeSpecials = FALSE)            ,
        UniqueAndNNull(obj$colony_id,removeSpecials = FALSE)           ,
        UniqueAndNNull(obj$metadata_group,removeSpecials = FALSE)      ,
        UniqueAndNNull(URL,removeSpecials = FALSE)
      ),
      sep = ' = '
    )
    message0('[',
             paste(if (is.null(skip)) {
               ReadMeTxt
             } else{
               ReadMeTxt[-skip]
             },
             collapse =  ']~>['),
             ']\n')
  } else{
    ReadMeTxt = 'Not found!'
  }

  return(ReadMeTxt)
}




# Replace in list
replaceInList <- function (x, FUN, ...)
{
  if (is.list(x)) {
    for (i in seq_along(x)) {
      x[i] <- list(replaceInList(x[[i]], FUN, ...))
    }
    x
  }
  else
    FUN(x, ...)
}



FinalJsonBobectCreator = function(FinalList,
                                  null = 'null',
                                  na = 'null' ,
                                  auto_unbox = TRUE,
                                  SpecialString = '==!!(:HAMED:)!!==',
                                  rep = 3) {
  message0('Forming the JSON object ...')
  FinalList = replaceInList(
    FinalList,
    FUN  = function(x) {
      if (is.null(x)) {
        x = SpecialString
      } else{
        x
      }
    }
  )
  for (i in 1:rep) {
    JsonObj  = jsonlite::toJSON(
      x = FinalList,
      null = null,
      na   = na,
      digits = NA,
      auto_unbox = ifelse(i == rep, auto_unbox, FALSE)
    )
    FinalList = jsonlite::fromJSON(txt =   JsonObj)
  }
  JsonObj = gsub(
    pattern = paste0('"', SpecialString, '"'),
    replacement = '{}',
    x = JsonObj,
    ignore.case = FALSE,
    fixed = TRUE
  )
  # Remove new lines! e.g. in Corneal Reflex
  JsonObj = gsub(pattern = '\n',
                 replacement = ' ',
                 x = JsonObj)
  return(JsonObj)
}



UnzipAndfilePath = function(file, quiet = TRUE, order = TRUE) {
  if (!grepl(x = basename(file),
             pattern = '.zip',
             fixed = TRUE)) {
    message0('It is not a zip file:')
    message0('\t', file)
    return(file)
  }
  message0('Unziping file ...')
  # get the file url
  forecasturl = file
  # create a temporary directory
  td = file.path0(
    tempdir(),
    paste0(sample(LETTERS, 4, replace = TRUE), collapse = ''),
    check = FALSE,
    IncludedFileName = FALSE,
    create = TRUE
  )
  # create the placeholder file
  tf =  tempfile(
    pattern = paste0(sample(LETTERS, 4, replace = TRUE), collapse = ''),
    tmpdir = td,
    fileext = ".zip"
  )
  if (file.exists(forecasturl)) {
    file.copy(from = forecasturl,
              to =  tf,
              overwrite = TRUE)
  } else{
    # download into the placeholder file
    download.file(forecasturl, tf, quiet = quiet)
  }
  # get the name of the first file in the zip archive
  fname = unzip(tf, list = TRUE)$Name
  # unzip the file to the temporary directory
  unzip(tf,
        files = fname,
        exdir = td,
        overwrite = TRUE)
  # fpath is the full path to the extracted file
  fpath = file.path(td, fname)
  if (order) {
    fpath = fpath[order(nchar(fpath), fpath)]
  }
  message0(
    'Successfully unziped. List of files:\n',
    paste0('\t\t', 1:length(fpath), ' ==> ', fpath, collapse = '\n')
  )
  return(fpath)
}

# Create the relative path from the full path
relativePath = function(path, reference) {
  if (is.null(path))
    return(NULL)
  ######
  p = gsub(
    pattern = reference,
    replacement = '/',
    x = path,
    fixed = TRUE
  )
  p = gsub(
    pattern = '//',
    replacement = '/',
    x = p,
    fixed = TRUE
  )
  return(p)
}




# Closest points for the mutants
closest.time.value  = function(x,
                               time1,
                               y,
                               time2,
                               yind,
                               TimeOnly = TRUE,
                               maxIter  = 200) {
  requireNamespace('abind')
  message0('Resampling. Creating the index ...')
  x.df = data.frame(x.val = x, x.time = time1)
  y.df = data.frame(y.val = y,
                    y.time  = time2,
                    y.index = yind)
  message0('Resampling. Finding the distance of values from the mutants ...')
  output2  = ol = lapply(1:nrow(x.df), function(i) {
    tt <-
      cbind(x.df[i,],
            lapply(x.df[i,]$x.val, function(v) {
              diff <- abs(y.df$y.val - v)
              y.df$dist.V = diff
              out <- y.df
            }),
            ind = i,
            row.names = NULL)
    tt$dist.T <- abs(tt$x.time - tt$y.time)
    tt$totalD  = tt$dist.V + tt$dist.T
    if (TimeOnly) {
      tt = tt[sample(1:nrow(tt)),]
    } else{
      tt = tt[order(tt$totalD)  ,]
      tt = tt[order(tt$dist.V)  ,]
    }
    tt   = tt[order(tt$dist.T)  ,]
  })
  dol = 1
  message0('Resampling. Refining ...')
  counter = 1
  while (sum(dol) > 0 && counter < maxIter) {
    ol  = lapply(
      X = output2,
      FUN = function(x) {
        if (!is.null(x)  && nrow(x) > 0) {
          x[1,]
        } else{
          NULL
        }
      }
    )
    ol2  = abind(ol, along = 1)
    dol  = duplicated(ol2[, 'y.index'])
    if (sum(dol)) {
      output2[dol] = lapply(
        output2[dol],
        FUN = function(x) {
          x[-1, , drop = FALSE]
        }
      )
    }
    counter = counter + 1
  }
  ####################
  message0('Resampling. Creating output ...')
  return(as.data.frame(abind(ol[!unlist(lapply(
    ol,
    FUN = function(x) {
      is.null(x) || length(x) < 1
    }
  ))], along = 1)))
}

Factor2CharAndSubstitution = function(df,
                                      Ind,
                                      column,
                                      ArtifLabel,
                                      neutralise,
                                      mutInd,
                                      baselines) {
  df[, column]         = as.character(df[, column])
  df[Ind, column]      = ArtifLabel
  if (neutralise)      {
    df[mutInd, column] = baselines
  }
  df[, column]         = as.factor(df[, column])
  return(df)
}

mimicControls = function(df                             ,
                         depVariable = 'data_point'     ,
                         ArtifLabel = 'ArtificialMutant',
                         mutLabel = 'experimental'      ,
                         baselines = 'control'          ,
                         sex       = 'sex'              ,
                         neutralise = TRUE              ,
                         removeMutants = FALSE          ,
                         resample  = FALSE              ,
                         minSampRequired                ,
                         SexGenResLevels                ,
                         indicator       = NULL         ,
                         plot            = FALSE        ,
                         CutIfLongerThanMutants = TRUE) {
  if (!is.null(indicator))
    set.seed(indicator)
  requireNamespace('stringi')
  df.bckOrg    = df
  note         = list(shift = NULL, new_mutant_indices = NULL, original_mutant_indices = NULL)
  shift        = Indices = NULL
  #####
  # remove zero frequency categories
  df_rzeros = RemoveZeroFrequencyCategories(
    x = df,
    minSampRequired = minSampRequired,
    depVar = depVariable,
    totalLevels = SexGenResLevels
  )
  df                             = df_rzeros$x
  note$removed_categories_detail = df_rzeros$note
  ###########
  if (is.null(df)                ||
      length(df) < 1             ||
      nrow  (df) < 1             ||
      nlevels(df[, sex]) < 1     ||
      !is.numeric(df[, depVariable])) {
    note$dataset_status = 'Empty dataset'
    return(list(
      df = df          ,
      ctv = NULL       ,
      input = df.bckOrg,
      note  = note
    ))
  }
  # Definind id is crutial otherwise it may end up final duplicated points
  # id must be numeric otherwise will change the entire output to factor! Stu.. R
  # Do not include 0 to have the same length always 021=21!
  df$id_d   = as.numeric(stri_rand_strings(nrow(df), length = 10, pattern = '[1-9]'))
  df$cTime  = Date2Integer(df$date_of_experiment)
  mutInd    = (df$biological_sample_group == mutLabel)
  conInd    = (df$biological_sample_group == baselines)
  if (!sum(mutInd)  ||
      !sum(conInd)) {
    note$dataset_status          = 'Missing controls or mutants'
    note$original_mutant_indices =  mutInd
    return(list(
      df    = df.bckOrg,
      ctv   = NULL     ,
      input = df.bckOrg,
      note  = note
    ))
  }
  ##################
  for (l in levels(df[, sex])) {
    mutIndSex    = (df$biological_sample_group == mutLabel  &
                      df[, sex] == l)
    conIndSex    = (df$biological_sample_group == baselines &
                      df[, sex] == l)
    ###
    y     = df[conIndSex, depVariable]
    time2 = df[conIndSex, 'cTime']
    yid   = df[conIndSex, 'id_d']
    x     = df[mutIndSex, depVariable]
    lx    = length(x)
    if (is.null(shift)) {
      message0(l,', Generating new SHIFT in progress ...')
      shift = ifelse(resample, round(runif(
        1                                             ,
        min(df$cTime)[1] - min(df$cTime[mutIndSex])[1],
        max(df$cTime)[1] - max(df$cTime[mutIndSex])[1]
      )), 0)
    }
    message0('\t ~> ', l, '\t ~~> Shift = ', round(shift))
    time1 = df$cTime[mutIndSex] + shift
    ####
    ctv = closest.time.value(
      x     = x      ,
      time1 = time1  ,
      y     = y      ,
      time2 = time2  ,
      yind  = yid
    )
    if (!is.null(ctv) && CutIfLongerThanMutants && nrow(ctv) > lx) {
      message0 ('The length of the simulated data is longer than the input. Cutting in progress ....')
      ctv = ctv[1:lx, ]
    }
    Indices = c(Indices, ctv$y.index)
  }
  Ind                            =  df$id_d %in% Indices
  note$new_mutant_indices        =  which(Ind)
  note$shift                     = shift
  note$original_mutant_indices   =  mutInd
  note$dataset_status            = 'No problem detected'
  # Replace the biological_sample_group and colony_id
  df$colony_id[Ind] = paste0(unique(na.omit(df$colony_id)), collapse = '')
  df = Factor2CharAndSubstitution(
    df = df,
    Ind = Ind,
    column = 'biological_sample_group',
    ArtifLabel = ArtifLabel,
    neutralise = neutralise,
    mutInd = mutInd,
    baselines = baselines
  )
  ####
  if (plot) {
    plot(
      time2,
      y,
      main = paste0('Shift = ', shift)     ,
      xlim = c(min(ctv$y.time, ctv$x.time) ,
               max(ctv$y.time, ctv$x.time)),
      type = 'p'
    )
    points(df$cTime     ,
           df[, depVariable],
           #col = as.integer(df$biological_sample_group) - 1,
           col = as.integer(as.factor(interaction(df$biological_sample_group,df$sex)))
           )
    points(
      df$cTime[mutInd],
      df[mutInd, depVariable],
      pch = paste0(1:length(mutInd)),
      cex = 1,
      lwd = 1,
      col = 2
    )
    points(
      time1,
      x,
      pch = paste0(1:length(x)),
      cex = 1.5,
      lwd = 2,
      col = 3
    )
    points(
      ctv$y.time,
      ctv$y.val,
      pch = paste0(1:length(ctv$y.val)),
      cex = 2,
      lwd = 3,
      col = 4
    )
    legend(
      'top',
      c('TRUE', 'Shifted Mut', 'Closest'),
      fill = c(2, 3, 4),
      horiz = TRUE
    )
  }
  # Must be in the end
  if (removeMutants) {
    df = df[!mutInd,]
  }

  df$id_d = df$cTime = NULL
  return(list(
    df = droplevels(df),
    ctv = ctv,
    input = df.bckOrg,
    note = note
  ))
}
