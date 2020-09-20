file = commandArgs(trailingOnly = TRUE)
library('RPostgreSQL')
library('data.table')
########################### Annotation pipeline #################################
##############################
library(data.table)
library(jsonlite)
library(rlist)
library(Tmisc)
##############################
GetMethodStPa <- function(x) {
  if (is.null(x)) {
    return("UNK")
  }
  if (grepl(pattern = "Fisher", x)) {
    r <- "FE"
  } else if (grepl(pattern = "(Mixed Model)|(Linear Model)", x)) {
    r <- "MM"
  } else if (grepl(pattern = "Reference Range", x)) {
    r <- "RR"
  } else {
    r <- "UNK"
  }
  return(r)
}

DirectionTagMM = function(x,
                          pvalue = NULL,
                          threshold = .0001) {
  tag = if (is.null(pvalue)) {
    'NoPValueAvailable'
  } else if (pvalue > threshold) {
    paste0('Test not significant at the level of ',
           threshold,
           ' (pvalue = ',
           pvalue,
           ')')
  }  else if (is.null(x) || length(x) < 1) {
    'NoEffectCalculated'
  } else if (x > 0 && pvalue < threshold) {
    'INCREASED'
  } else if (x < 0 && pvalue < threshold) {
    'DECREASED'
  } else if (x == 0 && pvalue < threshold) {
    'STEADY'
  } else{
    'Something is gone wrong. Please report this case to the IMPC team.'
  }
  return(tag)
}

DirectionTagFE = function(x,
                          threshold = .0001,
                          group = c('ABNORMAL', 'INCREASED', 'DECREASED')) {
  default = c('ABNORMAL', 'INCREASED', 'DECREASED')
  if (length(group) < 1) {
    message(
      'No group for the so called ABNORMAL category provided, swith to the backup plan:\n\t ',
      default
    )
    group = default
  }
  tag = if (is.null(x) || length(x) < 1) {
    'NoEffectCalculated'
  } else if (x < threshold) {
    group
  } else {
    paste0('NORMAL - Test not significant at the level of ',
           threshold,
           ' (pvalue = ',
           x,
           ')')
  }
  return(tag)
}

listM = function(x, What2Attach = list('OVERALL' = list('MPTERM' = NA))) {
  if (is.null(x))
    return(NULL)
  al = as.list(x)
  names(al) = al
  al[[names(al)]] = What2Attach
  return(al)
}


numbers_only <- function(x) {
  r = !grepl("\\D", x)
  return(r)
}

multiGrepl = function (x = NULL, pattern = NULL, ...)
{
  if (is.null(x)) {
    return(NULL)
  }
  if (is.null(pattern)) {
    return(TRUE)
  }
  r <- rep(TRUE, length(x))
  for (p in pattern) {
    r <- r & grepl(pattern = p, x = x, ...)
  }
  return(r)
}

returnWhatBasedOnThreshold = function(x = NULL,
                                      threshold = .0001,
                                      Return = 'ABNORMAL',
                                      ReturnNull = 'IgnoreThisCaseAtALL') {
  if (is.null(x)      ||
      is.null(Return) ||
      is.null(threshold) ||
      length(x) < 1      ||
      length(Return) < 1 ||
      length(threshold) < 1) {
    return(ReturnNull)
  } else if (is.numeric(x) && x < threshold) {
    return(Return)
  } else{
    return(ReturnNull)
  }
}

GenotypeTag = function(obj,
                       parameter_stable_id = NULL,
                       threshold = 10 ^ -4,
                       expDetailsForErrorOnly = NULL) {
  if (is.null(obj))
    return('NotAnalysed')
  method = GetMethodStPa(x = obj$`Applied method`)
  message('\t The analysis method = ', method)
  message('\t The decision threshold = ', threshold)
  if (method %in% 'MM') {
    tag = list(
      UNSPECIFIED = listM(
        DirectionTagMM(
          x      = obj$`Genotype estimate`$Value,
          pvalue = obj$`Genotype p-value`,
          threshold = threshold
        )
      ),
      UNSPECIFIED = listM(
        returnWhatBasedOnThreshold(
          x = obj$`Genotype p-value`,
          threshold = threshold,
          Return = 'ABNORMAL'
        )
      ),
      UNSPECIFIED = listM(
        returnWhatBasedOnThreshold(
          x = obj$`Genotype p-value`,
          threshold = threshold,
          Return = 'INFERRED'
        )
      ),
      FEMALE        = listM(
        DirectionTagMM(
          x      = obj$`Sex FvKO estimate`$Value,
          pvalue = obj$`Sex FvKO p-value`,
          threshold = threshold
        )
      ),
      FEMALE = listM(
        returnWhatBasedOnThreshold(
          x = obj$`Sex FvKO p-value`,
          threshold = threshold,
          Return = 'ABNORMAL'
        )
      ),
      FEMALE = listM(
        returnWhatBasedOnThreshold(
          x = obj$`Sex FvKO p-value`,
          threshold = threshold,
          Return = 'INFERRED'
        )
      ),
      MALE      = listM(
        DirectionTagMM(
          x      = obj$`Sex MvKO estimate`$Value,
          pvalue = obj$`Sex MvKO p-value`,
          threshold = threshold
        )
      ),
      MALE = listM(
        returnWhatBasedOnThreshold(
          x = obj$`Sex MvKO p-value`,
          threshold = threshold,
          Return = 'INFERRED'
        )
      ),
      MALE = listM(
        returnWhatBasedOnThreshold(
          x = obj$`Sex MvKO p-value`,
          threshold = threshold,
          Return = 'INFERRED'
        )
      )
    )
  } else if (method %in% 'FE') {
    ###########################################################################
    fmodels = obj$`Additional information`$Analysis$`Further models`$category
    if (is.null(fmodels))
      return(NULL)

    if (length(names(fmodels)) == 1 &&
        names(fmodels) == 'Complete table') {
      fmodels$Genotype$`Complete table` = fmodels$`Complete table`
      fmodels$`Complete table` = NULL
    }
    #fmodels$Genotype$`Complete table`
    AllCombinations = lapply(fmodels, function(x) {
      if (length(fmodels) > 0) {
        lapply(x, function(y) {
          DirectionTagFE(x = y$p.value, threshold = threshold)
        })
      } else{
        DirectionTagFE(x = x$p.value, threshold = threshold)
      }
    })
    if (is.null(AllCombinations))
      return(NULL)
    #### Make the list as sequence of names attached with dot (.)
    AllCombinations1 = unlist(AllCombinations)
    if (!is.null(AllCombinations1)){
      AllCombinations1 = AllCombinations1[grepl('Complete table',names(AllCombinations1))]
    }
    if (is.null(AllCombinations1))
      return(NULL)
    #### Keep only genotype analysis
    AllCombinations2 = AllCombinations1[grepl(pattern = 'Genotype', x =
                                                names(AllCombinations1))]
    #### The abnormal case is made from the data
    for (i in seq_along(AllCombinations2)) {
      ############# step 1
      names(AllCombinations2)[i] = gsub(
        pattern = 'Complete table',
        replacement = paste0(
          AllCombinations2[i],
          ifelse(
            AllCombinations2[i] == 'ABNORMAL',
            '.OVERALL.MPTERM',
            '.MPTERM'
          )
        ),
        x = names(AllCombinations2)[i]
      )
      ############# step 2
      gind = !multiGrepl(pattern = c('Genotype', 'OVERALL'),
                         x = names(AllCombinations2)[i])
      if (gind) {
        strs = unlist(strsplit(
          names(AllCombinations2)[i],
          split = '.',
          fixed = TRUE
        ))
        u = unique(c(strs[1],
                     AllCombinations2[i],
                     strs[-1]))
        names(AllCombinations2)[i] = paste(u,
                                           collapse = '.',
                                           sep = '.')
      }
      ############# step 3
      names(AllCombinations2)[i] = gsub(
        pattern = 'Genotype.',
        replacement = '',
        x = names(AllCombinations2)[i]
      )
      ############# step 4
      names(AllCombinations2)[i] = toupper(names(AllCombinations2)[i])
    }
    ############# step 5
    AllCombinations2 = AllCombinations2[grepl(
      pattern = 'MALE|FEMALE|OVERALL',
      x = names(AllCombinations2),
      ignore.case = TRUE
    )]
    ############# step 6
    if (length(AllCombinations2) < 1)
      return(NULL)
    nam = names(AllCombinations2)
    for (i in 0:1000) {
      sb = substr(nam, start = nchar(nam) - i, nchar(nam))
      for (j in 1:length(sb)) {
        if (numbers_only(sb[j]) && i < nchar(sb[j]))
          nam[j] = substr(nam[j], start = 0, stop = nchar(nam[j]) - i -
                            1)
      }
    }
    ############# step 7
    gp7 = grepl(pattern = '.MPTERM',
                x = nam,
                fixed = TRUE)
    nam[!gp7] = paste0(nam[!gp7], '.MPTERM')

    ############# step 8
    if (!is.null(parameter_stable_id)) {
      controlCat = read.delim('CategoryRemapping.tsv', sep = '\t')
      conCat = controlCat[controlCat$IMPRESS.ID %in% parameter_stable_id,]
      control = conCat[conCat$CATEGORY.1 == 0, 2]
      if (length(control) > 0) {
        for (c in control) {
          nam =   gsub(
            pattern = paste0('.', toupper(c), '.'),
            replacement = '.',
            x = nam,
            fixed = TRUE
          )
        }
        nam = gsub(
          pattern = '..',
          replacement = '.',
          x = nam,
          fixed = TRUE
        )
      }
    }
    ############# Finally!
    names(AllCombinations2) = nam
    tag = AllCombinations2

  } else if (method %in% 'RR') {
    ###########################################################################
    fmodels = obj$`Additional information`$Analysis$`Further models`
    if (is.null(fmodels))
      return(NULL)

    AllCombinations = lapply(fmodels, function(x) {
      lapply(x, function(y) {
        DirectionTagFE(x = y$Result$p.value, threshold = threshold)
      })
    })

    if (is.null(AllCombinations))
      return(NULL)
    #### Make the list as sequence of names attached with dot (.)
    AllCombinations1 = unlist(AllCombinations)
    ####
    AllCombinations1 = AllCombinations1[!(grepl(
      pattern = 'Low.',
      x = names(AllCombinations1),
      fixed = TRUE
    ) & AllCombinations1 %in% 'INCREASED')]
    AllCombinations1 = AllCombinations1[!(grepl(
      pattern = 'High.',
      x = names(AllCombinations1),
      fixed = TRUE
    ) & AllCombinations1 %in% 'DECREASED')]

    #### Keep only genotype analysis
    AllCombinations2 = AllCombinations1[grepl(pattern = 'Genotype.Genotype|Genotype.Genotype_Male|Genotype.Genotype_Female', x =
                                                names(AllCombinations1))]
    #### The abnormal case is made from the data
    for (i in seq_along(AllCombinations2)) {
      names(AllCombinations2)[i]  = gsub(
        pattern = 'data_point.Genotype.Genotype_Male',
        paste0('MALE.', AllCombinations2[i], '.OVERALL.MPTERM'),
        x = names(AllCombinations2)[i]
      )
      names(AllCombinations2)[i]  = gsub(
        pattern = 'data_point.Genotype.Genotype_Female',
        paste0('FEMALE.', AllCombinations2[i], '.OVERALL.MPTERM'),
        x = names(AllCombinations2)[i]
      )
      names(AllCombinations2)[i]  = gsub(
        pattern = 'data_point.Genotype.Genotype',
        paste0('UNSPECIFIED.', AllCombinations2[i], '.OVERALL.MPTERM'),
        x = names(AllCombinations2)[i]
      )
      names(AllCombinations2)[i] = gsub(
        pattern = '..',
        replacement = '.',
        x = names(AllCombinations2)[i],
        fixed = TRUE
      )
      ############# step 4
      names(AllCombinations2)[i] = toupper(names(AllCombinations2)[i])
    }
    ############# step 6
    if (length(AllCombinations2) < 1)
      return(NULL)
    nam = names(AllCombinations2)
    for (i in 0:1000) {
      sb = substr(nam, start = nchar(nam) - i, nchar(nam))
      for (j in 1:length(sb)) {
        if (numbers_only(sb[j]) && i < nchar(sb[j]))
          nam[j] = substr(nam[j], start = 0, stop = nchar(nam[j]) - i -
                            1)
      }
    }

    ############# Finally!
    names(AllCombinations2) = nam #gsub(pattern = 'LOW.|HIGH.',replacement = '',x = nam)
    tag = AllCombinations2
  } else{
    tag = NULL
    write(
      x = paste(
        head(expDetailsForErrorOnly, 18),
        sep = '\t',
        collapse = '\t'
      ),
      file = 'ErrorneousCases.tsv.err',
      ncolumns = 5000
    )
  }
  return(tag)
}


unScrewProcedure = function(x) {
  r = unlist(strsplit(x = x, split = '~'))
  return(r)
}


merge.two = function(l1, l2) {
  if (is.null(l1) && is.null(l2))
    return(NULL)
  if (is.null(l1) && !is.null(l2))
    return(l1)
  if (!is.null(l1) && is.null(l2))
    return(l2)
  #l = c(l1, l2[!names(l2) %in% names(l1)])
  l = c(l1, l2)
  return(l)
}
###########################

removeAbnormalIfIncDecDetected = function(x = NULL,
                                          method = 'MM',
                                          active = 1) {
  if (!active)
    return(x)
  if (length(x) > 0 &&
      !method %in% 'RR' &&
      any(grepl(
        x = names(x),
        pattern = ('\\.INCREASED\\.|\\.DECREASED\\.')
      ))) {
    message('\t INCREASED/DECREASED detected. ABNORMAL will be ignored')
    x = x[!grepl(x = names(x),
                 pattern = ('.ABNORMAL.'),
                 fixed = TRUE)]
  }
  return(x)
}

MatchTheRestHalfWithTheFirstOne = function(x) {
  if (length(x) < 1)
    return(x)
  #x = x[!is.na(x)]
  mid = length(x) / 2
  x[(mid + 1):(2 * mid)] = x[1:mid]
  return(x)
}

DecIncDetector = function(x) {
  if (length(x) < 1)
    return(NA)
  r = c()
  if (any(grepl(pattern = 'DECREASED', x = names(x))))
    r = c(r, 'DECREASED')
  if (any(grepl(pattern = 'INCREASED', x = names(x))))
    r = c(r, 'INCREASED')
  if (any(grepl(pattern = 'ABNORMAL', x = names(x))))
    r = c(r, 'ABNORMAL')
  if (any(grepl(pattern = 'INFERRED', x = names(x))))
    r = c(r, 'INFERRED')

  if (length(r) < 1)
    r = NA
  if (length(r) > 1) {
    if (any(grepl(pattern = 'ABNORMAL', x = r)))
      r = 'ABNORMAL'
    if (any(grepl(pattern = 'INFERRED', x = r)))
      r = 'INFERRED'
  }
  return(r)
}

DecIncDetectorRR = function(x) {
  if (length(x) < 1)
    return(NA)
  r = c()
  ###################################### MALE
  ###################################### LOW
  r1 = multiGrepl(pattern = c('(DECREASED)|(INCREASED)', 'LOW', 'MALE'),
                  x = names(x))
  if (sum(r1) > 1) {
    rlowM = c('ABNORMAL')
  } else if (sum(r1) == 1) {
    if (grepl(pattern = 'INCREASED', x = names(x[r1])))
      rlowM = c('INCREASED')
    if (grepl(pattern = 'DECREASED', x = names(x[r1])))
      rlowM = c('DECREASED')
  } else{
    rlowM = NA
  }
  ###################################### High
  r2 = multiGrepl(pattern = c('(DECREASED)|(INCREASED)', 'HIGH', 'MALE'),
                  x = names(x))
  if (sum(r2) > 1) {
    rhighM = c('ABNORMAL')
  } else if (sum(r2) == 1) {
    if (grepl(pattern = 'INCREASED', x = names(x[r2])))
      rhighM = c('INCREASED')
    if (grepl(pattern = 'DECREASED', x = names(x[r2])))
      rhighM = c('DECREASED')
  } else{
    rhighM = NA
  }
  ######################################
  r3 = na.omit(c(rlowM, rhighM))
  if (length(r3) > 0) {
    if (length(r3) == 1) {
      r = c(r, na.omit(r3))
    } else{
      r = c(r, 'ABNORMAL')
    }
  }
  ###################################### FEMALE
  ###################################### LOW
  r4 = multiGrepl(pattern = c('(DECREASED)|(INCREASED)', 'LOW', 'FEMALE'),
                  x = names(x))
  if (sum(r4) > 1) {
    rlowF = c('ABNORMAL')
  } else if (sum(r4) == 1) {
    if (grepl(pattern = 'INCREASED', x = names(x[r4])))
      rlowF = c('INCREASED')
    if (grepl(pattern = 'DECREASED', x = names(x[r4])))
      rlowF = c('DECREASED')
  } else{
    rlowF = NA
  }
  ######################################
  r5 = multiGrepl(
    pattern = c('(DECREASED)|(INCREASED)', 'HIGH', 'FEMALE'),
    x = names(x)
  )
  if (sum(r5) > 1) {
    rhighF = c('ABNORMAL')
  } else if (sum(r5) == 1) {
    if (grepl(pattern = 'INCREASED', x = names(x[r5])))
      rhighF = c('INCREASED')
    if (grepl(pattern = 'DECREASED', x = names(x[r5])))
      rhighF = c('DECREASED')
  } else{
    rhighF = NA
  }
  ######################################
  r6 = na.omit(c(rlowF, rhighF))
  if (length(r6) > 0) {
    if (length(r6) == 1) {
      r = c(r, na.omit(r6))
    } else{
      r = c(r, 'ABNORMAL')
    }
  }
  ######################################
  ###################################### ABNORMAL
  r7 = multiGrepl(
    pattern = c('(DECREASED)|(INCREASED)', 'LOW', 'ABNORMAL'),
    x = names(x)
  )
  if (sum(r7) > 1) {
    rlowA = c('ABNORMAL')
  } else if (sum(r7) == 1) {
    if (grepl(pattern = 'INCREASED', x = names(x[r7])))
      rlowA = c('INCREASED')
    if (grepl(pattern = 'DECREASED', x = names(x[r7])))
      rlowA = c('DECREASED')
  } else{
    rlowA = NA
  }
  ######################################
  r8 = multiGrepl(
    pattern = c('(DECREASED)|(INCREASED)', 'HIGH', 'ABNORMAL'),
    x = names(x)
  )
  if (sum(r8) > 1) {
    rhighA = c('ABNORMAL')
  } else if (sum(r8) == 1) {
    if (grepl(pattern = 'INCREASED', x = names(x[r8])))
      rhighA = c('INCREASED')
    if (grepl(pattern = 'DECREASED', x = names(x[r8])))
      rhighA = c('DECREASED')
  } else{
    rhighA = NA
  }
  ######################################
  r9 = na.omit(c(rlowA, rhighA))
  if (length(r9) > 0) {
    if (length(r9) == 1) {
      r = c(r, na.omit(r9))
    } else{
      r = c(r, 'ABNORMAL')
    }
  }
  ######################################
  abpattern = inferpattern = FALSE

  if (length(x) > 0) {
    abpattern = grepl(pattern = 'ABNORMAL', x = names(x))
    inferpattern = grepl(pattern = 'INFERRED', x = names(x))
  }

  if (length(r) < 1 && !(any(abpattern) || any(inferpattern)))
    r = NA
  ######################################
  if (length(r) < 1 && (any(abpattern) || any(inferpattern)))
    if (any(grepl(pattern = 'ABNORMAL', x = names(x))))
      r = 'ABNORMAL'
  if (any(grepl(pattern = 'INFERRED', x = names(x))))
    r = 'INFERRED'
  ######################################
  if (length(r) > 1) {
    if (any(grepl(pattern = 'ABNORMAL', x = r)))
      r = 'ABNORMAL'
    if (any(grepl(pattern = 'INFERRED', x = r)))
      r = 'INFERRED'
  }
  ######################################
  return(r)
}

detectLevel = function(x, level) {
  if (is.null(x) || length(x) < 1)
    return(FALSE)
  r = grepl(pattern = level, x = x)
  return(any(r))
}

NullOrvalueReturn = function(x, list) {
  if (length(x) < 1)
    return(NULL)
  else
    return(list)
}

fA = function(x, pasteterms = TRUE) {
  if (length(x) > 0)
    x = x[!duplicated(x)]
  if (length(x) > 1) {
    x2 = x[!grepl(pattern = c('(MALE)|(FEMALE)'), names(x)) &
             grepl(pattern = c('ABNORMAL'), names(x))]
    x3 = x[grepl(pattern = c('ABNORMAL'), names(x))]
    if (length(x2) > 1 && length(x3) > 0 && pasteterms)
      x = paste(x3, sep = '~', collapse = '~')
    else if (length(x2) > 1 && length(x3) < 1 && pasteterms)
      x = paste(x2, sep = '~', collapse = '~')
    else if (pasteterms)
      x = paste(x, sep = '~', collapse = '~')
    else
      x = x
  }
  return(x)
}

fM = function(x, pasteterms = TRUE) {
  if (length(x) > 0)
    x = x[!duplicated(x)]
  if (length(x) > 1) {
    x2 = x[!grepl(pattern = c('ABNORMAL'), names(x)) &
             grepl(pattern = c('MALE'), names(x))]
    x3 = x[grepl(pattern = c('ABNORMAL'), names(x))]
    if (length(x2) > 1 && length(x3) > 0 && pasteterms)
      x = paste(x3, sep = '~', collapse = '~')
    else if (length(x2) > 1 && length(x3) < 1 && pasteterms)
      x = paste(x2, sep = '~', collapse = '~')
    else if (pasteterms)
      x = paste(x, sep = '~', collapse = '~')
    else
      x = x
  }
  return(x)
}

fF = function(x, pasteterms = TRUE) {
  if (length(x) > 0)
    x = x[!duplicated(x)]
  if (length(x) > 1) {
    x2 = x[!grepl(pattern = c('ABNORMAL'), names(x)) &
             grepl(pattern = c('FEMALE'), names(x))]
    x3 = x[grepl(pattern = c('ABNORMAL'), names(x))]
    if (length(x2) > 1 && length(x3) > 0 && pasteterms)
      x = paste(x3, sep = '~', collapse = '~')
    else if (length(x2) > 1 && length(x3) < 1 && pasteterms)
      x = paste(x2, sep = '~', collapse = '~')
    else if (pasteterms)
      x = paste(x, sep = '~', collapse = '~')
    else
      x = x
  }
  return(x)
}

bselect = function(x) {
  if(length(x)<1)
    return(x)
  nx = names(x)
  if (length(nx) > 1 && any(grepl(pattern = 'ABNORMAL', x = nx)))
    x = x[grepl(pattern = 'ABNORMAL', x = nx)][1]
  if (length(nx) > 1 && any(grepl(pattern = 'INFERRED', x = nx)))
    x = x[grepl(pattern = 'INFERRED', x = nx)][1]
  return(x)
}

MaleFemaleAbnormalCategories = function(x, method = 'AA') {
  fgrep = grepl(pattern = 'FEMALE', names(x), fixed = TRUE)
  mgrep = grepl(pattern = 'MALE', names(x), fixed = TRUE) & !fgrep
  agrep = grepl(pattern = '(ABNORMAL)|(INFERRED)|(OVERAL)', names(x)) &
    !fgrep & !mgrep

  if (method %in% 'RR') {
    fevent = DecIncDetectorRR(x[fgrep])
    mevent = DecIncDetectorRR(x[mgrep])
    oevent = DecIncDetectorRR(x[agrep])
  } else{
    fevent = DecIncDetector(x[fgrep])
    mevent = DecIncDetector(x[mgrep])
    oevent = DecIncDetector(x[agrep])
  }
  # if males and females are in cross direction set abnormal otherwise select male/females
  fmevents = c(fevent, mevent)
  if (all(!is.na(fmevents)) && fmevents[1] == fmevents[2])
    oevent = mevent
  if (is.na(fmevents[1]) && !is.na(fmevents[2]))
    oevent = fmevents[2]
  if (!is.na(fmevents[1]) && is.na(fmevents[2]))
    oevent = fmevents[1]
  #############################

  #############################
  if (method %in% 'RR') {
    MPTERM = list(
      NullOrvalueReturn(
        x[agrep],
        list(
          'term_id' = ifelse(length(fA(x[agrep], pasteterms = FALSE)) >
                               1, fA(x[agrep], pasteterms = FALSE)[1], fA(x[agrep])),
          event = oevent,
          sex = "not_considered",
          otherPossibilities = ifelse(length(fA(x[agrep], pasteterms = FALSE)) >
                                        1, fA(x[agrep]), '')
        )
      ),
      NullOrvalueReturn(
        x[fgrep],
        list(
          'term_id' = ifelse(length(fF(x[fgrep], pasteterms = FALSE)) >
                               1, fF(x[fgrep], pasteterms = FALSE)[1], fF(x[fgrep])),
          event = fevent,
          sex = "female",
          otherPossibilities = ifelse(length(fF(x[fgrep], pasteterms = FALSE)) >
                                        1, fF(x[fgrep]), '')
        )
      ),
      NullOrvalueReturn(
        x[mgrep],
        list(
          'term_id' = ifelse(length(fM(x[mgrep], pasteterms = FALSE)) >
                               1, fM(x[mgrep], pasteterms = FALSE)[1], fM(x[mgrep])),
          event = mevent,
          sex = "male",
          otherPossibilities = ifelse(length(fM(x[mgrep], pasteterms = FALSE)) >
                                        1, fM(x[mgrep]), '')
        )
      )
    )
  } else{
    MPTERM = list(
      NullOrvalueReturn(
        x[agrep],
        list(
          'term_id' = bselect(x[agrep]),
          event = oevent,
          sex = "not_considered"
        )
      ),
      NullOrvalueReturn(x[fgrep], list(
        'term_id' = bselect(x[fgrep]),
        event = fevent,
        sex = "female"
      )),
      NullOrvalueReturn(x[mgrep], list(
        'term_id' = bselect(x[mgrep]),
        event = mevent,
        sex = "male"
      ))
    )
  }
  return(MPTERM)
}
##############################
f = function(statpacket = NULL,
             level = 10 ^ -4) {
  ulistTag3 = MPTERMS = NA
  message('Running the annotation pipeline')
  if (!statpacket$V2 %in% 'Successful') {
    message('Not a successfull StatPackage!')
    return(invisible(list(
      MPTERM = ulistTag3, statpacket = statpacket
    )))
  }
  pipeline = statpacket$V15
  procedure = statpacket$V3
  parameter = statpacket$V6
  json      =  jsonlite::fromJSON(statpacket$V20)
  method   = GetMethodStPa(json$Result$`Vector output`$`Normal result`$`Applied method`)
  ##################################################################
  Gtag = GenotypeTag(
    obj = json$Result$`Vector output`$`Normal result`,
    parameter_stable_id = statpacket$V6,
    threshold = level
  )
  ##################################################################
  message('\t Reading the index file ...')
  load('mp_chooser_20200520.json.Rdata')
  message(
    '\t~>',
    paste(
      pipeline,
      procedure,
      parameter,
      sep = '~>',
      collapse = '~>'
    ),
    '\n\t~>',
    json$Result$Details$`Gene page URL`
  )

  ################################
  b = a[[unScrewProcedure(pipeline)[1]]]
  c = b[[unScrewProcedure(procedure)[1]]]
  d = c[[unScrewProcedure(parameter)[1]]]
  if (is.null(d)) {
    message('No annotation available by IMPC. See https://www.mousephenotype.org/impress/')
    return(invisible(list(
      MPTERM = ulistTag3, statpacket = statpacket
    )))
  }
  ################################
  if (length(d) > 0 &&
      length(Gtag) > 0) {
    ################################
    ulistTag  = unlist(Gtag)
    ulistD    = unlist(d)
    names(ulistD)  = toupper(names(ulistD))
    ulistD         = removeAbnormalIfIncDecDetected(ulistD, method = method, active = FALSE)
    if (length(ulistD) < 1) {
      message('No annotation available by IMPC. See https://www.mousephenotype.org/impress/')
      return(invisible(list(
        MPTERM = ulistTag3, statpacket = statpacket
      )))
    }
    ################################
    ulistTag2 = ulistTag
    names(ulistTag2) = gsub(
      pattern = 'MALE|FEMALE',
      x = names(ulistTag2),
      replacement = 'UNSPECIFIED'
    )

    ulistTag3 = merge.two(ulistTag2, ulistTag)
    for (name in names(ulistTag3)) {
      # print(name)
      splN = unlist(strsplit(name, split = '.', fixed = TRUE))
      splN = splN[!splN %in% c('LOW', 'HIGH')]
      splnI = multiGrepl(pattern = splN,
                         x = names(ulistD),
                         fixed = TRUE)
      ulistTag3[names(ulistTag3) %in% name] = ifelse(is.null(ulistD[splnI]), 'CanNotFindMPTerm', head(ulistD[splnI], 1))
    }
    #print(ulistD)
    ################################
    if (length(ulistTag3) < 1)
      return(invisible(list(
        MPTERM = ulistTag3, statpacket = statpacket
      )))
    ##########################
    ulistTag3 = MatchTheRestHalfWithTheFirstOne(ulistTag3)
    ulistTag3 =  ulistTag3[!duplicated(names(ulistTag3))]
    ##########################
    if (length(ulistTag3) < 1)
      return(invisible(list(
        MPTERM = ulistTag3, statpacket = statpacket
      )))
    ##########################
    ulistTag3 = ulistTag3[!is.na(ulistTag3)]
    ################################
    MPTERMS = MaleFemaleAbnormalCategories(x = ulistTag3, method = method)
    MPTERMS = rlist::list.clean(MPTERMS)
    ################################
    if (length(ulistTag3) > 0 &&
        method %in% 'RR') {
      lIncDec = multiGrepl(pattern = c('INCREASED|DECREASED'),
                           x = names(ulistTag3))
      lAbn    = multiGrepl(pattern = c('ABNORMAL'),
                           x = names(ulistTag3))
      if (sum(lIncDec) >= 2 && sum(lAbn) > 0) {
        ulistTag3 = ulistTag3[lAbn][1]
        names(ulistTag3) = gsub(pattern = '(LOW\\.)|(\\.HIGH)',
                                replacement = '',
                                names(ulistTag3))
      }
    }
    # Do not use unique!
    if (length(ulistTag3) > 0) {
      ulistTag3 =  ulistTag3[!duplicated(names(ulistTag3))]
    }
  }
  if (!is.null(MPTERMS)     &&
      length(ulistTag3) > 0 &&
      length(na.omit(ulistTag3)) > 0) {
    json$Result$Details$MPTERM = MPTERMS
    statpacket$V20 = DRrequiredAgeing:::FinalJson2ObjectCreator(FinalList = json)
  }
  return(invisible(list(
    MPTERM = MPTERMS, statpacket = statpacket
  )))
}

#################################################################################
Write2Postg = function(df) {
  library(RPostgreSQL)
  library(data.table)
  drv <- dbDriver("PostgreSQL")
  con <- dbConnect(
    drv,
    dbname = "test",
    host = "hh-yoda-07-07",
    port = 5432,
    user = "impc",
    password = 'impc'
  )

  dbBegin(conn = con)
  r = dbWriteTable(
    conn = con,
    name = "dr12withmptermsv719082020_16",
    value = df,
    append = TRUE,
    row.names = FALSE
  )
  dbCommit(conn = con)
  dbDisconnect (conn = con)
  return(r)
}

randomIdGenerator = function(l = 10) {
  r=paste0(sample(
    x = 0:9,
    size = l,
    replace = TRUE
  ), collapse = '')
  return(r)
}

# file = 'https://www.ebi.ac.uk/~hamedhm/windowing/DR11/jobs/Results_DR11V2OpenStats/RBRC/IMPC_CSD/IMPC_CSD_050_002/RIKEN_Fbn2_A06_/homozygote/978482d700f332623fc65f3e729dea93/output_Successful.tsv'
# file ='https://www.ebi.ac.uk/~hamedhm/windowing/DR11/jobs/Results_DR11V2OpenStats/ICS/ESLIM_008/ESLIM_008_001_014/EPD0060_2_H09/heterozygote/d41d8cd98f00b204e9800998ecf8427e/output_Successful.tsv'
# file = 'https://www.ebi.ac.uk/~hamedhm/windowing/DR12/jobs/Results_DR12V1OpenStats/ICS/ESLIM_008/ESLIM_008_001_014/EPD0060_2_H09/heterozygote/d41d8cd98f00b204e9800998ecf8427e/output_Successful.tsv'
# download.file(url = file, destfile = 'delme.tsv')
# file = 'delme.tsv'
###########################################
flist = readLines(con = file[1])
lflist = length(flist)
id = 1
rnd <-
  DRrequired:::RandomRegardSeed(
    n = 1,
    decimal = 0,
    stringOutput = 1,
    max = 2500,
    round = 1
  )
for (i in 1:lflist) {
  cat('\r', i, '/', lflist)
  status = FALSE
  file = flist[i]
  cat("\n", i, "/", lflist, "~>", file, "")
  if (file.exists(file) &&
      (grepl(pattern = 'NotProcessed', x = file) ||
       grepl(pattern = 'Successful'  , x = file))) {
    df = fread(
      file = file,
      header = FALSE,
      sep = '\t',
      quote = "",
      stringsAsFactors = FALSE
    )
    if (ncol(df) != 20)
      next
    ###################
    r = f(statpacket = df, level = .0001)
    #write(paste(r$statpacket[1,19],sep = '\t',collapse = '\t'),file = 'd:/delme.tsv',ncolumns = 10^4)
    df = r$statpacket
    ###################
    df = cbind(as.numeric(paste0(randomIdGenerator(),format(Sys.time(), "%H%S"))), df)
    names(df)  =  c(
      "id",
      "analysisTime",
      "status",
      "procedure_group",
      "procedure_stable_id",
      "procedure_name",
      "parameter_stable_id",
      "parameter_name",
      "phenotyping_center",
      "allele_symbol",
      "allele_name",
      "allele_accession_id",
      "gene_symbol",
      "gene_accession_id",
      "pipeline_name",
      "pipeline_stable_id",
      "strain_accession_id",
      "metadata_group",
      "zygosity",
      "colony_id",
      "statpacket"
    )
    status  = Write2Postg(df = df[1, ])
    id = id + 1
    #print(file)
  }
  if (!dir.exists("log")) {
    dir.create("log")
  }
  write(
    paste(status, file, sep = "\t-->", collapse = "\t-->"),
    file = paste0("./log/", rnd, "__Log.log"),
    append = TRUE
  )
  gc()
}

############################## batch maker
# lf = list.files(path = getwd(),pattern = 'split_index',full.names = TRUE,recursive = FALSE,include.dirs = FALSE)
# err = paste0(' -e ',dirname(lf),'/err/',basename(lf))
# out = paste0(' -o ',dirname(lf),'/out/',basename(lf))
# batch = paste0('bsub -M 12000 ',err,out,' Rscript loader.R "',basename(lf),'"')
# write(batch,'jobs.bch')
