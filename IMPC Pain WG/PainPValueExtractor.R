args = commandArgs(trailingOnly = TRUE)
##########################################
VariableNames = function(x = NULL) {
  colNames = c(
    "Centre",
    "Gene_symbol",
    "Colony_id",
    "Zygosity",
    "Control strain",
    "Mutant strain",
    "Metadata_group",
    "Age_in_weeks (w1~w2~w3)",
    "Test0 vs Test1 pvalue",
    "Test0 vs Test2 pvalue",
    "Test1 vs Test2 pvalue",
    "Applied method",
    "Gender included in analysis",
    "Residual variances homogeneity",
    "Batch included",
    "Sexual dymorphism detected?",
    "Genotype p-value",
    "Genotype estimate",
    "Overal - Genotype effect size",
    "Time0h - Genotype effect size",
    "Time24h - Genotype effect size",
    "Time48h - Genotype effect size",
    "Type of Genotype effect size",
    "Sex p-value",
    "Sex estimate",
    "Overal - Sex effect size",
    "Time0- Sex effect size",
    "Time24h - Sex effect size",
    "Time48h - Sex effect size",
    "Type of Sex effect size",
    "Sex FvKO p-value",
    "Sex FvKO estimate",
    "Overal - Sex FvKO effect size",
    "Time0h - Sex FvKO effect size",
    "Time24h - Sex FvKO effect size",
    "Time48h - Sex FvKO effect size",
    "Type of Sex FvKO effect size",
    "Sex MvKO p-value",
    "Sex MvKO estimate",
    "Overal - Sex MvKO effect size",
    "Time0h - Sex MvKO effect size",
    "Time24h - Sex MvKO effect size",
    "Time48h- Sex MvKO effect size",
    "Type of Sex MvKO effect size",
    "Weight p-value",
    "Weight estimate",
    "Weight standard error",
    "Weight effect size",
    "Type of Weigh effect size",
    "Summary statistics"
  )
  return(colNames)
}

qvalue0 = function(x, ...) {
  library("qvalue")
  qx = rep(NA, length(x))
  if (length(x[!is.na(x)]) > 0) {
    qx[!is.na(x)] = p.adjust(x[!is.na(x)], method = 'fdr')
  } else{
    message('Less than 1 pvalue!')
  }
  return(qx)
}

unlist0 = function(x, active = TRUE) {
  x2 = unlist(x)
  if ((is.null(x) || is.null(x2)) && active) {
    return(NA)
  } else{
    return(x2)
  }
}
is.nullorNA = function(x) {
  if (is.null(x) || all(x == '')) {
    return(NA)
  } else{
    return(x)
  }
}


SelectAnalysis = function(object) {
  summStats = lapply(object$`Additional information`$Data$`Summary statistics`, function(x) {
    x$Count
  })
  r = c(
    ######
    paste(
      unlist0(object$'Applied method'),
      paste0(
        ' Final (optimised) model: ',
        unlist0(
          object$`Additional information`$Analysis$`Model setting`$`Final model`$Formula
        )
      ),
      collapse = '~',
      sep = '~'
    ),
    unlist0(
      object$'Additional information'$Analysis$'Gender included in analysis'
    ) ,
    ######
    unlist0(object$'Residual variances homogeneity')  ,
    unlist0(object$'Batch included')                  ,
    OpenStats:::pasteComma(unlist0(object$'Genotype contribution'$'Sexual dimorphism detected')),
    # pvalue/Estimates/Effect sizes
    unlist0(object$'Genotype p-value')                ,
    unlist0(object$'Genotype estimate'$Value)         ,
    unlist0(object$'Genotype effect size'$Value)      ,
    unlist0(
      object$`Additional information`$Analysis$`Effect sizes`$Genotype_Test0$Value
    )  ,
    unlist0(
      object$`Additional information`$Analysis$`Effect sizes`$Genotype_Test1$Value
    ) ,
    unlist0(
      object$`Additional information`$Analysis$`Effect sizes`$Genotype_Test2$Value
    ) ,
    unlist0(object$'Genotype effect size'$Type)       ,
    #
    unlist0(object$'Sex p-value')                     ,
    unlist0(object$'Sex estimate'$Value)              ,
    unlist0(object$'Sex effect size'$Value)           ,
    unlist0(
      object$`Additional information`$Analysis$`Effect sizes`$Sex_Test0$Value
    )  ,
    unlist0(
      object$`Additional information`$Analysis$`Effect sizes`$Sex_Test1$Value
    ) ,
    unlist0(
      object$`Additional information`$Analysis$`Effect sizes`$Sex_Test2$Value
    ) ,
    unlist0(object$'Sex effect size'$Type)            ,
    #
    unlist0(object$'Sex FvKO p-value')               ,
    unlist0(object$'Sex FvKO estimate'$Value)         ,
    unlist0(object$'Sex FvKO effect size'$Value)      ,
    unlist0(
      object$`Additional information`$Analysis$`Effect sizes`$`Genotype_Female Test0`$Value
    )   ,
    unlist0(
      object$`Additional information`$Analysis$`Effect sizes`$`Genotype_Female Test1`$Value
    )  ,
    unlist0(
      object$`Additional information`$Analysis$`Effect sizes`$`Genotype_Female Test2`$Value
    )  ,
    unlist0(object$'Sex FvKO effect size'$Type)       ,
    #
    unlist0(object$'Sex MvKO p-value')                ,
    unlist0(object$'Sex MvKO estimate'$Value)         ,
    unlist0(object$'Sex MvKO effect size'$Value)      ,
    unlist0(
      object$`Additional information`$Analysis$`Effect sizes`$`Genotype_Male Test0`$Value
    )  ,
    unlist0(
      object$`Additional information`$Analysis$`Effect sizes`$`Genotype_Male Test1`$Value
    )  ,
    unlist0(
      object$`Additional information`$Analysis$`Effect sizes`$`Genotype_Male Test2`$Value
    )  ,
    unlist0(object$'Sex MvKO effect size'$Type)       ,
    # Bodyweight
    unlist0(object$'Weight p-value')          ,
    unlist0(object$'Weight estimate'$Value)   ,
    unlist0(object$'Weight standard error')   ,
    unlist0(object$'Weight effect size'$Value)   ,
    unlist0(object$'Weight effect size'$Type)   ,
    # Others
    paste(
      names(summStats),
      summStats,
      sep = ': ',
      collapse = ' | '
    )
  )
  return(r)
}


createValues = function(r0, r1) {
  x = c(
    unlist(r0[1, 1:11])              ,
    ################# VectorOutput Results
    SelectAnalysis(r1)
  )
  return(x)
}
########## Main function
f2 = function(file,ofname = 'values2.tsv') {
  library(RJSONIO)
  library(DRrequired)
  library(data.table)
  library(DBI)
  if (file.exists(paste0('./resultF/', ofname)))
    unlink(paste0('./resultF/', ofname))
  ####
  rs = read.delim(
    file,
    header = FALSE,
    sep = '\t',
    check.names = FALSE,
    quote = '',
    stringsAsFactors = FALSE
  )
  rs = rs[!apply(rs == "", 1, all),]
  n  = nrow(rs)
  ####
  for(i in 1:n){
    r0 = rs[i,]
    r1 =
      tryCatch(
        expr =
          jsonlite::fromJSON(r0[, 12], flatten = TRUE),
        warning = function(w) {
          write(file,
                file = paste0('Error_inJSON_', Sys.Date(), '.txt'),
                append = TRUE)
          return(NULL)
        },
        error = function(e) {
          write(file,
                file = paste0('Error_inJSON_', Sys.Date(), '.txt'),
                append = TRUE)
          return(NULL)
        }
      )
    ######
    if (is.null(r1))
      next
    ###### only MM's
    #if (!is.null(method) && method  %in% 'MM') {
    message(paste(i, '|',  n ,'~>',round(i/n*100),'%'))
    x = createValues(r0, r1)
    fpath = paste0('./resultF/', r1$Result$Details$`Response type`, '/', ofname)
    write(
      x = paste(x, collapse = '\t'),
      file = DRrequiredAgeing:::file.path0(
        fpath,
        create = TRUE,
        check  = FALSE,
        IncludedFileName = TRUE
      ),
      append = TRUE,
      ncolumns = 10 ^ 4
    )
  }
  dfNotCorrected = read.delim(
    file   = fpath,
    header = FALSE,
    sep    = '\t',
    check.names = FALSE
  )
  names(dfNotCorrected) = VariableNames()
  table(dfNotCorrected$Centre)
  counter = 1
  dfCorrected = NULL
  for (c in unique(dfNotCorrected$Centre)) {
    for (cs in unique(dfNotCorrected$`Control strain`)) {
      message(' ~> ', c, ' - ', cs)
      sset = subset(x = dfNotCorrected ,
                    dfNotCorrected$Centre             %in% c &
                      dfNotCorrected$`Control strain` %in% cs)
      sset$`Genotype p-value corrected` = qvalue0(sset$`Genotype p-value`)
      sset$`Sex p-value corrected`      = qvalue0(sset$`Sex p-value`)
      sset$`Sex FvKO p-value corrected` = qvalue0(sset$`Sex FvKO p-value`)
      sset$`Sex MvKO p-value corrected` = qvalue0(sset$`Sex MvKO p-value`)
      sset$`Age included in analysis?`  = grepl(pattern = 'age',x = sset$`Applied method`)
      if (counter < 2) {
        dfCorrected = sset
        counter = counter + 1
      } else{
        dfCorrected = rbind(dfCorrected, sset)
      }
    }
  }
  write.table(
    x = dfCorrected,
    file = paste0(fpath, '_corrected.tsv'),
    sep = '\t',
    row.names = FALSE
  )


}

f2(file = args[1],ofname = args[2])
