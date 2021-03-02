args = commandArgs(trailingOnly = TRUE)
library(data.table)
library(jsonlite)
library(digest)
####################
unlist0 = function(x, active = TRUE) {
  x2 = unlist(x)
  if ((is.null(x) || is.null(x2)) && active) {
    return(NA)
  } else{
    return(paste(x2, sep = '~', collapse = '~'))
  }
}

is.nullorNA = function(x) {
  if (is.null(x) || all(x == '')) {
    return(NA)
  } else{
    return(x)
  }
}

SelectOthers = function(object) {
  paste0(is.nullorNA(c(
    unlist(object$`Concurrent control selection`  )
  )))
}

SelectMISC = function(object) {
  ReferenceGene = c('Nxn', 'Rnf10', 'Ap4e1', 'Prkab1', 'Dnase1l2', 'Dbn1')
  rReference = object$V12   %in% ReferenceGene
  rIgnorome  = object$V11  %in% DRrequiredAgeing:::ignoromeGenes()
  rBehaviour = object$V6   %in% DRrequiredAgeing:::BehviourParamters()
  return(c(rReference, rIgnorome, rBehaviour))
}

URLdecode0 = function(x) {
  if (is.null(x))
    return(NA)
  return(URLdecode(URL = x))
}

MakeURL = function(object, objectJSON) {
  r = URLencode(
    paste0(
      "https://wwwdev.ebi.ac.uk/mi/impc/dev/phenotype-archive/media/images/windowing/?alleleSymbol=",
      object$V9,
      "&colonyID=",
      object$V19,
      "&procedure= | ",
      object$V3,
      "&parameter= | ",
      object$V6,
      "&center=",
      object$V8,
      "&zygosity=",
      object$V18,
      "&metadata=",
      object$V17
    )
  )
  r = c(r, URLdecode0(unlist(
    objectJSON$Result$Details$`Exported raw data file`
  )))
  return(r)
}

length0 = function(x) {
  r = length(x)
  if (r < 1)
    return(NA)
  return(r)
}

SelectWindowingParameters = function(object) {
  r = c(
    ###### Windowing parameters
    unlist0(object$`Window parameters`$l[1])             ,
    unlist0(object$`Window parameters`$k[1])             ,
    unlist0(object$`Window parameters`$`Min obs required in the window`),
    unlist0(object$`Window parameters`$`The number of DOE in the window`) ,
    unlist0(object$`Window parameters`$Threshold)                   ,
    unlist0(length0(object$`Window parameters`$`Window weights`))             ,
    unlist0(object$`Window parameters`$`Total obs or weight in the window`)
  )
  return(r)
}

SelectAnalysis = function(object) {
  r = c(
    ######
    unlist0(object$`Applied method`)                ,
    unlist0(object$`Classification tag`$`Classification tag`) ,
    #unlist0(object$formula)                                  ,
    unlist0(object$`Residual variances homogeneity`)  ,
    unlist0(object$`Batch included`)                  ,
    # Sexual dymorphism
    unlist0(object$`Genotype contribution`$`Sexual dimorphism detected`$Criteria)   ,
    unlist0(object$`Genotype contribution`$`Sexual dimorphism detected`$Note)   ,
    # effect size
    unlist0(object$`Genotype estimate`$Value)               ,
    unlist0(object$`Sex FvKO estimate`$Value)               ,
    unlist0(object$`Sex MvKO estimate`$Value)               ,
    # Genotype
    unlist0(object$`Genotype p-value`)                  ,
    unlist0(object$`Genotype standard error`)           ,
    # sexDim
    unlist0(object$`Sex FvKO p-value`)                  ,
    unlist0(object$`Sex MvKO p-value`)                  ,
    unlist0(object$`Sex FvKO standard error`)           ,
    unlist0(object$`Sex MvKO standard error`)           ,
    # sex
    unlist0(object$`Sex p-value`)                       ,
    unlist0(object$`Sex estimate`$Value)                ,
    unlist0(object$`Sex standard error`)                ,
    # Bodyweight
    unlist0(object$`Weight p-value`)                    ,
    unlist0(object$`Weight estimate`$Value)             ,
    unlist0(object$`Weight standard error`)
  )
  return(r)
}

SelectAnalysisFE = function(object) {
  r = c(
    ######
    unlist0(object$`Applied method`)                ,
    unlist0(object$`Classification tag`$`Classification tag`) ,
    #unlist0(object$formula)                                  ,
    unlist0(object$`Residual variances homogeneity`)  ,
    unlist0(object$`Batch included`)                  ,
    # Sexual dymorphism
    unlist0(object$`Genotype contribution`$`Sexual dimorphism detected`$Criteria)   ,
    unlist0(object$`Genotype contribution`$`Sexual dimorphism detected`$Note)   ,
    # effect size
    unlist0(as.list(object$`Genotype estimate`$`Complete table`)$Value)               ,
    unlist0(as.list(object$`Sex FvKO estimate`$`Complete table`)$Value)               ,
    unlist0(as.list(object$`Sex MvKO estimate`$`Complete table`)$Value)               ,
    # Genotype
    unlist0(object$`Genotype p-value`$`Complete table`)                      ,
    unlist0(object$`Genotype standard error`)           ,
    # sexDim
    unlist0(object$`Sex FvKO p-value`$`Complete table`)                  ,
    unlist0(object$`Sex MvKO p-value`$`Complete table`)                  ,
    unlist0(object$`Sex FvKO standard error`)           ,
    unlist0(object$`Sex MvKO standard error`)           ,
    # sex
    unlist0(object$`Sex p-value`$`Complete table`)                       ,
    unlist0(as.list(object$`Sex estimate`$`Complete table`)$Value)                ,
    unlist0(object$`Sex standard error`)                ,
    # Bodyweight
    unlist0(object$`Weight p-value`)                    ,
    unlist0(object$`Weight estimate`$Value)             ,
    unlist0(object$`Weight standard error`)
  )
  return(r)
}

tableCount = function(Gen,
                      Sex,
                      levels = c('experimental.male',
                                 'experimental.female',
                                 'control.male',
                                 'control.female')) {
  t = table(interaction(Gen, Sex, sep = '.'))
  t = as.data.frame(t)
  ll = length(levels)
  r = c(rep(NA, ll))
  for (i in 1:ll)
    r[i] =   ifelse(levels[i] %in% t[, 1], t[levels[i] == t[, 1], 2], NA)
  return(r)
}


outputNames = function(){
  c1 = c(
    "Analysis date",
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
    "colony_id"
  )

  c2 = c(
    "Applied method",
    "Classification tag",
    "Residual variances homogeneity",
    "Batch included",
    "Sexual dimorphism detected",
    "Sexual dimorphism detected note",
    "Genotype estimate",
    "Sex FvKO estimate",
    "Sex MvKO estimate",
    "Genotype p-value",
    "Genotype standard error",
    "Sex FvKO p-value",
    "Sex MvKO p-value",
    "Sex FvKO standard error",
    "Sex MvKO standard error",
    "Sex p-value",
    "Sex estimate",
    "Sex standard error",
    "Weight p-value",
    "Weight estimate",
    "Weight standard error"
  )

  c3 = c(
    "Concurrent control selection",
    "is referenc gene",
    "is ignorome gene",
    "is behaviour gene",
    "variation_in_respone_after_preprocess",
    "variation_in_respone_before_preprocess",
    "NGenotype p-value",
    "WGenotype p-value"

  )

  c4 = c('MP_both',
         'MP_male',
         'MP_female')
  c5 = c(
    c1,
    "response Type",
    "applied Method",
    "windowing l",
    "windowing k",
    "Min obs required in the window",
    "the number of DateOfExperiment in the window",
    "windowing Threshold",
    "Total number of window weights",
    "sum of weights",
    paste0('N_', c2),
    paste0('W_', c2),
    c3,
    paste0('N_', c4),
    paste0('W_', c4),
    'URL',
    'Data URL',
    'Total KO male',
    'Total KO female',
    'Total WT male',
    'Total WT female'
  )
  return(c5)
}

makeIndexColumn = function(x) {
  library(digest)
  r = apply(x,1,function(y){
    digest(paste(y,sep = '~',collapse = '~'))
  })
  return(r)
}


replaceNA = function(x, replaceby = 'No MP term assigned') {
  x[is.na(x)] = replaceby
  return(x)
}

DRSummary = function(x) {
  df = data.table::fread(file = 'AllUnidimensionals.tsv',
                         header = FALSE,
                         sep = '\t')
  names(df) = outputNames()
  print(table(df$status))
  print(nrow(df))
  print(length(na.omit(df$N_MP_both)))
  print(length(na.omit(df$N_MP_male)))
  print(length(na.omit(df$N_MP_female)))

  print(length(na.omit(df$W_MP_both)))
  print(length(na.omit(df$W_MP_male)))
  print(length(na.omit(df$W_MP_female)))


  df1 = df[!is.na(df$N_MP_both) | !is.na(df$W_MP_both),]
  table(df1$N_MP_both == df1$W_MP_both,useNA = 'always')

  df2 = df[!is.na(df$N_MP_male) | !is.na(df$W_MP_male),]
  table(df2$N_MP_male == df2$W_MP_male,useNA = 'always')

  df3 = df[!is.na(df$N_MP_female) | !is.na(df$W_MP_female),]
  table(df3$N_MP_female == df3$W_MP_female,useNA = 'always')


}

calculateDiscrepancy = function(x,y){
  r = (is.na(x) & !is.na(y)) | (!is.na(x) & is.na(y))
  return(sum(r))
}

DRSummaryAcross = function(x) {
  df12 = data.table::fread(file = '/homes/hamedhm/impc_statistical_pipeline/IMPC_DRs/flatten_observations_dr12.0_15092020/SP/jobs/ExtractPvalues/resultF/data_point_of_type_unidimensional/AllUnidimensionals.tsv',
                           header = FALSE,
                           sep = '\t')

  df11 = data.table::fread(file = '/homes/hamedhm/impc_statistical_pipeline/IMPC_DRs/flatten_observations_dr11.0_16092020/SP/jobs/ExtractPvalues/resultF/data_point_of_type_unidimensional/AllDR11Unidimensionals.tsv',
                           header = FALSE,
                           sep = '\t')


  names(df11)=outputNames()
  names(df12)=outputNames()



  df11$id = makeIndexColumn(df11[,2:19])
  df12$id = makeIndexColumn(df12[,2:19])

  mall = merge(df11,df12,by = 'id',all = TRUE,suffixes = c('_dr11','_dr12'))
  m = merge(df11,df12,by = 'id',all = FALSE,suffixes = c('_dr11','_dr12'))

  table(m$N_MP_both_dr11==m$N_MP_both_dr12,useNA = 'always')
  calculateDiscrepancy(m$N_MP_both_dr11,m$N_MP_both_dr12)


  table(m$N_MP_male_dr11==m$N_MP_male_dr12,useNA = 'always')
  calculateDiscrepancy(m$N_MP_male_dr11,m$N_MP_male_dr12)

  table(m$N_MP_female_dr11==m$N_MP_female_dr12,useNA = 'always')
  calculateDiscrepancy(m$N_MP_female_dr11,m$N_MP_female_dr12)

  # window

  table(m$W_MP_both_dr11==m$W_MP_both_dr12,useNA = 'always')
  calculateDiscrepancy(m$W_MP_both_dr11,m$W_MP_both_dr12)


  table(m$W_MP_male_dr11==m$W_MP_male_dr12,useNA = 'always')
  calculateDiscrepancy(m$W_MP_male_dr11,m$W_MP_male_dr12)

  table(m$W_MP_female_dr11==m$W_MP_female_dr12,useNA = 'always')
  calculateDiscrepancy(m$W_MP_female_dr11,m$W_MP_female_dr12)


}



########## Main function
f = function(start, end, file = 'Index_DR101_V1.txt') {
  if (is.na(end))
    end = start
  ofname = paste0('R', '_', start, '-', end, '_pval.tsv')
  if (file.exists(paste0('./resultF/', ofname)))
    unlink(paste0('./resultF/', ofname))
  library(RJSONIO)
  library(DRrequiredAgeing)
  library(data.table)
  ####
  df = readLines(file)
  for (i in start:end) {
    cat(i, '|')
    file = df[i]
    if (!(
      grepl(
        pattern = 'Successful',
        x = file,
        fixed = TRUE
      ) ||
      grepl(
        pattern = 'NotProcessed',
        x = file,
        fixed = TRUE
      ) ||
      grepl(
        pattern = 'Failed',
        x = file,
        fixed = TRUE
      )
    ))
    next
    if (!file.exists(file)) {
      if (file.exists(dirname(file),    'output_Successful.tsv'))
        file = file.path(dirname(file), 'output_Successful.tsv')
      else if (file.exists(dirname(file), 'output_NotProcessed.tsv'))
        file = file.path(dirname(file),   'output_NotProcessed.tsv')
      else
        write(
          file,
          file =  paste0('DoesNotExistFiles_', Sys.Date(), '.txt'),
          append = TRUE
        )
      next
    }

    r0 = fread(
      file = file,
      header = FALSE,
      sep = '\t',
      quote = "",
      stringsAsFactors = FALSE
    )
    rN = DRrequiredAgeing:::annotationChooser(statpacket = r0,
                                              level = .0001)
    rW = DRrequiredAgeing:::annotationChooser(
      statpacket = r0,
      level = .0001,
      resultKey = 'Windowed result',
      TermKey = 'WMPTERM'
    )



    if (ncol(r0) != 20) {
      write(file,
            file = paste0('Error_Overal_', Sys.Date(), '.txt'),
            append = TRUE)
      next
    }
    r1 =
      tryCatch(
        expr =
          fromJSON(r0$V20, flatten = TRUE),
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
    #method = r1$result$detail$applied_method
    if (is.null(r1))
      next


    method = r1$Result$Details$`Applied method`
    message(paste(i, '|',  end , ':', file))
    if (is.null(method) || method  %in% c('MM', 'RR')) {
      x = c(
        unlist(r0[1, 1:19]),
        unlist0(r1$Result$Details$`Response type`) ,
        unlist0(r1$Result$Details$`Applied method`),
        ################# Window parameters
        # Have you changed that for new structure of l and k output???
        SelectWindowingParameters(object = r1$Result$Details)                 ,
        ################# VectorOutput Results
        SelectAnalysis(r1$Result$`Vector output`$`Normal result`)             ,
        SelectAnalysis(r1$Result$`Vector output`$`Windowed result`)           ,
        #SelectAnalysis(r1$Result$`Vector output`$`Full model result`)         ,
        #SelectAnalysis(r1$Result$`Vector output`$`Full model windowed result`),
        ################# Other results
        SelectOthers(r1$Result$Details)                                      ,
        ################ Ignorome/Reference/Behaviour
        SelectMISC(r0)                                                       ,
        ##### Variation in response
        unlist0(
          r1$Result$Details$'Variation in respone after preprocessing'[1]
        )   ,
        unlist0(
          r1$Result$Details$'Variation in respone before preprocessing'[1]
        )   ,
        ##### Pvals
        unlist0(
          r1$Result$`Vector output`$`Normal result`$`Genotype p-value`
        )      ,
        unlist0(
          r1$Result$`Vector output`$`Windowed result`$`Genotype p-value`
        )    ,
        #unlist0(r1$Result$`Vector output`$`Full model result`$`Genotype p-value`)  ,
        #unlist0(r1$Result$`Vector output`$`Full model windowed result`$`Genotype p-value`),
        ##### MP TERM
        DRrequiredAgeing:::StratifiedMPTerms(rN),
        DRrequiredAgeing:::StratifiedMPTerms(rW),
        ##### URL
        MakeURL(r0, r1),
        tableCount(
          Gen = r1$Result$Details$Original_biological_sample_group,
          Sex = r1$Result$Details$Original_sex
        )
      )
    } else{
      # FE only
      x = c(
        unlist(r0[1, 1:19]),
        unlist0(r1$Result$Details$`Response type`) ,
        unlist0(r1$Result$Details$`Applied method`),
        ################# Window parameters
        # Have you changed that for new structure of l and k output???
        SelectWindowingParameters(object = r1$Result$Details)                 ,
        ################# VectorOutput Results
        SelectAnalysisFE(r1$Result$`Vector output`$`Normal result`)             ,
        SelectAnalysisFE(r1$Result$`Vector output`$`Windowed result`)           ,
        #SelectAnalysis(r1$Result$`Vector output`$`Full model result`)         ,
        #SelectAnalysis(r1$Result$`Vector output`$`Full model windowed result`),
        ################# Other results
        SelectOthers(r1$Result$Details)                                      ,
        ################ Ignorome/Reference/Behaviour
        SelectMISC(r0)                                                       ,
        ##### Variation in response
        unlist0(
          r1$Result$Details$'Variation in respone after preprocessing'[1]
        )   ,
        unlist0(
          r1$Result$Details$'Variation in respone before preprocessing'[1]
        )  ,
        ##### Pvals
        unlist0(
          r1$Result$`Vector output`$`Normal result`$`Genotype p-value`$`Complete table`
        )      ,
        unlist0(
          r1$Result$`Vector output`$`Windowed result`$`Genotype p-value`$`Complete table`
        )    ,
        #unlist0(r1$Result$`Vector output`$`Full model result`$`Genotype p-value`)  ,
        #unlist0(r1$Result$`Vector output`$`Full model windowed result`$`Genotype p-value`),
        ##### MP TERM
        DRrequiredAgeing:::StratifiedMPTerms(rN),
        DRrequiredAgeing:::StratifiedMPTerms(rW),
        ##### URL
        MakeURL(r0, r1),
        tableCount(
          Gen = r1$Result$Details$Original_biological_sample_group,
          Sex = r1$Result$Details$Original_sex
        )
      )
    }
    write(
      x = paste(x, collapse = '\t'),
      file = DRrequiredAgeing:::file.path0(
        paste0(
          './resultF/',
          r1$Result$Details$`Response type`,
          '/',
          ofname
        ),
        create = TRUE,
        check  = FALSE,
        IncludedFileName = TRUE
      ),
      append = TRUE,
      ncolumns = 10 ^ 4
    )
    #}
  }
}

factor2number = function(x) {
  return(suppressWarnings(as.numeric(as.character(x))))
}

qvalueEstimator = function(x) {
  x$`NGenotype q-value` = p.adjust(factor2number(x$`NGenotype p-value`))
  x$`N_Sex FvKO q-value` = p.adjust(factor2number(x$`N_Sex FvKO p-value`))
  x$`N_Sex MvKO q-value` = p.adjust(factor2number(x$`N_Sex MvKO p-value`))

  x$`WGenotype q-value`  = p.adjust(factor2number(x$`WGenotype p-value`))
  x$`W_Sex FvKO q-value` = p.adjust(factor2number(x$`W_Sex FvKO p-value`))
  x$`W_Sex MvKO q-value` = p.adjust(factor2number(x$`W_Sex MvKO p-value`))

  return(x)
}



qvaluesGenerator = function(df,filterdfparameter=NULL) {
  df = as.data.frame(df)
  names(df) = outputNames()
  # df = df[df$status == 'Successful' &
  #           df$`applied Method` %in% c('MM', 'FE'), ]
  if (nrow(df) < 1)
    return(NULL)

  if (!is.null(filterdfparameter))
    df = droplevels(subset(df, df$parameter_stable_id == filterdfparameter))

  if (nrow(df) < 1)
    return(NULL)

  counter = 1
  d = NULL
  for (centre in unique(df$phenotyping_center)) {
    df1 = subset(df, df$phenotyping_center == centre)
    for (procedure in unique(df1$procedure_stable_id)) {
      df2 = subset(df1, df1$procedure_stable_id == procedure)
      for (parameter in unique(df2$parameter_stable_id)) {
        df3 = subset(df2, df2$parameter_stable_id == parameter)
        for (zygosity in unique(df3$zygosity)) {
          df4 = subset(df3, df3$zygosity == zygosity)
          for (strain in unique(df4$strain_accession_id)) {
            df5 = subset(df4, df4$strain_accession_id == strain)
            for (metadata in unique(df5$metadata_group)) {
              df6 = subset(df5, df5$metadata_group == metadata)
              df6 = droplevels(df6)

              df6$`NGenotype p-value`[df6$status == 'NotProcessed' & df6$'Variation in respone after preprocessing' == FALSE ]=1
              df6$`N_Sex FvKO p-value`[df6$status == 'NotProcessed' & df6$'Variation in respone after preprocessing' == FALSE ]=1
              df6$`N_Sex MvKO p-value`[df6$status == 'NotProcessed' & df6$'Variation in respone after preprocessing' == FALSE ]=1
              df6$`WGenotype p-value`[df6$status == 'NotProcessed' & df6$'Variation in respone after preprocessing' == FALSE ]=1
              df6$`W_Sex FvKO p-value`[df6$status == 'NotProcessed' & df6$'Variation in respone after preprocessing' == FALSE ]=1
              df6$`W_Sex MvKO p-value`[df6$status == 'NotProcessed' & df6$'Variation in respone after preprocessing' == FALSE ]=1


              if (counter == 1) {
                d =  qvalueEstimator(df6)
              } else{
                d = rbind(d, qvalueEstimator(df6))
              }
              counter = counter + 1
            }

          }
        }
        cat('\r-->', counter)
      }
    }
  }
  if(!is.null(d) && nrow(d)>0){
    #d = d[, colSums(is.na(d)) < nrow(d), drop = FALSE]
    d = d[!duplicated(d),]
  }
  return(d)

}



qvalue2AllZips = function(path = getwd()) {
  files = list.files(
    path = path,
    pattern = '.zip',
    all.files = TRUE,
    full.names = TRUE,
    include.dirs = FALSE
  )
  for (file in files) {
    df = NULL
    df = data.table::fread(
      file = DRrequiredAgeing::UnzipAndfilePath(file = file)[1],
      header = FALSE,
      sep = '\t',
      stringsAsFactors = TRUE
    )
    df = qvaluesGenerator(df)
    write.csv(df,
              file = paste0(basename(file), '.csv'),
              row.names = FALSE)
  }
}

parameter2qvalue = function(parameter, file) {
  df  = data.table::fread(
    file = DRrequiredAgeing::UnzipAndfilePath(file = file)[1],
    header = FALSE,
    sep = '\t',
    stringsAsFactors = TRUE
  )

  dir = 'QvalueResults' # DRrequiredAgeing:::RemoveSpecialChars(basename(file))

  if(!dir.exists(dir))
    dir.create(dir,recursive = TRUE)

  d = qvaluesGenerator(df = as.data.frame(df), filterdfparameter = parameter)
  if(!is.null(d) || nrow(d)>0){
    write.csv(
      d,
      file = paste0(
        dir,
        '/',
        basename(file),
        '_',
        DRrequiredAgeing:::randomIdGenerator(l = 16),
        '.csv'
      ),
      row.names = FALSE
    )
  }else{
    return(NULL)
  }

}

makejobs = function(path = getwd()) {
  files = list.files(
    path = path,
    pattern = '.tsv',
    all.files = TRUE,
    full.names = TRUE,
    include.dirs = FALSE
  )
  for (file in files) {
    df  = data.table::fread(
      file = DRrequiredAgeing::UnzipAndfilePath(file = file)[1],
      header = FALSE,
      sep = '\t',
      stringsAsFactors = TRUE
    )
    names(df) = outputNames()

    parameters = unique(df$parameter_stable_id)

    if (!dir.exists('err'))
      dir.create('err')

    if (!dir.exists('out'))
      dir.create('out')

    bf = basename(file)
    n = length(parameters)
    jobs = paste0 (
      'bsub -M 8000 -e err/err',
      bf,
      1:n,
      ' -o out/out',
      bf,
      1:n,
      ' "Rscript ExtractPVals.R ',
      parameters,
      ' ',
      file,
      '"'
    )
    write(
      x = paste(jobs,sep = '\n',collapse = '\n'),
      file = 'Jobs.bch',
      ncolumns = 10000
    )
  }

}

# THE LINES BELOW MUST NOT CHANGE AT ALL (used in the statspipeline 1line command)
#makejobs()
#qvalue2AllZips()
#parameter2qvalue(args[1], args[2])
ignore.my.name = f(start =  as.numeric(args[1]), end = as.numeric(args[2]),file = args[3])





# After getting the table
# library(plyr)
# TotalN = length(OnlyChanges_P_Value_results_2018_08_08_10_59_03$Procedure_group)
# cdata <-
# 	ddply(
# 		OnlyChanges_P_Value_results_2018_08_08_10_59_03,
# 		c("Centre","Procedure_group"),
# 		summarise,
# 		N        = length(`Gain (windowing significant but not normal)`),
# 		CentreContribution = paste0(round(N / TotalN * 100,1), '%'),
# 		nGain    = sum(`Gain (windowing significant but not normal)`),
# 		nLoss    = sum(`Loss (Normal significant but not windowing)`),
# 		#GainLossRatio = round(nGain / nLoss, 2),
# 		Gain2TotalInCentre    = paste0(round(nGain / N * 100), '%'),
# 		Loss2TotalInCentre    = paste0(round(nLoss / N * 100), '%')
#
# 	)
# cdata[order(cdata$N),]
