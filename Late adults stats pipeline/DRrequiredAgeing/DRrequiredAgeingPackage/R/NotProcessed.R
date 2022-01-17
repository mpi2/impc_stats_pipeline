# all in small case and separated by underscore
NotProcessedOutput = function(args, writeOutputToDB = FALSE) {
  requireNamespace("DBI")
  bsg3.5        = droplevels0(args$n3.5$biological_sample_group)
  bsg3.5.2      = droplevels0(args$n3.5.2$biological_sample_group)
  cid3.5.2      = NA2LabelInFactor(args$n3.5.2$colony_id)
  cid3.5        = NA2LabelInFactor(args$n3.5$colony_id)
  n3.5.2OnlyKO  = droplevels0(subset(args$n3.5.2,args$n3.5.2$biological_sample_group %in% 'experimental'))
  ######## 1 LIST
  NotProcessedLogics = list(
    'Is exception'   =  args$isException,
    'Empty dataset after preprocessing'  = is.null(args$n3.5.2),
    'Empty response after preprocessing' = !length(na.omit(args$n3.5.2[, args$depVar])) > 0,
    'Variation in respone after preprocessing' =   list(
      'Criteria result' = ifelse(
        args$depVariable$accepted,
        NonZeroVariation(args$n3.5.2[, args$depVar]),
        'Not numeric or factor response'
      ),
      'Threshold' = 0,
      'Model'     = 'R base unique() function'
    ),
    'Variation in respone before preprocessing' =   list(
      'Criteria result' = columnLevelsVariationRadio(dataset = args$n3.5.2, columnName = args$depVar) > 0.005,
      'Value'           = columnLevelsVariationRadio(dataset = args$n3.5.2, columnName = args$depVar)        ,
      'Threshold'       = 0.005,
      'Model'           = 'PhenStat variation function'
    ),
    'Both mutants and controls in the data after preprocessing' = list(
      'Criteria result' = length(unique(bsg3.5.2)) > 1 ,
      'Levels'          = if (length(unique(bsg3.5.2)) > 0) {
        unique(bsg3.5.2)
      } else{
        'No colony found'
      }
    ),
    'Min number of observations required in each Genotype level in the raw data before preprocessing' = list(
      'Criteria result'   = min0(table(bsg3.5)) >= args$minSampRequired,
      'Threshold'         = args$minSampRequired,
      'Stage'             = 'Before preprocessing',
      'Min number of observations in data'   = min0(table(bsg3.5))
    ),
    'Min number of observations required in each Genotype level in the raw data after preprocessing' = list(
      'Criteria result'   = min0(table(bsg3.5.2)) >= args$minSampRequired,
      'Threshold'         = args$minSampRequired,
      'Stage'             = 'After preprocessing',
      'Min observations in data' = min0(table(bsg3.5.2))
    ),
    'Max number of mutants in GenotypeSex interaction table after preprocessing' = list(
      'Criteria result'   = max0(table(n3.5.2OnlyKO$biological_sample_group,n3.5.2OnlyKO$sex)) > 1,
      'Threshold'         = 1,
      'Stage'             = 'After preprocessing',
      'Max number of mutants in GenotypeSex table' = max0(table(n3.5.2OnlyKO$biological_sample_group,n3.5.2OnlyKO$sex))
    ),
    'Total number of colonies before preprocessing' = list(
      'Criteria result' = length(RepBlank(unique(cid3.5), match = c('', NA, 'NA'))) > 1,
      'Threshold'       = 2,
      'Colonies'        = RepBlank(unique(cid3.5), match = c('', NA, 'NA'))
    ),
	'Total number of colonies after preprocessing' = list(
      'Criteria result' = length(RepBlank(unique(cid3.5.2), match = c('', NA, 'NA'))) > 1,
      'Threshold'       = 2,
      'Colonies'        = if (length(RepBlank(unique(cid3.5.2), match = c('', NA, 'NA'))) > 0) {
        RepBlank(unique(cid3.5.2), match = c('', NA, 'NA'))
      } else{
        'No colony found'
      }
    ),
    'Min observations required for RR method' =
       RR_thresholdCheck(
        data = args$n3.5.2,
        depVar = args$depVar,
        parameter = args$parameter,
        methodmap = args$methodmap
      )
  )

  ### 2 Experiment detail
  experiment_details      = list(
    ####
    StatPacketCreationDate= as.character(Sys.time())                                            , #0
    status                = 'NotProcessed'                                                      , #1
    procedure_group       = UniqueAndNNull(args$n3.5$procedure_group,removeSpecials = FALSE)    , #2
    procedure_stable_id   = UniqueAndNNull(args$n3.5$procedure_stable_id,removeSpecials = FALSE), #3
    procedure_name        = UniqueAndNNull(args$n3.5$procedure_name,removeSpecials = FALSE)     , #4
    parameter_stable_id   = args$parameter                                                      , #5
    parameter_name        = UniqueAndNNull(args$n3.5$parameter_name,removeSpecials = FALSE)     , #6
    phenotyping_center    = args$center                                                         , #7
    allele_symbol         = UniqueAndNNull(args$n3.5$allele_symbol,removeSpecials = FALSE)      , #8
    allele_name           = UniqueAndNNull(args$n3.5$allelic_composition,removeSpecials = FALSE) , #8_1
    allele_accession_id   = UniqueAndNNull(args$n3.5$allele_accession_id,removeSpecials = FALSE), #9
    gene_symbol           = UniqueAndNNull(args$n3.5$gene_symbol,removeSpecials = FALSE)        , #10
    gene_accession_id     = UniqueAndNNull(args$n3.5$gene_accession_id,removeSpecials = FALSE)  , #11
    pipeline_name         = UniqueAndNNull(args$n3.5$pipeline_name,removeSpecials = FALSE)      , #12
    pipeline_stable_id    = UniqueAndNNull(args$n3.5$pipeline_stable_id,removeSpecials = FALSE) , #13
    strain_accession_id   = args$strain               , #14
    metadata_group        = args$meta                 , #15
    zygosity              = args$zyg                  , #16
    colony_id             = args$colony                 #17
  )

  ### 3 JSON
  message0('Forming the list before applying JSON transformation ...')
  args$note$'Experiment detail'  = experiment_details[-1]

  # Add statpacket ids
  datasignFormula = '~colony_id+data_point+category+discrete_point+phenotyping_center+procedure_group+procedure_stable_id+parameter_stable_id+allele_accession_id+gene_accession_id+pipeline_stable_id+zygosity+metadata_group+sex+biological_sample_group+strain_accession_id+weight+date_of_experiment+date_of_birth+time_point+text_value'
  args$note$statpacket_raw_id       = OpenStats:::dataSignature(formula = paste0(datasignFormula,'+DRversion'),data = args$n3.5,digits = 10)
  args$note$statpacket_universal_id = digest::digest(OpenStats:::dataSignature(formula = paste0(datasignFormula,'+DRversion'),data = args$n3.5,digits = 10))
  args$note$statpacket_stable_id    = digest::digest(OpenStats:::dataSignature(formula = datasignFormula,data = args$n3.5,digits = 10))


  listDetails                    = list('Details' = sortList(c(
    NotProcessedLogics,
    args$note
  )))
  listVectorOutput        = list('Vector output' = NULL)
  FinalList               = list('Result' = c(listVectorOutput, cleanNULLkeys(listDetails)))
  JsonObj                 = FinalJson2ObjectCreator(FinalList = FinalList)
  ######## 3 CSV
  optFail =   c(
    unlist(experiment_details),
    base64(x =
             JsonObj,
           active = args$encode)
  )
  # Write to the SQLite DB
  if (writeOutputToDB)
  WriteToDB(optFail,
            dbname = paste0(experiment_details$gene_accession_id[1], "_DR10SQLite.db"))

    return(optFail)
}



