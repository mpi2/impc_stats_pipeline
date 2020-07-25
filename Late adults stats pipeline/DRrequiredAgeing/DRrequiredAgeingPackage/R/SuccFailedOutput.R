# all in small case and separated by underscore
SuccessfulOutput = function(args, writeOutputToDB = FALSE, b64Encode = TRUE, plainOutput = FALSE) {
  requireNamespace("DBI")
  ### 1 Experiment detail
  experiment_details      = list(
    ####
    StatPacketCreationDate= as.character(Sys.time())                                             , #0
    status                =  ifelse(
      !NullOrError(args$c.ww0$NormalObj)      &&
        NullOrError(args$c.ww0$NormalObj$messages),
      'Successful',
      'Failed'
    )                                                                                            , #1
    procedure_group       = UniqueAndNNull(args$n3.5$procedure_group,removeSpecials = FALSE)     , #2
    procedure_stable_id   = UniqueAndNNull(args$n3.5$procedure_stable_id,removeSpecials = FALSE) , #3
    procedure_name        = UniqueAndNNull(args$n3.5$procedure_name,removeSpecials = FALSE)      , #4
    parameter_stable_id   = UniqueAndNNull(args$parameter,removeSpecials = FALSE)                , #5
    parameter_name        = UniqueAndNNull(args$n3.5$parameter_name,removeSpecials = FALSE)      , #6
    phenotyping_center    = UniqueAndNNull(args$center,removeSpecials = FALSE)                   , #7
    allele_symbol         = UniqueAndNNull(args$n3.5$allele_symbol,removeSpecials = FALSE)       , #8
    allele_name           = UniqueAndNNull(args$n3.5$allelic_composition,removeSpecials = FALSE) , #8_1
    allele_accession_id   = UniqueAndNNull(args$n3.5$allele_accession_id,removeSpecials = FALSE) , #9
    gene_symbol           = UniqueAndNNull(args$n3.5$gene_symbol,removeSpecials = FALSE)         , #10
    gene_accession_id     = UniqueAndNNull(args$n3.5$gene_accession_id,removeSpecials = FALSE)   , #11
    pipeline_name         = UniqueAndNNull(args$n3.5$pipeline_name,removeSpecials = FALSE)       , #12
    pipeline_stable_id    = UniqueAndNNull(args$n3.5$pipeline_stable_id,removeSpecials = FALSE)  , #13
    strain_accession_id   = UniqueAndNNull(args$strain, removeSpecials = FALSE)              , #14
    metadata_group        = UniqueAndNNull(args$meta, removeSpecials = FALSE )               , #15
    zygosity              = UniqueAndNNull(args$zyg, removeSpecials = FALSE  )               , #16
    colony_id             = UniqueAndNNull(args$colony, removeSpecials = FALSE )               #17
  )


  #args$note$`Raw data summary statistics_2`=dictionary2listConvert(args$note$`Raw data summary statistics`)
  #args$note$`OpenStatsList object summary statistics_2`=dictionary2listConvert(args$note$`OpenStatsList object summary statistics`)


  ######## 2 JSON
  message0('Forming the list before applying JSON transformation ...')
  args$note$'Experiment detail' = experiment_details[-1]
  listDetails                   = list('Details' = sortList(args$note))
  listVectorOutput              = list('Vector output' = args$c.ww.vec$list)
  FinalList                     = list('Result' = c(listVectorOutput, listDetails))
  JsonObj                       = FinalJson2ObjectCreator(FinalList = cleanNULLkeys(FinalList))

  ####### Just to get names - useful in post processing data (like in pvalue extraction)
  if (plainOutput)
    return(experiment_details)
  ######## 2 CSV
  outP =   c(unlist(experiment_details),
             base64(x =
                      JsonObj,
                    active = b64Encode && args$encode))

  # Write to the SQLite DB
  if (writeOutputToDB)
  WriteToDB(outP,
            dbname = paste0(experiment_details$gene_accession_id[1], "_DR10SQLite.db"))

  return(outP)
}

