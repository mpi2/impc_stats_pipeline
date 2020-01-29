# all in small case and separated by underscore
SuccessfulOutput = function(args) {
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
    procedure_group       = args$procedure                                                       , #2
    procedure_stable_id   = UniqueAndNNull(args$n3.5$procedure_stable_id,removeSpecials = FALSE) , #3
    procedure_name        = UniqueAndNNull(args$n3.5$procedure_name,removeSpecials = FALSE)      , #4
    parameter_stable_id   = args$parameter                                                       , #5
    parameter_name        = UniqueAndNNull(args$n3.5$parameter_name,removeSpecials = FALSE)      , #6
    phenotyping_center    = args$center                                                          , #7
    allele_symbol         = UniqueAndNNull(args$n3.5$allele_symbol,removeSpecials = FALSE)       , #8
    allele_accession_id   = UniqueAndNNull(args$n3.5$allele_accession_id,removeSpecials = FALSE) , #9
    gene_symbol           = UniqueAndNNull(args$n3.5$gene_symbol,removeSpecials = FALSE)         , #10
    gene_accession_id     = UniqueAndNNull(args$n3.5$gene_accession_id,removeSpecials = FALSE)   , #11
    pipeline_name         = UniqueAndNNull(args$n3.5$pipeline_name,removeSpecials = FALSE)       , #12
    pipeline_stable_id    = UniqueAndNNull(args$n3.5$pipeline_stable_id,removeSpecials = FALSE)  , #13
    strain_accession_id   = args$strain               , #14
    metadata_group        = args$meta                 , #15
    zygosity              = args$zyg                  , #16
    colony_id             = args$colony                 #17
  )

  ######## 2 JSON
  message0('Forming the list before applying JSON transformation ...')
  args$note$'Experiment detail' = experiment_details[-1]
  listDetails                   = list('Details' = sortList(args$note))
  listVectorOutput              = list('Vector output' = args$c.ww.vec$list)
  FinalList                     = list('Result' = c(listVectorOutput, listDetails))
  JsonObj                       = FinalJson2ObjectCreator(FinalList = cleanNULLkeys(FinalList))

  ######## 2 CSV
  outP =   c(
    unlist(experiment_details),
    base64(x =
             JsonObj,
           active = args$encode)
  )

  # Write to the SQLite DB
  WriteToDB(outP,
            dbname = paste0(experiment_details$gene_accession_id[1], "_DR10SQLite.db"))

  return(outP)
}

