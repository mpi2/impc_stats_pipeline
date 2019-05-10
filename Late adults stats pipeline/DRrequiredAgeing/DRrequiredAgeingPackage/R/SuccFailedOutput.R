# all in small case and separated by underscore
SuccessfulOutput = function(args) {

  ### 1 Experiment detail
  experiment_details      = list(
    ####
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
    #allele_symbol         = UniqueAndNNull(args$n3.5$allele_symbol,removeSpecials = FALSE)      , #8
    allele_accession_id   = UniqueAndNNull(args$n3.5$allele_accession_id,removeSpecials = FALSE) , #8
    gene_symbol           = UniqueAndNNull(args$n3.5$gene_symbol,removeSpecials = FALSE)         , #9
    gene_accession_id     = UniqueAndNNull(args$n3.5$gene_accession_id,removeSpecials = FALSE)   , #10
    pipeline_name         = UniqueAndNNull(args$n3.5$pipeline_name,removeSpecials = FALSE)       , #11
    pipeline_stable_id    = UniqueAndNNull(args$n3.5$pipeline_stable_id,removeSpecials = FALSE)  , #12
    strain_accession_id   = args$strain               , #13
    metadata_group        = args$meta                 , #14
    zygosity              = args$zyg                  , #15
    colony_id             = args$colony               , #16
    reserved              = 'NA'                        #17
  )

  ######## 2 JSON
  message0('Forming the list before applying JSON transformation ...')
  args$note$experiment_details = experiment_details
  listDetails                  = list(details = sortList(args$note))
  listVectorOutput             = list(vectoroutput = args$c.ww.vec$list)
  FinalList                    = list(result = c(listVectorOutput, listDetails))
  JsonObj                      = FinalJsonBobectCreator(FinalList = FinalList)

  ######## 2 CSV
  outP =   c(
    unlist(experiment_details),
    base64(x =
             JsonObj,
           active = args$encode)
  )
  return(outP)
}