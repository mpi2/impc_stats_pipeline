# all in small case and separated by underscore
NotProcessedOutput = function(args) {
  bsg3.5        = args$n3.5$biological_sample_group
  bsg3.5.2      = args$n3.5.2$biological_sample_group
  cid3.5.2      = args$n3.5.2$colony_id
  n3.5.2OnlyKO  = subset(args$n3.5.2,args$n3.5.2$biological_sample_group %in% 'experimental')
  ######## 1 LIST
  NotProcessedLogics = list(
    is_exception   =  args$isException,
    empty_dataset_after_preprocess  = is.null(args$n3.5.2),
    empty_response_after_preprocess = !length(na.omit(args$n3.5.2[, args$depVar])) > 0,
    variation_in_respone_after_preprocess =   ifelse(
      args$depVariable$accepted,
      NonZeroVariation(args$n3.5.2[, args$depVar]),
      'Not numeric or factor response'
    ),
    both_mut_and_control_after_preprocess = list(
      criteria_result = length(unique(bsg3.5.2)) > 1 ,
      levels          = unique(bsg3.5.2)
    ),
    min_onbs_in_each_group_raw_data_before_preprocess = list(
      criteria_result   = min0(table(bsg3.5)) >= args$minSampRequired,
      threshold         = args$minSampRequired,
      stage             = 'before_preprocessing',
      min_obs_in_data   = min0(table(bsg3.5))
    ),
    min_onbs_in_each_group_processed_data_after_preprocess = list(
      criteria_result   = min0(table(bsg3.5.2)) >= args$minSampRequired,
      threshold         = args$minSampRequired,
      stage             = 'after_preprocessing',
      min_obs_in_data = min0(table(bsg3.5.2))
    ),
    max_mutants_in_genotype_sex_table_after_preprocess = list(
      criteria_result   = max0(table(n3.5.2OnlyKO$biological_sample_group,n3.5.2OnlyKO$sex)) > 1,
      threshold         = 1,
      stage             = 'after_preprocessing',
      max_mutants_in_genotype_sex_table = max0(table(n3.5.2OnlyKO$biological_sample_group,n3.5.2OnlyKO$sex))
    ),
    the_num_colonies_after_preprocess = list(
      #criteria_result = length(unique(args$n3.5.12$colony_id)) > 1,
      criteria_result = length(RepBlank(unique(cid3.5.2), match = c('', NA, 'NA'))) > 1,
      threshold       = 2,
      colonies        = RepBlank(unique(cid3.5.2), match = c('', NA, 'NA'))
    ),
    min_onbs_required_for_rr = list(
      criteria_result = RR_thresholdCheck(
        data = args$n3.5.2,
        depVar = args$depVar,
        parameter = args$parameter,
        methodmap = args$methodmap
      )
    )
  )

  ### 2 Experiment detail
  experiment_details      = list(
    ####
    status                = 'NotProcessed'                                                      , #1
    procedure_group       = args$procedure                                                      , #2
    procedure_stable_id   = UniqueAndNNull(args$n3.5$procedure_stable_id,removeSpecials = FALSE), #3
    procedure_name        = UniqueAndNNull(args$n3.5$procedure_name,removeSpecials = FALSE)     , #4
    parameter_stable_id   = args$parameter                                                      , #5
    parameter_name        = UniqueAndNNull(args$n3.5$parameter_name,removeSpecials = FALSE)     , #6
    phenotyping_center    = args$center                                                         , #7
    allele_symbol         = UniqueAndNNull(args$n3.5$allele_symbol,removeSpecials = FALSE)      , #8
    gene_symbol           = UniqueAndNNull(args$n3.5$gene_symbol,removeSpecials = FALSE)        , #9
    gene_accession_id     = UniqueAndNNull(args$n3.5$gene_accession_id,removeSpecials = FALSE)  , #10
    pipeline_name         = UniqueAndNNull(args$n3.5$pipeline_name,removeSpecials = FALSE)      , #11
    pipeline_stable_id    = UniqueAndNNull(args$n3.5$pipeline_stable_id,removeSpecials = FALSE) , #12
    strain_accession_id   = args$strain               , #13
    metadata_group        = args$meta                 , #14
    zygosity              = args$zyg                  , #15
    colony_id             = args$colony               , #16
    reserved              = 'NA'                        #17
  )

  ### 3 JSON
  message0('Forming the list before applying JSON transformation ...')
  args$note$experiment_details = experiment_details
  listDetails                  = list(details = sortList(c(
    NotProcessedLogics,
    args$note
  )))
  listVectorOutput        = list(vectoroutput = NULL)
  FinalList               = list(result = c(listVectorOutput, listDetails))
  JsonObj                 = FinalJsonBobectCreator(FinalList = FinalList)
  ######## 3 CSV
  optFail =   c(
    unlist(experiment_details),
    base64(x =
             JsonObj,
           active = args$encode)
  )

  return(optFail)
}
