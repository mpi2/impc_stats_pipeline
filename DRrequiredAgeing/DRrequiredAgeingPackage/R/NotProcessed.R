# all in small case and separated by underscore
NotProcessedOutput = function(args) {
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
      criteria_result = length(unique(args$n3.5.2$biological_sample_group)) > 1 ,
      levels          = unique(args$n3.5.2$biological_sample_group)
    ),
    min_onbs_in_each_group_raw_data_before_preprocess = list(
      criteria_result   = min0(table(args$n3.5$biological_sample_group)) >= args$minSampRequired,
      threshold         = args$minSampRequired,
      stage             = 'before_preprocessing',
      min_obs_in_data   = min0(table(args$n3.5$biological_sample_group))
    ),
    min_onbs_in_each_group_processed_data_after_preprocess = list(
      criteria_result   = min0(table(args$n3.5.2$biological_sample_group)) >= args$minSampRequired,
      threshold         = args$minSampRequired,
      stage             = 'after_preprocessing',
      min_obs_in_data = min0(table(args$n3.5.2$biological_sample_group))
    ),
    the_num_colonies_after_preprocess = list(
      #criteria_result = length(unique(args$n3.5.12$colony_id)) > 1,
      criteria_result = length(RepBlank(unique(args$n3.5.2$colony_id), match = c('', NA, 'NA'))) > 1,
      threshold       = 2,
      colonies        = RepBlank(unique(args$n3.5.2$colony_id), match = c('', NA, 'NA'))
    )
  )
  ### 2 JSON
  listDetails             = list(details = c(NotProcessedLogics, args$note))
  listVectorOutput        = list(vectoroutput = NULL)
  FinalList               = list(result = c(listVectorOutput, listDetails))
  JsonObj                 = FinalJsonBobectCreator(FinalList = FinalList)

  ######## 3 CSV
  optFail =   c(
    ####
    'NotProcessed'                                  , #1
    args$procedure                                  , #2
    UniqueAndNNull(args$n3.5$procedure_stable_id)   , #3
    UniqueAndNNull(args$n3.5$procedure_name)        , #4
    args$parameter                                  , #5
    UniqueAndNNull(args$n3.5$parameter_name)        , #6
    args$center                                     , #7
    UniqueAndNNull(args$n3.5$allele_symbol)         , #8
    UniqueAndNNull(args$n3.5$gene_symbol)           , #9
    UniqueAndNNull(args$n3.5$gene_accession_id)     , #10
    UniqueAndNNull(args$n3.5$pipeline_name)         , #11
    UniqueAndNNull(args$n3.5$pipeline_stable_id)    , #12
    args$strain               , #13
    args$meta                 , #14
    args$zyg                  , #15
    args$colony               , #16
    'NA'                      , #17
    base64(x =
             JsonObj,
           active = args$encode)
  )

  return(optFail)
}
