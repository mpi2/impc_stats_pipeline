# all in small case and separated by underscore
SuccessfulOutput = function(args) {
  ######## 1 JSON
  listDetails             = list(details = args$note)
  listVectorOutput        = list(vectoroutput = args$c.ww.vec$list)
  FinalList               = list(result = c(listVectorOutput, listDetails))
  JsonObj                 = FinalJsonBobectCreator(FinalList = FinalList)
  # Sort the list based on the size of its elements
  #JsonObj$details = JsonObj$details[order(sapply(JsonObj$details,length),decreasing=T)]

  ######## 2 CSV
  outP =   c(
    ifelse(
      NullOrError(args$c.ww0$NormalObj$value),
      'Failed',
      'Successful'
    ),
    args$procedure                                                         ,
    UniqueAndNNull(args$n3.5$procedure_stable_id,removeSpecials = FALSE)   ,
    UniqueAndNNull(args$n3.5$procedure_name,removeSpecials = FALSE)        ,
    args$parameter                                                         ,
    UniqueAndNNull(args$n3.5$parameter_name,removeSpecials = FALSE)        ,
    args$center                                                            ,
    UniqueAndNNull(args$n3.5$allele_symbol,removeSpecials = FALSE)         ,
    UniqueAndNNull(args$n3.5$gene_symbol,removeSpecials = FALSE)           ,
    UniqueAndNNull(args$n3.5$gene_accession_id,removeSpecials = FALSE)     ,
    UniqueAndNNull(args$n3.5$pipeline_name,removeSpecials = FALSE)         ,
    UniqueAndNNull(args$n3.5$pipeline_stable_id,removeSpecials = FALSE)    ,
    args$strain                                     ,
    args$meta                                       ,
    args$zyg                                        ,
    args$colony                                     ,
    'NA'                                            ,
    base64(x =
             JsonObj,
           active = args$encode)
  )
  return(outP)
}