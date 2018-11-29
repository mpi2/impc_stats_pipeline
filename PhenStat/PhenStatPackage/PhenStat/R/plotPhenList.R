plot.PhenList =  function(x                    ,
                          depVariable = 'Value',
                          graphingName = NULL  ,
                          outputMessages = TRUE,
                          type = NULL,...       )
{
  phenList  = x
  if (is.null(depVariable))
    stop('Please define the dependent variable!')

  ## Test: depVariable is numeric
  x                = getDataset(phenList)
  columnOfInterest = x[, c(depVariable)]
  if (!is.numeric(columnOfInterest)){
    message('This function is only defined for the numeric responses.')
  }else{
    funList <-
      unclass(lsf.str(envir = asNamespace("PhenStat"), all = TRUE))
    if (!is.null(type) && !(type %in% c(funList)))
      stop ('"type" is not recognised, see the manual')

    atLeastOne = FALSE
    if (is.null(type) || 'boxplotSexGenotype' %in% type) {
      boxplotSexGenotype (
        phenList,
        depVariable = depVariable,
        graphingName = graphingName,
        outputMessages = outputMessages
      )
      atLeastOne = TRUE
    }

    if (is.null(type) || 'boxplotSexGenotypeBatchAdjusted' %in% type) {
      boxplotSexGenotypeBatchAdjusted (
        phenList,
        depVariable = depVariable,
        graphingName = graphingName,
        outputMessages = outputMessages
      )
      atLeastOne = TRUE
    }

    if (is.null(type) ||
        'boxplotSexGenotypeWeightBatchAdjusted' %in% type) {
      if('Weight' %in%  names(phenList@datasetPL)){
        boxplotSexGenotypeWeightBatchAdjusted(
          phenList,
          depVariable = depVariable,
          graphingName = graphingName,
          outputMessages = outputMessages
        )
      }else{
        message('Column Weight does not found.')
      }
      atLeastOne = TRUE
    }

    if (is.null(type) ||  'boxplotSexGenotypeBatch' %in% type) {
      boxplotSexGenotypeBatch (
        phenList,
        depVariable = depVariable,
        graphingName = graphingName,
        outputMessages = outputMessages
      )
      atLeastOne = TRUE
    }

    if (is.null(type) || 'scatterplotSexGenotypeBatch' %in% type) {
      scatterplotSexGenotypeBatch (
        phenList,
        depVariable = depVariable,
        graphingName = graphingName,
        outputMessages = outputMessages
      )
      atLeastOne = TRUE
    }

    if (is.null(type) || 'scatterplotGenotypeWeight' %in% type) {
      scatterplotGenotypeWeight (
        phenList,
        depVariable = depVariable,
        graphingName = graphingName,
        outputMessages = outputMessages
      )
      atLeastOne = TRUE
    }
    if (!atLeastOne)
      message('No plot found for this dataset.')
  }
}
