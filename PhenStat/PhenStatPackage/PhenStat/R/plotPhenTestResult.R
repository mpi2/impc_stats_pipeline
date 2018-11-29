plot.PhenTestResult = function (x                    ,
                                graphingName = NULL  ,
                                outputMessages = TRUE,
                                type = NULL,...    ) {

  phenTestResult = x
  if (is.null(phenTestResult@depVariable))
    stop('Wrong dataset is imported!')


  funList <-
    unclass(lsf.str(envir = asNamespace("PhenStat"), all = TRUE))
  if (!is.null(type) && !(type %in% c(funList)))
    stop ('"type" is not recognised, see the manual')

  atLeastOne = FALSE
  if ((is.null(type) ||
       'boxplotSexGenotypeResult'  %in% type) &&
      phenTestResult@method == 'MM') {
    boxplotSexGenotypeResult (phenTestResult,
                              graphingName = graphingName,
                              outputMessages = outputMessages)
    atLeastOne = TRUE
  }

  if ((is.null(type) ||
       'scatterplotSexGenotypeBatchResult'  %in% type) &&
      phenTestResult@method == 'MM') {
    scatterplotSexGenotypeBatchResult (phenTestResult,
                                       graphingName = graphingName,
                                       outputMessages = outputMessages)
    atLeastOne = TRUE
  }


  if ((is.null(type) ||
       'scatterplotGenotypeWeightResult'  %in% type) &&
      phenTestResult@method == 'MM') {
    if('Weight' %in%  names(phenTestResult@analysedDataset)){
      scatterplotGenotypeWeightResult (phenTestResult,
                                       graphingName = graphingName,
                                       outputMessages = outputMessages)
    }else{
      message('Column Weight does not found.')
    }
    atLeastOne = TRUE
  }

  if ((is.null(type) ||
       'qqplotGenotype'  %in% type) &&
      phenTestResult@method == 'MM') {
    qqplotGenotype (phenTestResult, outputMessages = outputMessages)
    atLeastOne = TRUE
  }

  if ((is.null(type) ||
       'plotResidualPredicted'  %in% type) &&
      phenTestResult@method == 'MM') {
    plotResidualPredicted (phenTestResult, outputMessages = outputMessages)
    atLeastOne = TRUE
  }

  if ((is.null(type) ||
       'qqplotRandomEffects'  %in% type) &&
      phenTestResult@method == 'MM') {
    qqplotRandomEffects (phenTestResult, outputMessages = outputMessages)
    atLeastOne = TRUE
  }

  if ((is.null(type) ||
       'boxplotResidualBatch'  %in% type) &&
      phenTestResult@method == 'MM') {
    boxplotResidualBatch (phenTestResult, outputMessages = outputMessages)
    atLeastOne = TRUE
  }

  if ((is.null(type) ||
       'qqplotRotatedResiduals'  %in% type) &&
      phenTestResult@method == 'MM') {
    qqplotRotatedResiduals (phenTestResult, outputMessages = outputMessages)
    atLeastOne = TRUE
  }

  if ((is.null(type) ||
       'categoricalBarplot'  %in% type) &&
      phenTestResult@method == 'FE') {
    categoricalBarplot (phenTestResult, outputMessages = outputMessages)
    atLeastOne = TRUE
  }

  if (!atLeastOne)
    message('No plot found for this dataset.')
}
