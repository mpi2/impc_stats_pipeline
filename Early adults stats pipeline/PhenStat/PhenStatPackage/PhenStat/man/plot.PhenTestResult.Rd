\name{plot.PhenTestResult}
\alias{plot.PhenTestResult}
\title{Plot phenTestResult object}
\description{
  Plot method for objects of class "phenTestResult".
}
\usage{
\method{plot}{PhenTestResult}(
    x                    ,
    graphingName = NULL  ,
    outputMessages = TRUE,
    type = NULL,
    ...
  )
}
%-maybe also 'usage' for other objects documented here.
\arguments{
  \item{
    x
  }{
    instance of the PhenTestResult class
    mandatory argument
  }
  \item{
    graphingName
  }{
    character string defining the name to be used on the graph for the dependent variable
  }
  \item{
    outputMessages
  }{
    flag:"FALSE" value to suppress output messages
    "TRUE" value to show output messages
    default value TRUE
  }
  \item{
    type
  }{
    a vector of names. Select one or more from the list below to only get those plots.
    \itemize{
      \item boxplotSexGenotypeResult. See documentations for {
        {
          boxplotSexGenotypeResult
        }
      }
      \item scatterplotSexGenotypeBatchResult. See documentations for {
        {
          scatterplotSexGenotypeBatchResult
        }
      }
      \item scatterplotGenotypeWeightResult. See documentations for {
        {
          scatterplotGenotypeWeightResult
        }
      }
      \item plotResidualPredicted. See documentations for {
        {
          plotResidualPredicted
        }
      }
      \item qqplotRandomEffects. See documentations for {
        {
          qqplotRandomEffects
        }
      }
      \item boxplotResidualBatch. See documentations for {
        {
          boxplotResidualBatch
        }
      }
      \item qqplotRotatedResiduals. See documentations for {
        {
          qqplotRotatedResiduals
        }
      }
      \item qqplotGenotype. See documentations for {
        {
          qqplotGenotype
        }
      }
      \item categoricalBarplot. See documentations for {
        {
          categoricalBarplot
        }
      }
    }
  }
  \item{\dots}{
  }
}

\seealso{
  \code{\link{testDataset}}
}

\author{
  Hamed Haselimashhadi < hamedhm@eb.ac.uk >
}


\examples{
  example(testDataset)
}

\value{
    Empty.
}
