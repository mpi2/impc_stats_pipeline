\name{plot.PhenList}
\alias{plot.PhenList}
\title{Plot Phenlist object}
\description{
  Plot method for objects of class "PhenList".
}
\usage{
\method{plot}{PhenList}(
    x                    ,
    depVariable = 'Value',
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
    instance of the PhenList class
    mandatory argument
  }
  \item{
    depVariable
  }{
    a character string defining the dependent variable of interest
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
      \item boxplotSexGenotype. See documentations for {
        {
          boxplotSexGenotype
        }
      }
      \item boxplotSexGenotypeBatchAdjusted. See documentations for {
        {
          boxplotSexGenotypeBatchAdjusted
        }
      }
      \item boxplotSexGenotypeWeightBatchAdjusted. See documentations for {
        {
          boxplotSexGenotypeWeightBatchAdjusted
        }
      }
      \item scatterplotSexGenotypeBatch. See documentations for {
        {
          scatterplotSexGenotypeBatch
        }
      }
      \item scatterplotGenotypeWeight. See documentations for {
        {
          scatterplotGenotypeWeight
        }
      }
      \item boxplotSexGenotypeBatch. See documentations for {
        {
          boxplotSexGenotypeBatch
        }
      }
    }
  }
 \item{\dots}{
  }
}
\author{
  Hamed Haselimashhadi < hamedhm@eb.ac.uk >
}

\seealso{
  \code{\link{PhenList}},
}

\examples{
  example(testDataset)
}
\value{
    Empty.
}
