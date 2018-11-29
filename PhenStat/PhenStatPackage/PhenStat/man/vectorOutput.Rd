\name{vectorOutput}
\alias{vectorOutput}
\title{Method "vectorOutput"}
\description{
    Wrapper for the model fitting results. Returns model fitting and testing results in a vector format.
    Assumes that modelling results are stored in the \code{\linkS4class{PhenTestResult}} object
    (output from function \code{\link{testDataset}}).
}
\usage{
vectorOutput(
	phenTestResult,
	phenotypeThreshold = 0.01,
	othercolumns = NULL,
	quote  = '"'
)
}
\arguments{
    \item{phenTestResult}{instance of the \code{\linkS4class{PhenTestResult}} class that comes from the function \code{\link{testDataset}}; mandatory argument}
    \item{phenotypeThreshold}{a numerical value defining the threshold to use in classificationTag; default value 0.01}
    \item{othercolumns}{a vector of column names that needs to be included in the function output}
    \item{quote}{The charachter that is used to encompass the JSON keys in summary statistics. The default is quotation}
}
\value{
    Returns a vector with model fitting results: model output summary and other values
}

\references{
    Karp N, Melvin D, Sanger Mouse Genetics Project, Mott R (2012): Robust and Sensitive Analysis of Mouse Knockout Phenotypes. \emph{PLoS ONE} \bold{7}(12): e52410. doi:10.1371/journal.pone.0052410

    West B, Welch K, Galecki A (2007): Linear Mixed Models: A practical guide using statistical software  \emph{New York: Chapman & Hall/CRC} 353 p.
}
\author{Natalja Kurbatova, Natasha Karp, Jeremy Mason}
\seealso{\code{\linkS4class{PhenTestResult}}}
\examples{
    # Mixed Model framework
    file <- system.file("extdata", "test1.csv", package="PhenStat")
    test <- PhenStat:::PhenList(dataset=read.csv(file,na.strings = '-'),
            testGenotype="Sparc/Sparc")
    result <- PhenStat:::testDataset(test,
            depVariable="Lean.Mass")
    PhenStat:::vectorOutput(result)

    # Fishe Exact Test framework
    file <- system.file("extdata", "test_categorical.csv", package="PhenStat")
    test2 <- PhenStat:::PhenList(dataset=read.csv(file,na.strings = '-'),
            testGenotype="Aff3/Aff3")
    result2 <- PhenStat:::testDataset(test2,
            depVariable="Thoracic.Processes",
            method="FE")
    PhenStat:::vectorOutput(result2)
}
