\name{analysisResults}
\alias{analysisResults}
\title{Method "analysisResults"}
\docType{methods}
\description{
    Method \code{analysisResults} returns analysis results
}
\section{Methods}{
    \describe{
        \item{\code{signature(obj = "PhenTestResult")}}{
            %%  Returns results of the analysis
        }
    }}
\keyword{methods}

\examples{
    # Mixed Model framework
    file <- system.file("extdata", "test1.csv", package="PhenStat")
    test <- PhenStat:::PhenList(dataset=read.csv(file,na.strings = '-'),
                                testGenotype="Sparc/Sparc")

    result <- PhenStat:::testDataset(test,
                                     depVariable="Lean.Mass")
    r = PhenStat:::analysisResults(obj = result)
}

\value{
    The function returns the analysis result.
}
