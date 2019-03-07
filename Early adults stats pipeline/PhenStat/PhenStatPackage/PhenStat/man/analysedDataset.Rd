\name{analysedDataset}
\alias{analysedDataset}
\title{Method "analysedDataset"}
\docType{methods}
\description{
    Method \code{analysedDataset} returns subset of the original dataset
    that was analysed
}
\section{Methods}{
    \describe{
        \item{\code{signature(obj = "PhenTestResult")}}{
            %%  Returns analysed dataset
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
    head(PhenStat:::analysedDataset(obj = result))

}
\value{
    This function returns the analysed dataset.
}
