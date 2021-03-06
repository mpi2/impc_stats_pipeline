\name{testDataset}
\alias{testDataset}
\title{Method "testDataset"}
\description{
    Statistical analysis manager function in PhenStat package.

    Firstly, it performs the checks for dependent variable in the dataset.
    Some checks are generic, some depends on selected statistical framework/method.

    Secondly, if checks are passed it runs the stages of selected framework to analyse dependent variable.

    For instance, runs startModel and finalModel for the MM framework if the argument "callAll" is set to TRUE.
    If framework contains only one stage like in Fisher Exact Test case then runs that one stage regardless "callAll" value.

}
\usage{
testDataset(
          phenList = NULL                ,
          depVariable = NULL             ,
          equation = "withWeight"        ,
          outputMessages = TRUE          ,
          pThreshold = 0.05              ,
          method = "MM"                  ,
          modelWeight = NULL             ,
          callAll = TRUE                 ,
          keepList = NULL                ,
          dataPointsThreshold = 4        ,
          RR_naturalVariation = 95       ,
          RR_controlPointsThreshold = 60 ,
          transformValues = FALSE        ,
          useUnfiltered = FALSE          ,
          threshold = 10 ^ -18           ,
          check = 1                      ,
          upper = 5
)
}
\arguments{
    \item{phenList}{instance of the \code{\linkS4class{PhenList}} class; mandatory argument}
    \item{depVariable}{a character string defining the dependent variable of interest; mandatory argument}
    \item{equation}{a character string defining the equation to use.
        Possible values "withWeight" (default), "withoutWeight"}
    \item{outputMessages}{flag: "FALSE" value to suppress output messages; "TRUE" value to show output messages;
        default value TRUE}
    \item{pThreshold}{a numerical value for the p-value threshold used to determine which fixed effects
        to keep in the model, default value 0.05}
    \item{method}{a character string ("MM", "FE", "TF" or "RR") defining the method to use for model building;
        default value "MM" for mixed model}
    \item{modelWeight}{
        a vector of possible weights (the same length as the response) that are imposed on the weighted MM. Only works when method="MM". The default is NULL that leads to normal (non-weighted) MM.
    }
    \item{callAll}{flag: "FALSE" value to run only the first stage of the selected framework; "TRUE" value (default)
        to run all stages of framework one after another}
    \item{keepList}{a logical vector defining the significance of different model effects: keep_batch, keep_equalvar,
        keep_weight, keep_sex, keep_interaction; default value NULL}
    \item{dataPointsThreshold}{threshold for the number of data points in the MM framework; default value is 4
        ; minimal value is 2}
    \item{RR_naturalVariation}{threshold for the variation ranges in the RR framework; default value is 95
        ; minimal value is 60}
    \item{RR_controlPointsThreshold}{threshold for the number of control data points in the RR framework
        ; default value is 60; minimal value is 40}
    \item{transformValues}{flag: "FALSE" value to suppress transformation; "TRUE" value to perform transformation if needed;
        default value FALSE}
    \item{useUnfiltered}{flag: "FALSE" value to use filtered dataset; "TRUE" value to use unfiltered dataset;
        default value FALSE}
    \item{threshold}{The threshold for the model weights to be considered as zero. It only works when method="MM" and modelWeight is not null. Default is 10^-18.}
    \item{check}{Only useful when modelWeight is included. Select between 0, 1, 2 to impose different weighting strategies on the Linear Mixed model. check=1 (default) selects the weights that are greater than the threshold (above) and keeps the batches that include more than one single sample. check=2 keeps only the weights that are greater than the threshold but ignores the single sample batches. check=0 to disable the function.}
    \item{upper}{Single integer: maximum allowed number of levels for the Batch in the TF framework. Default 5}
}
\value{
    Returns results stored in instance of the \code{\linkS4class{PhenTestResult}} class
}

\references{
    Karp N, Melvin D, Sanger Mouse Genetics Project, Mott R (2012): Robust and Sensitive Analysis of Mouse Knockout Phenotypes. \emph{PLoS ONE} \bold{7}(12): e52410. doi:10.1371/journal.pone.0052410

    West B, Welch K, Galecki A (2007): Linear Mixed Models: A practical guide using statistical software  \emph{New York: Chapman & Hall/CRC} 353 p.
}
\author{Natalja Kurbatova, Natasha Karp, Jeremy Mason}
\seealso{ \code{\linkS4class{PhenList}},
          \code{\link{plot.PhenTestResult}},
          \code{\link{plot.PhenList}}
          }
\examples{
    # Mixed Model framework
    file <- system.file("extdata", "test1.csv", package="PhenStat")
    test <- PhenStat:::PhenList(dataset=read.csv(file,na.strings = '-'),
            testGenotype="Sparc/Sparc")
    plot(test,depVariable="Lean.Mass")

    result <- PhenStat:::testDataset(test,
            depVariable="Lean.Mass")
    # print out formula that has been created
    PhenStat:::analysisResults(result)$model.formula.genotype
    summary(result)
    plot(result)






    # Mixed Model framework with user defined effects
    user_defined_effects <- c(keep_batch=TRUE,
            keep_equalvar=TRUE,
            keep_weight=TRUE,
            keep_sex=TRUE,
            keep_interaction=TRUE)
    result3 <- PhenStat:::testDataset(test,
            depVariable="Lean.Mass",
            keepList=user_defined_effects)
    # print out formula that has been created
    PhenStat:::analysisResults(result3)$model.formula.genotype
    summary(result3)
    plot(result3)




    # Fisher Exact Test framework
    file <- system.file("extdata", "test_categorical.csv", package="PhenStat")
    test2 <- PhenStat:::PhenList(dataset=read.csv(file,na.strings = '-'),
            testGenotype="Aff3/Aff3")
    plot(test2,depVariable="Thoracic.Processes")
    result2 <- PhenStat:::testDataset(test2,
            depVariable="Thoracic.Processes",
            method="FE")
    summary(result2)
    plot(result2)
}
