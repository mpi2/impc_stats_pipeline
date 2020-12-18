\name{mainAgeing}
\alias{mainAgeing}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The driver function for the International Mouse Phenotyping Consortium (IMPC) statistical pipeline
}
\description{
This is the driver function to execute the IMPC statistical pipeline. The function allows several input parameters to manage the different aspects of the statistical engine such as performing data preprocessing and packaging, statistical method, plots, simulations, optimization, outlier neutralization, soft windowing etc. \cr

The output of the pipeline is a dierctory structure as follows: centre_name/procedure_stable_id/parameter_stable_ir/colony_id/zygosity/metadata/<STATPACKET>|<RawData>.\cr

StatPacket refers to a tab-delimited file containing the experiment details and a JSON object containing the statistical results.\cr

The StatPackets come in three varieties: Successful, Not_processed and Failed. \cr

Succesful: Contains all statistical results\cr
Not_processed: For know reason that in explained in the file, this dataset is not processed by the IMPC statistical pipeline\cr
Failed: In this case, the statistical pipeline failed to process this data.
}
\usage{
mainAgeing(
  file = NULL                                    ,
  sep = ','                                      ,
  na.strings = 'NA'                              ,
  normalisedPhenlist = FALSE                     ,
  subdir = 'Results'                             ,
  seed = 123456                                  ,
  MethodOfReadingCategoricalCategories  = 'file' ,
  # `file`, `solr` or `update`
  OverwriteExistingFiles  = FALSE                ,
  ignoreSkipList          = FALSE                ,
  onlyFillNotExisitingResults = FALSE            ,
  WhiteListMethods  = NULL                       ,
  # Carefully use this option!
  # It can remove the entire result file (some colonies in a single output file)
  # Only for multicore
  activateMulticore       = TRUE       ,
  coreRatio = 5 / 14                   ,
  concurrentControlSelect = FALSE      ,
  MultiCoreErrorHandling  = 'pass'     ,
  inorder                 = FALSE      ,
  verbose                 = TRUE       ,
  # OpenStats
  MMOptimise              = c(1, 1, 1, 1, 1, 1) ,
  FERROptimise            = c(TRUE, TRUE)        ,
  FERRrep                 = 1500       ,
  equation                = 'auto'     ,
  # Only for simulations
  simulation              = FALSE      ,
  Simulation.iteration    = 1          ,
  # Only for windowing
  activeWindowing = TRUE               ,
  sensitivity = c(1, 1, 1, 0)          ,
  pvalThreshold = c(0, 0, 0, 0)        ,
  check = 2                            ,
  direction = c(1, 1)                  ,
  weightORthreshold = 'weight'         ,
  predFunction = function(m) {
    predict(m)
  }                                    ,
  residFunction = function(m) {
    resid(m)
  }                                    ,
  messages = TRUE                      ,
  threshold = sqrt(.Machine$double.eps) * 10,
  outdelim = '\t'                      ,
  debug = TRUE                         ,
  encode = FALSE                       ,
  noSpaceAllowed = TRUE                ,
  plotWindowing = TRUE                 ,
  storeplot = TRUE                     ,
  virtualDrive   = FALSE               ,
  checkNamesForMissingColNames = TRUE  ,
  # Just for the OpenStats manuscript
  measureExecutionTime   = FALSE       ,
  # Raw data
  storeRawData           = TRUE        ,
  compressRawData        = TRUE        ,
  writeOutputToDB        = FALSE       ,
  # Only for Batch generator
  BatchProducer          =  FALSE      ,
  cpu = 4                              ,
  memory = 9000                        ,
  nMax                   = 10000       ,
  ChunkSize              = 24          ,
  MinColoniesInChunks    = 32          ,
  controlSize            = 1500        ,
  ### Just for outlier detection
  outlierDetection       = FALSE       ,
  ### Just for Ageing Batch GEnerator  ,
  combineEAandLA        = FALSE        ,
  solrBaseURL           = 'http://hx-noah-74-10:8090' ,
  ### Just for debugging
  superDebug             = FALSE                      ,
  subBreakColumns        = NULL,
  extraBatchParameters   = NULL,
  ...
)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{
  The input data file. It can be parquet files from Apache Spark or any delimited file such as CSV, TSV etc
}
  \item{sep}{
  The delimiter for the input data file. The default is comma for the CSV files. Does not apply to the parquet input files.
}
  \item{na.strings}{
  The proxy for the NA in the input data. The default is set to "NA". See "na.string" parameter in the read.csv function.
}
  \item{normalisedPhenlist}{
  Setting this parameter to TRUE normalizes the dependent variable and the (body)Weight prior to applying the statistical methods. The normalization results in Z scores (value-mean of column)/sd of the column. The default FALSE
}
  \item{subdir}{
  The directory of the output results. All results will be stored in this directory. Default "Results"
}
  \item{seed}{
  Seed to be able to reproduce the results. Default 123456
}
  \item{MethodOfReadingCategoricalCategories}{
  The method that reads the levels for the categorical data from IMPReSS. an example is "left eye missing", "right eye missing", "normal", "abnormal". allowed options are "file" and "update". Should always use the "file" for the high-throughput pipelines otherwise very time-consuming. You can always update the this file by running: `DRrequiredAgeing:::updateImpress(TRUE,TRUE,TRUE)`
}
  \item{OverwriteExistingFiles}{
  If set to TRUE then all resulting StatPackets will be overwritten if already exist. Default FALSE

}
  \item{ignoreSkipList}{
  Setting to TRUE ignores the skip list and the statistical methods will apply to all datasets. Skip list refers to a list of parameters that are flagged to skip from the statistical analysis (majority because annotation=FALSE in IMPReSS). Default FALSE. See the following link for the skip list: https://git.io/JIL24
}
  \item{onlyFillNotExisitingResults}{
  Setting to TRUE only applies the methods to the datasets where have no result already in the directory of results. This option is useful when the pipeline gets interrupted and you require to run the pipeline for the missing results. Default FALSE
}
  \item{WhiteListMethods}{
  Choose between "MM", "FE" and "RR". If set to any, then the other methods and their corresponding datasets will be ignored. Default NULL (all method will be applied)
}
  \item{activateMulticore}{
  Set to TRUE to activate multicore functioning of the pipeline. See "foreach" R package. Not recommended unless unlimited computing power/memory available. Default FALSE.
}
  \item{coreRatio}{
  A number in [0,1] interval. The core ratio for the multicore analysis. For example, if set to 0 then only 1 core will be used and if set to 1 then all cores for the computing machine will be used. The formula is ceiling(cores available on the machine*coreRatio). Default is 5/14.
}
  \item{concurrentControlSelect}{
  Only applies to the continues variables. Activate the concurrent control selection strategy. This will search for the same day controls as mutants. Min control required is 4 for each sex. Default FALSE.
}
  \item{MultiCoreErrorHandling}{
  See ".errorhandling" in the "foreach" R package. Options are "pass" and "stop". Default "pass"
}
  \item{inorder}{
  Does the order of processes need to be preserved? if set to TRUE then one analysis must be finish before the other one starts. See ".inorder" in "foreach" R package. Default FALSE
}
  \item{verbose}{
  Shows the progress of the multicore function. See "verbose" in the "foreach" R package. Default TRUE
}
  \item{MMOptimise}{
  See "MM_optimise" in the "OpenStatsAnalysis" function manual from "OpenStats" R package.
}
  \item{FERROptimise}{
  See "FERR_FullComparisions" in the "OpenStatsAnalysis" function manual from "OpenStats" R package.
}
  \item{FERRrep}{
  See "FERR_rep" in the "OpenStatsAnalysis" function manual from "OpenStats" R package.
}
  \item{equation}{
  Choose between "auto", "withWeight" and "withoutWeight". Setting to "auto" (default) will choose the best model for the IMPC data.
}
  \item{simulation}{
  Turn the pipeline into the simulation pipeline. That means that the mutants are replaced with the nearest controls and a random shift to left/right. This is used to create synthetic mutants and testing FDR. Default FALSE
}
  \item{Simulation.iteration}{
  The number of iterations in the simulation. Default 1 (must be kept 1 even if the simulation=FALSE)
}
  \item{activeWindowing}{
  Activates the SoftWindowing pipeline. This means that the output StatPackets will be contained the non-windowed and windowed results. Default FALSE
}
  \item{sensitivity}{
  See "sensitivity" parameter from the "SmoothWin" function in the "SmoothWin" R package.
}
  \item{pvalThreshold}{
  See "pvalThreshold" parameter from the "SmoothWin" function in the "SmoothWin" R package.
}
  \item{check}{
  See "check" parameter from the "SmoothWin" function in the "SmoothWin" R package.
}
  \item{direction}{
  See "direction" parameter from the "SmoothWin" function in the "SmoothWin" R package.
}
  \item{weightORthreshold}{
  See "weightORthreshold" parameter from the "SmoothWin" function in the "SmoothWin" R package.
}
  \item{predFunction}{
  See "predFunction" parameter from the "SmoothWin" function in the "SmoothWin" R package.
}
  \item{residFunction}{
  See "residFunction" parameter from the "SmoothWin" function in the "SmoothWin" R package.
}
  \item{messages}{
  See "messages" parameter from the "SmoothWin" function in the "SmoothWin" R package.
}
  \item{threshold}{
  See "threshold" parameter from the SmoothWin function in the SmoothWin R package.
}
  \item{outdelim}{
  Separator for the StatPacket fields. Default is set to tab "\\t"
}
  \item{debug}{
  Set to TRUE to see the details of the progress of the pipeline. Default TRUE.
}
  \item{encode}{
  Setting to TRUE compresses the output StatPackets (in-memory zip compression) and encodes the results to Base64 in the StatPacket output. Useful if disk space is a concern. Default FALSE
}
  \item{noSpaceAllowed}{
  If set to TRUE then any character outside "[^0-9A-Za-z/:]" range in the resulting directory structure (file structure in the outputs) will be replaced by underscore (_). Default TRUE
}
  \item{plotWindowing}{
  Set to TRUE to see the soft windowing plots. Default TRUE

}
  \item{storeplot}{
  Set to TRUE to store the soft windowing plots (pdf) in the output result directory. Note that this option only works on Linux/Unix. Default TRUE
}
  \item{virtualDrive}{
  Only applies to the Windows OS. Sometimes Windows OS does not allow long paths. Setting this option to TRUE will create a virtual directory that soft linked to the current working directory. It is normally a new virtual partition named X. Default FALSE
}
  \item{checkNamesForMissingColNames}{
  Equivalent to "check.names" in the "read.csv" function. Default TRUE
}
  \item{measureExecutionTime}{
  Setting this option to TRUE will export the execution time of the function in the procedure and parameter levels stored in the (/OpenStatsProcedureTime/a random name_procedure name) and (/OpenStatsParameterTime/a random name_parameter name). default FALSE
}
  \item{storeRawData}{
  Set to TRUE to store the raw input data in the output directory. Default TRUE
}
  \item{compressRawData}{
  Set to TRUE to store a compressed zip version of the raw data files. Default TRUE
}
  \item{writeOutputToDB}{
  Setting this option to TRUE will write the output of the statistical pipeline to the individual sqlitDB in the /db/ directory. This can stretch the execution time significantly. Default FALSE
}
  \item{BatchProducer}{
  This parameter activates the preprocessing step. In this step, no statistical analysis will be performed on the data but the data packages will be created. The big datasets will be breaking into small packages. The new directory structure with the data and the LSF jobs will be created. The resources such as CPU/memory allocated to the LSF jobs are specified by the parameters below. Default FALSE.
}
  \item{cpu}{
  For the preprocessing only. How many CPUs needed to be allocated for each LSF job? Default 4
}
  \item{memory}{
  For the preprocessing only. How much memory needed to be allocated for each LSF job? Default 9 GB
}
  \item{nMax}{
  For the preprocessing only. The maximum number of samples allowed in a package. Default 10K
}
  \item{ChunkSize}{
  For the preprocessing only. The number of colonies distributed in a data package? Formula =  (#colonies/ChunkSize). Default 24
}
  \item{MinColoniesInChunks}{
  If the total number of colonies in a packege exceed "MinColoniesInChunks", the break it into more packages. Default 32
}
  \item{controlSize}{
  If the total number of controls exceeds "controlSize" then make single colony jobs (one single colony per package). Default 1500
}
  \item{outlierDetection}{
Set to TRUE to activate the outlier detection algorithm. This algorithm removes the outliers in the controls. Default FALSE,
}
  \item{combineEAandLA}{
  Old parameter. ignore it
}
  \item{solrBaseURL}{
  Old parameter. ignore it
}
  \item{superDebug}{
  More detailed debug than the normal debug=TRUE. It exports a Rdata file in the output directory that contains all details of the process. Default FALSE
}
  \item{subBreakColumns}{
  ignore this parameter.
}
  \item{extraBatchParameters}{
  Old parameter. ignore it
}
  \item{\dots}{
  not applicable.
}
}

\value{
The output of this function is a set of files called StatPackets in a subdirectory specified by "subdir" parameter.
}
\references{
https://academic.oup.com/bioinformatics/article/36/5/1492/5583690 \cr
https://www.biorxiv.org/content/10.1101/2020.05.13.091157v1 \cr
http://bioconductor.org/packages/release/bioc/html/OpenStats.html \cr
https://cran.r-project.org/web/packages/SmoothWin/index.html \cr
}
\author{
Hamed Haseli Mashhadi <hamedhm@ebi.ac.uk>
}
\note{
a detailed description of the IMPC statistical pipeline is available from https://git.io/JILe4
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{UnzipAndfilePath}}
}
