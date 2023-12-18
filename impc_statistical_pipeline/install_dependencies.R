##########################################
# In this script the orders are important
##########################################

#options(repos=structure(c(CRAN="YOUR FAVORITE MIRROR")))
options(repos = c(CRAN = "https://cran.rstudio.com/"))
update.packages(ask = FALSE);
##########################################
# Install the driver devtools package
#########################################
R_REMOTES_NO_ERRORS_FROM_WARNINGS="false"
options(warn=1)

if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages(
    "devtools",
    repos = "https://cran.rstudio.com/"
  )
}
require("devtools")
require("remotes")


if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cran.rstudio.com/", quiet = FALSE)
}


# Installer script

install.packages.auto <- function(x, v) {
  message('\n\n\t-> Installing package: ', x, ' version = ', v)
  avpackages = rownames(available.packages())
  ps = installed.packages()
  exist = FALSE
  for (p in ps[, 1]) {
    if (p == x) {
      if (packageVersion(p) == v || 1 == 1) {
        exist = TRUE
        message('\t --> package is already installed.')
        break
      }
    }
  }

  if (!exist) {
    if (x %in% avpackages) {
      tryCatch({
        remotes::install_version(
          package = x,
          version = v,
          repos = "https://cran.rstudio.com/",
          quiet = FALSE,
          force = FALSE,
          upgrade = 'never',
          dependencies = TRUE,
          type = "source"
        )
        # devtools::install_version(
        #   "lmerTest",
        #   version = "3.1-2",
        #   repo = "lib.ugent.be/CRAN",
        #   dependencies = TRUE,
        #   type = "source"
        # )
      },
      error = function(cond) {
        message('****** --->', cond)
        return(NA)
      },
      warning = function(cond) {
        message('****** --->', cond)
      })
    } else{
      eval(parse(text = sprintf(
        "BiocManager::install(\"%s\",ask=FALSE)", x
      )))
    }
  }

  #require(x)
}


##########################################
# Install packages that need to be installed from github
##########################################
if (!requireNamespace("data.table", quietly = TRUE)) {
  devtools::install_github("Rdatatable/data.table", upgrade = "never",force = FALSE,dependencies =TRUE)
}

if (!requireNamespace("latticeExtra", quietly = TRUE)) {
  devtools::install_github("cran/latticeExtra",upgrade = "never",force = FALSE,dependencies = TRUE)
}

if (!requireNamespace("miniparquet", quietly = TRUE)) {
  devtools::install_github("hannesmuehleisen/miniparquet",upgrade = "never",force = FALSE,dependencies = TRUE)
}

if (!requireNamespace("Hmisc", quietly = TRUE)) {
  install.packages('Hmisc', repos = "https://cran.rstudio.com/")
}
if (!requireNamespace("Cairo", quietly = TRUE)) {
  install.packages('Cairo', repos = "https://cran.rstudio.com/")
}
if (!requireNamespace("rwebhdfs", quietly = TRUE)) {
  devtools::install_github(c("saurfang/rwebhdfs"),upgrade = 'always',force = TRUE)
}

##########################################
# install packages
##########################################
packages <- c(
 "littler"      ,"0.3.12"  ,
 "RcppGSL"      ,"0.3.8"  ,
 "AICcmodavg"   ,"2.3.1"  ,
 "nloptr"      , "1.2.2.2",
 "car"         , "3.0.10" ,
 "RJSONIO"     , "1.3.1.4",
 "base64enc"   , "0.1.3"  ,
 "doParallel"  , "1.0.16" ,
 "parallel"    , "4.0.2"  ,
 "abind"       , "1.4.5"  ,
 "DBI"         , "1.1.1"  ,
 "plyr"        , "1.8.6"  ,
 "nortest"     , "1.0.4"  ,
 "pingr"       , "2.0.1"  ,
 "RPostgreSQL" , "0.6.2"  ,
 "quantreg"    , "5.82"   ,
 "car"        ,  "3.0.10" ,
 "RcppZiggurat", "0.1.6"  ,
 "tidyr"      ,  "1.1.2"  ,
 "methods"    ,  "4.0.2"  ,
 "jsonlite"   ,  "1.7.2"  ,
 "foreach"    ,  "1.5.1"  ,
 "MASS"       ,  "7.3.53" ,
 "survival"   ,  "3.2.7"  ,
 "RSQLite"    ,  "2.2.2"  ,
 "robustbase" ,  "0.93.7" ,
 "msgps"      ,  "1.3.1"  ,
 "corrplot"   ,  "0.84"   ,
 "Tmisc"      ,  "1.0.0"  ,
 "Hmisc"       , "4.4.1"  ,
 "summarytools", "0.9.8"  ,
 "lme4"        , "1.1.26" ,
 "stringi"     , "1.5.3"  ,
 "pingr"       , "2.0.1"  ,
 "nlme"        , "3.1.151",
 "base"        , "4.0.2"  ,
 "rlist"       , "0.4.6.1",
 "gtools"      , "3.8.2"  ,
 "rlang"       , "0.4.10" ,
 "logistf"     , "1.24"   ,
 "graph"       , "1.66.0" ,
 "digest"      , "0.6.27" ,
 "magick"      , "2.0"    ,
 "Rfast"       , "1.9.4"  ,
 "nloptr"      , "1.2.2.1",
 "tidyr"       , "1.0.2"
)

for (i in seq(1,length(packages),by=2)) {
  install.packages.auto(packages[i],packages[i+1])
}

##########################################
# Key packages that must be updated
# even if they are in the list above
##########################################
# SmoothWin
devtools::install_github(
  repo = "mpi2/impc_stats_pipeline/SoftWindowing/SmoothWin/SmoothWinPackage/",
  dependencies = TRUE,
  upgrade = "never",
  force = TRUE,
  build = TRUE,
  quiet = FALSE
)

# OpenStats
devtools::install_github(
  # repo = 'mpi2/impc_stats_pipeline/Late adults stats pipeline/OpenStats/OpenStatsPackage/',
  repo = "mpi2/OpenStats",
  dependencies = TRUE,
  upgrade = "never",
  force = TRUE,
  build = TRUE,
  quiet = FALSE
)

# DRrequired
devtools::install_github(
  repo = "mpi2/impc_stats_pipeline/Early adults stats pipeline/DRrequired/DRrequiredPackage",
  dependencies = TRUE,
  upgrade = "never",
  force = TRUE,
  build = TRUE,
  quiet = FALSE
)

# DRrequiredAgeing
devtools::install_github(
  repo = "mpi2/impc_stats_pipeline/Late adults stats pipeline/DRrequiredAgeing/DRrequiredAgeingPackage",
  dependencies = TRUE,
  upgrade = "never",
  force = TRUE,
  build = TRUE,
  quiet = FALSE
)

# testing for errors
pts = c(
  packages[seq(1, length(packages), by = 2)],
  'SmoothWin',
  'OpenStats',
  'DRrequired',
  'DRrequiredAgeing',
  'latticeExtra',
  'Hmisc',
  'miniparquet',
  'data.table'
)
for(pt in pts) {
  if (!requireNamespace(pt, quietly = TRUE))
    message('oh no, the package ', pt, ' does not exists!')
}
