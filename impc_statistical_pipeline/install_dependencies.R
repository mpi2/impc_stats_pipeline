##########################################
# In this script the orders are important
##########################################

##########################################
# Install the driver devtools package
#########################################

if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools",repos = "https://cloud.r-project.org")
} else{
  require(devtools)
}

R_REMOTES_NO_ERRORS_FROM_WARNINGS="false"
options(warn=1)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager",repos = "https://cloud.r-project.org")
}


# Installer script

install.packages.auto <- function(x, v) {
  message('\t-> Installing package: ', x, ' version = ', v)
  avpackages = rownames(available.packages())
  ps = installed.packages()
  exist = FALSE
  for (p in ps[, 1]) {
    if (p == x) {
      if (packageVersion(p) == v || 1 == 1)
        exist = TRUE
    }
  }
  
  if (!exist) {
    if (x %in% avpackages) {
      tryCatch({
        remotes::install_version(
          package = x,
          version = v,
          repos = "https://cloud.r-project.org",
          quiet = TRUE,
          upgrade = 'never'
        )
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
  
  require(x)
}


##########################################
# Install packages that need to be installed from github
##########################################
if (!requireNamespace("data.table", quietly = TRUE)) {
  devtools::install_github("Rdatatable/data.table", upgrade = "never")
}

if (!requireNamespace("rcppgsl", quietly = TRUE)) {
  devtools::install_github("eddelbuettel/rcppgsl",upgrade = "never")
}

if (!requireNamespace("latticeExtra", quietly = TRUE)) {
  devtools::install_github("cran/latticeExtra",upgrade = "never")
}

if (!requireNamespace("miniparquet", quietly = TRUE)) {
  devtools::install_github("hannesmuehleisen/miniparquet",upgrade = "never")
}

if (!requireNamespace("Hmisc", quietly = TRUE)) {
  install.packages('Hmisc', repos = "https://cloud.r-project.org")
}

##########################################
# install packages
##########################################
packages <- c(
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
 "PhenStat"    , "2.18.0" ,
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

# PhenStat
devtools::install_github(
  repo = "mpi2/impc_stats_pipeline/Early adults stats pipeline/PhenStat/PhenStatPackage/PhenStat",
  dependencies = FALSE,
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
