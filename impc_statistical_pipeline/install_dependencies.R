##########################################
# In this script the orders are important
##########################################
options(repos = c(CRAN = "https://cloud.r-project.org/"))
# Installer script
install.packages.auto <- function(x) {
  x <- as.character(x)
  if (isTRUE(x %in% .packages(all.available = TRUE))) {
    # eval(parse(text = sprintf("require(\"%s\")", x)))
  } else {
    # update.packages(ask= FALSE) #update installed packages.
    eval(parse(
      text = sprintf(
        "install.packages(\"%s\", dependencies = TRUE,repos='https://cloud.r-project.org')",
        x
      )
    ))
  }
  if (isTRUE(x %in% .packages(all.available = TRUE))) {
    # eval(parse(text = sprintf("require(\"%s\")", x)))
  } else {
    # Is bioconductor installed?
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }

    eval(parse(text = sprintf(
      "BiocManager::install(\"%s\",ask=FALSE)", x
    )))
    eval(parse(text = sprintf("require(\"%s\")", x)))
  }
}

##########################################
# Install the driver devtools package
#########################################
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}else{
  require(devtools)
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}


##########################################
#### Install the proper version of some R packages
##########################################
if (!requireNamespace("Rfast", quietly = TRUE)) {
  install.packages(
    "https://cran.r-project.org/src/contrib/Archive/Rfast/Rfast_1.9.4.tar.gz",
    repos = NULL,
    type = "source"
  )
}

if (!requireNamespace("nloptr", quietly = TRUE)) {
  install.packages(
    "https://cran.r-project.org/src/contrib/Archive/nloptr/nloptr_1.2.2.1.tar.gz",
    repos = NULL,
    type = "source"
  )
}

if (!requireNamespace("tidyr", quietly = TRUE)) {
  install.packages(
    "https://cran.r-project.org/src/contrib/Archive/tidyr/tidyr_1.0.2.tar.gz",
    repos = NULL,
    type = "source"
  )
}

if (!requireNamespace("magick", quietly = TRUE)) {
  install.packages(
    "https://cran.r-project.org/src/contrib/Archive/magick/magick_2.0.tar.gz",
    repos = NULL,
    type = "source"
  )
}
#####


##########################################
# Install packages that need to be installed from the github
##########################################
if (!requireNamespace("data.table", quietly = TRUE)) {
  devtools::install_github("Rdatatable/data.table", upgrade = "never")
}


if (!requireNamespace("rcppgsl", quietly = TRUE)) {
  devtools::install_github("eddelbuettel/rcppgsl")
}

# devtools::install_github("cran/rjags")
if (!requireNamespace("latticeExtra", quietly = TRUE)) {
  devtools::install_github("cran/latticeExtra")
}

if (!requireNamespace("miniparquet", quietly = TRUE)) {
  devtools::install_github("hannesmuehleisen/miniparquet")
}


##########################################
# install packages
##########################################
packages <- c(
  # "data.table",
  "RcppGSL",
  "quantreg",
  "Hmisc",
  "AICcmodavg",
  "car",
  "summarytools",
  "nloptr",
  "RcppZiggurat",
  "lme4",
  # "Rfast",
  "car",
  "tidyselect",
  "tidyr",
  "PhenStat",
  "RJSONIO",
  "methods",
  # "SmoothWin",
  "stringi",
  "base64enc",
  "jsonlite",
  "pingr",
  "doParallel",
  "foreach",
  "nlme",
  "parallel",
  "MASS",
  "base",
  "abind",
  "OpenStats",
  "rlist",
  "DBI",
  "RSQLite",
  "gtools",
  "plyr",
  "robustbase",
  "rlang",
  "nortest",
  "msgps",
  "logistf",
  "pingr",
  "corrplot",
  "graph",
  "RPostgreSQL",
  "Tmisc",
  "digest"
)

for (package in packages) {
  install.packages.auto(package)
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


