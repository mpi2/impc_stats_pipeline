##########################################
# In this script the orders are important
##########################################
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
    source("http://bioconductor.org/biocLite.R")
    # biocLite(character(), ask=FALSE) #update installed packages.
    eval(parse(text = sprintf("biocLite(\"%s\",ask=FALSE)", x)))
    eval(parse(text = sprintf("require(\"%s\")", x)))
  }
}

##########################################
# Install the driver devtools package
#########################################
install.packages('devtools')
library(devtools)
source("https://bioconductor.org/biocLite.R")
biocLite("BiocInstaller", ask = FALSE)

##########################################
#### Install the proper version of some R packages
##########################################
install.packages(
  "https://cran.r-project.org/src/contrib/Archive/Rfast/Rfast_1.9.4.tar.gz",
  repos = NULL,
  type = "source"
)

install.packages(
  "https://cran.r-project.org/src/contrib/Archive/nloptr/nloptr_1.2.2.1.tar.gz",
  repos = NULL,
  type = "source"
)

install.packages(
  "https://cran.r-project.org/src/contrib/Archive/tidyr/tidyr_1.0.2.tar.gz",
  repos = NULL,
  type = "source"
)

install.packages(
  "https://cran.r-project.org/src/contrib/Archive/magick/magick_2.0.tar.gz",
  repos = NULL,
  type = "source"
)
#####


##########################################
# Install packages that need to be installed from the github
##########################################
# devtools::install_github("tidyverse/tidyr")
devtools::install_github("eddelbuettel/rcppgsl")
devtools::install_github("cran/rjags")
devtools::install_github("cran/latticeExtra")
devtools::install_github("hannesmuehleisen/miniparquet")


##########################################
# install packages
##########################################
packages <- c(
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
  "graph"
)

for (package in packages) {
  install.packages.auto(package)
}

##########################################
# Key packages that must be updated
# even if they are in the list above
##########################################
# PhenStat
install_github(
  repo = "mpi2/impc_stats_pipeline/Early adults stats pipeline/PhenStat/PhenStatPackage/PhenStat",
  dependencies = TRUE,
  upgrade = "never",
  force = TRUE,
  build = TRUE,
  quiet = FALSE
)

# DRrequired
install_github(
  repo = "mpi2/impc_stats_pipeline/Early adults stats pipeline/DRrequired/DRrequiredPackage",
  dependencies = TRUE,
  upgrade = "never",
  force = TRUE,
  build = TRUE,
  quiet = FALSE
)

# DRrequiredAgeing
install_github(
  repo = "mpi2/impc_stats_pipeline/Late adults stats pipeline/DRrequiredAgeing/DRrequiredAgeingPackage",
  dependencies = TRUE,
  upgrade = "never",
  force = TRUE,
  build = TRUE,
  quiet = FALSE
)

# SmoothWin
install_github(
  repo = "mpi2/impc_stats_pipeline/SoftWindowing/SmoothWin/SmoothWinPackage/",
  dependencies = TRUE,
  upgrade = "never",
  force = TRUE,
  build = TRUE,
  quiet = FALSE
)

# OpenStats
install_github(
  # repo = 'mpi2/impc_stats_pipeline/Late adults stats pipeline/OpenStats/OpenStatsPackage/',
  repo = "mpi2/OpenStats",
  dependencies = TRUE,
  upgrade = "never",
  force = TRUE,
  build = TRUE,
  quiet = FALSE
)