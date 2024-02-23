library(devtools)

# Check if the correct number of arguments is provided
if (length(commandArgs(trailingOnly = TRUE)) != 2) {
  stop("Please provide exactly two arguments: repository name and branch name.")
}

# Retrieve the command line arguments
repository <- commandArgs(trailingOnly = TRUE)[1]
branch <- commandArgs(trailingOnly = TRUE)[2]

# DRrequiredAgeing
install_github(
  repo = "marinak-ebi/impc_stats_pipeline/Late adults stats pipeline/DRrequiredAgeing/DRrequiredAgeingPackage",
  dependencies = FALSE,
  upgrade = "always",
  force = TRUE,
  build = TRUE,
  quiet = FALSE,
  ref = "debug-pipeline"
)

# SmoothWin
install_github(
  repo = "marinak-ebi/impc_stats_pipeline/SoftWindowing/SmoothWin/SmoothWinPackage/",
  dependencies = FALSE,
  upgrade = "always",
  force = TRUE,
  build = TRUE,
  quiet = FALSE,
  ref = "debug-pipeline"
)

# OpenStats
install_github(
  repo = "mpi2/OpenStats",
  dependencies = FALSE,
  upgrade = "always",
  force = TRUE,
  build = TRUE,
  quiet = FALSE
)
