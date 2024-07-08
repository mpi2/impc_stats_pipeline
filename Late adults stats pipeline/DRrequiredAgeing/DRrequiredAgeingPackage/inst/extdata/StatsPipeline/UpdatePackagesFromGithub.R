library(devtools)

# Check if the correct number of arguments is provided
if (length(commandArgs(trailingOnly = TRUE)) != 2) {
  stop("Please provide exactly two arguments: repository name and branch name.")
}

# Retrieve the command line arguments
repository <- commandArgs(trailingOnly = TRUE)[1]
branch <- commandArgs(trailingOnly = TRUE)[2]

# Update DRrequiredAgeing
install_github(
  repo = paste(repository,
               "/impc_stats_pipeline/Late adults stats pipeline/DRrequiredAgeing/DRrequiredAgeingPackage",
               sep = ""),
  dependencies = TRUE,
  upgrade = "always",
  force = TRUE,
  build = TRUE,
  quiet = FALSE,
  ref = branch
)

# Update SmoothWin
install_github(
  repo = paste(repository,
               "/impc_stats_pipeline/SoftWindowing/SmoothWin/SmoothWinPackage",
               sep = ""),
  dependencies = TRUE,
  upgrade = "always",
  force = TRUE,
  build = TRUE,
  quiet = FALSE,
  ref = branch
)

# Update OpenStats
install_github(
  repo = "mpi2/OpenStats",
  dependencies = TRUE,
  upgrade = "always",
  force = TRUE,
  build = TRUE,
  quiet = FALSE
)
