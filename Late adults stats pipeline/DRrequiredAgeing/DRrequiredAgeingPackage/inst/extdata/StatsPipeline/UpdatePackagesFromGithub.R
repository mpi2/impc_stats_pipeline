library(devtools)

# Check if the correct number of arguments is provided.
if (length(commandArgs(trailingOnly = TRUE)) != 3) {
  stop("Please provide exactly three arguments: repository name, branch name, and whether to update dependencies.")
}

# Retrieve the command line arguments.
repository <- commandArgs(trailingOnly = TRUE)[1]
branch <- commandArgs(trailingOnly = TRUE)[2]
dependencies <- commandArgs(trailingOnly = TRUE)[3]

# Update SmoothWin.
install_github(
  repo = paste(repository,
               "/impc_stats_pipeline/SoftWindowing/SmoothWin/SmoothWinPackage",
               sep = ""),
  dependencies = dependencies,
  upgrade = "always",
  force = TRUE,
  build = TRUE,
  quiet = FALSE,
  ref = branch
)

# Update OpenStats.
install_github(
  repo = "mpi2/OpenStats",
  dependencies = dependencies,
  upgrade = "always",
  force = TRUE,
  build = TRUE,
  quiet = FALSE
)

# Update DRrequiredAgeing.
install_github(
  repo = paste(repository,
               "/impc_stats_pipeline/Late adults stats pipeline/DRrequiredAgeing/DRrequiredAgeingPackage",
               sep = ""),
  dependencies = dependencies,
  upgrade = "always",
  force = TRUE,
  build = TRUE,
  quiet = FALSE,
  ref = branch
)
