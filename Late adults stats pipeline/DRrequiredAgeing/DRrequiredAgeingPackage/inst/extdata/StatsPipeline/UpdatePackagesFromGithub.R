library(devtools)
# DRrequiredAgeing
install_github(
  repo = "mpi2/impc_stats_pipeline/Late adults stats pipeline/DRrequiredAgeing/DRrequiredAgeingPackage",
  dependencies = FALSE,
  upgrade = "always",
  force = TRUE,
  build = TRUE,
  quiet = FALSE,
  ref = "dev"
)

# SmoothWin
install_github(
  repo = "mpi2/impc_stats_pipeline/SoftWindowing/SmoothWin/SmoothWinPackage/",
  dependencies = FALSE,
  upgrade = "always",
  force = TRUE,
  build = TRUE,
  quiet = FALSE,
  ref = "dev"
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
