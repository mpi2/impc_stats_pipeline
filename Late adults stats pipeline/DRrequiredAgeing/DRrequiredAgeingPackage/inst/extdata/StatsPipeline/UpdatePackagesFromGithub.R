library(devtools)
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
