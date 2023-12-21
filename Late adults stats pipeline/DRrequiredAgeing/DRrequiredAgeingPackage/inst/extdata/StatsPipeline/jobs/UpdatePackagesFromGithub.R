library(devtools)
# SmoothWin
install_github(
	repo = 'mpi2/impc_stats_pipeline/SoftWindowing/SmoothWin/SmoothWinPackage/',
	dependencies = FALSE,
	upgrade = 'always',
	force = TRUE,
	build = TRUE,
	quiet = FALSE
)

# OpenStats
install_github(
	#repo = 'mpi2/impc_stats_pipeline/Late adults stats pipeline/OpenStats/OpenStatsPackage/',
	repo = 'mpi2/OpenStats',
	dependencies = FALSE,
	upgrade = 'always',
	force = TRUE,
	build = TRUE,
	quiet = FALSE
)


# DRrequired
install_github(
	repo = 'mpi2/impc_stats_pipeline/Early adults stats pipeline/DRrequired/DRrequiredPackage',
	dependencies = FALSE,
	upgrade = 'always',
	force = TRUE,
	build = TRUE,
	quiet = FALSE
)

# DRrequiredAgeing
install_github(
	repo = 'mpi2/impc_stats_pipeline/Late adults stats pipeline/DRrequiredAgeing/DRrequiredAgeingPackage',
	dependencies = FALSE,
	upgrade = 'always',
	force = TRUE,
	build = TRUE,
	quiet = FALSE
)


