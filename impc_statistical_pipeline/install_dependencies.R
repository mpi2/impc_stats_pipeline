##########################################
# In this script the orders are important
##########################################

#options(repos=structure(c(CRAN="YOUR FAVORITE MIRROR")))
options(repos = c(CRAN = "https://cloud.r-project.org/"))
update.packages(ask = FALSE);
##########################################
# Install the driver devtools package
#########################################
R_REMOTES_NO_ERRORS_FROM_WARNINGS="false"
options(warn=1)

if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages(
    "devtools",
    repos = "https://cloud.r-project.org/",
    dependencies = TRUE,
    quiet = TRUE
  )
} else{
  require(devtools)
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org/", quiet = TRUE)
}
