# Recreating execution environment from scratch

```bash
# After becoming mi_stats

# 1. Deactivate existing environment
conda deactivate

# 2. Update conda
conda update -n base -c defaults conda

# 3. Remove environment
conda remove -n R2D2 --all

# 4. Create new environment
conda create --name R2D2

# 5. Activate new environment
conda activate R2D2

# 6. Install R and packages
conda install -c conda-forge r-base r-devtools libpq r-arrow r-mass r-rcppgsl r-magick r-matrix

# 7. Run main package installation script
export REMOTE="mpi2"
export BRANCH="dev"
wget -qO- https://raw.githubusercontent.com/${REMOTE}/impc_stats_pipeline/${BRANCH}/Late%20adults%20stats%20pipeline/DRrequiredAgeing/DRrequiredAgeingPackage/inst/extdata/StatsPipeline/UpdatePackagesFromGithub.R > UpdatePackagesFromGithub.R
Rscript UpdatePackagesFromGithub.R ${REMOTE} ${BRANCH} TRUE
rm UpdatePackagesFromGithub.R
```
