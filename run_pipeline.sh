#!/bin/bash
set -ex

# Check if the correct number of arguments are provided
if [ "$#" -ne 6 ]; then
    echo "Usage: $0 <VERSION> <REMOTE>"
    exit 1
fi

# Assign arguments to variables
VERSION="$1"
REMOTE="$2"
BRANCH="$3"
KOMP_PATH="$4"
PARQUET_FOLDER="$5"
MP_CHOOSER_FOLDER="$6"

# Create a working directory
mkdir --mode=775 ${KOMP_PATH}/impc_statistical_pipeline/IMPC_DRs/stats_pipeline_input_dr${VERSION}
cd ${KOMP_PATH}/impc_statistical_pipeline/IMPC_DRs/stats_pipeline_input_dr${VERSION}

# Copy the input parquet files and mp_chooser_json
cp ${PARQUET_FOLDER}/*.parquet ./
cp ${MP_CHOOSER_FOLDER}/part*.txt ./mp_chooser.json

# Convert the mp_chooser JSON file to Rdata
R -e "a = jsonlite::fromJSON('mp_chooser.json');save(a,file='mp_chooser.json.Rdata')"
export MP_CHOOSER_FILE=$(echo -n '"'; realpath mp_chooser.json.Rdata | tr -d '\n'; echo -n '"')

# Update packages to the latest version
echo "Update started"
cd ${KOMP_PATH}/impc_statistical_pipeline/IMPC_DRs/stats_pipeline_input_dr${VERSION}
wget https://raw.githubusercontent.com/${REMOTE}/impc_stats_pipeline/${BRANCH}/Late%20adults%20stats%20pipeline/DRrequiredAgeing/DRrequiredAgeingPackage/inst/extdata/StatsPipeline/UpdatePackagesFromGithub.R
Rscript UpdatePackagesFromGithub.R ${REMOTE} ${BRANCH} FALSE
rm UpdatePackagesFromGithub.R
echo "Update completed"

# Run Statistical Pipeline
job1_txt=$(sbatch \
    --job-name=statistical_pipeline \
    --time=30-00:00:00 \
    --mem=8G \
    -o ../stats_pipeline_logs/stats_pipeline_${VERSION}.log \
    -e ../stats_pipeline_logs/stats_pipeline_${VERSION}.err \
    --wrap="cd ${KOMP_PATH}/impc_statistical_pipeline/IMPC_DRs/stats_pipeline_input_dr${VERSION} && R -e 'DRrequiredAgeing:::StatsPipeline(DRversion=${VERSION})'")
job1_id=$(echo $job1_txt | cut -d" " -f4)

# Run annotation pipeline
sbatch \
    --job-name=annotation_pipeline \
    --dependency=afterok:"${job1_id}" \
    --time=3-00:00:00 \
    --mem=8G \
    -o ../../../../stats_pipeline_logs/annotation_pipeline_${VERSION}.log \
    -e ../../../../stats_pipeline_logs/annotation_pipeline_${VERSION}.err \
    --wrap="cd ${KOMP_PATH}/impc_statistical_pipeline/IMPC_DRs/stats_pipeline_input_dr${VERSION}/SP/jobs/Results_IMPC_SP_Windowed && R -e 'DRrequiredAgeing:::IMPC_HadoopLoad(prefix=${VERSION},transfer=FALSE,mp_chooser_file=${MP_CHOOSER_FILE})'"

