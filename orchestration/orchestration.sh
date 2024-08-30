#!/bin/bash

# Assign arguments to variables.
VERSION="$1"
REMOTE="$2"
BRANCH="$3"

# Function prints messages to logs.
function message0() {
    echo "$(date '+%Y-%m-%d %H:%M:%S.') $1"
}

# Function to wait until all jobs on SLURM complete.
function waitTillCommandFinish() {
    while true; do
        if ! (squeue --format="%A %.30j" | grep -q "impc_stats_pipeline_job"); then
            message0 "Done waiting for SLURM jobs to complete."
            break
        fi
        sleep 60
    done
}

# Statistical pipeline.
export start_time=$(date '+%Y-%m-%d %H:%M:%S')
message0 "Starting the IMPC statistical pipeline..."
mkdir SP compressed_logs
export input_path=$(realpath .)
export sp_results=$(realpath SP)
message0 "Parquet files path: ${input_path}"
message0 "Output path: ${sp_results}"
cd SP

message0 "Phase I. Convert parquet files into Rdata..."

message0 "Step 1. Create jobs"
step1_files=$(find .. -type f -name '*.parquet' -exec realpath {} \;)
for file in $step1_files; do
  echo "sbatch --job-name=impc_stats_pipeline_job --mem=10G --time=00:10:00 -e ${file}.err -o ${file}.log --wrap='Rscript Step2Parquet2Rdata.R $file'" >> jobs_step2_Parquet2Rdata.bch
done

message0 "Step 2. Read parquet files and create pseudo Rdata"
wget --quiet "https://github.com/${REMOTE}/impc_stats_pipeline/raw/${BRANCH}/Late%20adults%20stats%20pipeline/DRrequiredAgeing/DRrequiredAgeingPackage/inst/extdata/StatsPipeline/0-ETL/Step2Parquet2Rdata.R"
sbatch --job-name=impc_stats_pipeline_job --time=01:00:00 --mem=1G -o ../compressed_logs/step2_job_id.txt --wrap="bash jobs_step2_Parquet2Rdata.bch"
waitTillCommandFinish
rm Step2Parquet2Rdata.R
find ../ -type f -name '*.log' -exec zip -m ../compressed_logs/step2_logs.zip {} +
find ../ -type f -name '*.err' -exec zip -m ../compressed_logs/step2_logs.zip {} +

message0 "Step 3. Merging pseudo Rdata files into single file for each procedure - jobs creator"
dirs=$(find "${sp_results}/ProcedureScatterRdata" -maxdepth 1 -type d)
for dir in $dirs; do
  echo "sbatch --job-name=impc_stats_pipeline_job --mem=50G --time=01:30:00 -e ${dir}/step4_merge_rdatas.err -o ${dir}/step4_merge_rdatas.log --wrap='Rscript Step4MergingRdataFiles.R ${dir}'" >> jobs_step4_MergeRdatas.bch
done

message0 "Step 4. Merging psudo Rdata files into single files per procedure"
wget --quiet "https://github.com/${REMOTE}/impc_stats_pipeline/raw/${BRANCH}/Late%20adults%20stats%20pipeline/DRrequiredAgeing/DRrequiredAgeingPackage/inst/extdata/StatsPipeline/0-ETL/Step4MergingRdataFiles.R"
sbatch --job-name=impc_stats_pipeline_job --time=01:00:00 --mem=1G -o ../compressed_logs/step4_job_id.txt --wrap="bash jobs_step4_MergeRdatas.bch"
waitTillCommandFinish
rm Step4MergingRdataFiles.R
find . -type f -name '*.log' -exec zip -m ../compressed_logs/step4_logs.zip {} +
find . -type f -name '*.err' -exec zip -m ../compressed_logs/step4_logs.zip {} +






# Calculate total execution time.
end_time=$(date '+%Y-%m-%d %H:%M:%S')
start_seconds=$(date -j -f '%Y-%m-%d %H:%M:%S' "$start_time" '+%s')
end_seconds=$(date -j -f '%Y-%m-%d %H:%M:%S' "$end_time" '+%s')
duration_seconds=$((end_seconds - start_seconds))
duration_minutes=$(echo "scale=2; $duration_seconds / 60" | bc)
message0 "SP finished in ${duration_minutes}"
