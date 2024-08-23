#!/bin/bash

# Function prints messages to logs
function message0() {
    echo "$(date '+%Y-%m-%d %H:%M:%S.') $1"
}

# Statistical pipeline
export start_time=$(date '+%Y-%m-%d %H:%M:%S')
message0 "Starting the IMPC statistical pipeline..."
mkdir SP
export input_path=$(realpath .)
export sp_results=$(realpath SP)
message0 "Parquet files path: ${input_path}"
message0 "Output path: ${sp_results}"
cd SP

# Phase I: Prepare parquet files
message0 "Phase I. Convert parquet files into Rdata..."
message0 "Step 1. Create jobs"
step1_files=$(find .. -type f -name "*.parquet" -exec realpath {} \;)
for file in $step1_files; do
  echo "sbatch --job-name=impc_stats_pipeline_job --mem=10G --time=00:10:00 -e ${file}.err -o ${file}.log --wrap='Rscript Step2Parquet2Rdata.R $file'" >> jobs_step2_Parquet2Rdata.bch
done

# Calculate total execution time
end_time=$(date '+%Y-%m-%d %H:%M:%S')
start_seconds=$(date -j -f '%Y-%m-%d %H:%M:%S' "$start_time" '+%s')
end_seconds=$(date -j -f '%Y-%m-%d %H:%M:%S' "$end_time" '+%s')
duration_seconds=$((end_seconds - start_seconds))
duration_minutes=$(echo "scale=2; $duration_seconds / 60" | bc)
message0 "SP finished in ${duration_minutes}"
