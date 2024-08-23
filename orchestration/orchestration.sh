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

# Calculate total execution time
end_time=$(date '+%Y-%m-%d %H:%M:%S')
start_seconds=$(date -j -f '%Y-%m-%d %H:%M:%S' "$start_time" '+%s')
end_seconds=$(date -j -f '%Y-%m-%d %H:%M:%S' "$end_time" '+%s')
duration_seconds=$((end_seconds - start_seconds))
duration_minutes=$(echo "scale=2; $duration_seconds / 60" | bc)
message0 "SP finished in ${duration_minutes}"
