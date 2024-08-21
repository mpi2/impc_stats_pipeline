#!/bin/bash

# Function prints messages to logs
function message0() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') $1"
}

# Statistical pipeline
message0 "Starting the IMPC statistical pipeline..."
mkdir SP
export input_path=$(realpath .)
export sp_results=$(realpath SP)
message0 "Parquet files path: ${input_path}"
message0 "Output path: ${sp_results}"

cd SP
