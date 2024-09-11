#!/bin/bash
set -e

# Assign arguments to variables.
VERSION="$1"
REMOTE="$2"
BRANCH="$3"
WINDOWING_PIPELINE="$4"

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
        sleep 5
    done
}

# Function to fetch a specific R script.
function fetch_script() {
    file_name=$(basename $1)
    wget -O ${file_name} --quiet "https://github.com/${REMOTE}/impc_stats_pipeline/raw/${BRANCH}/Late%20adults%20stats%20pipeline/DRrequiredAgeing/DRrequiredAgeingPackage/inst/extdata/StatsPipeline/$1"
}

# Statistical pipeline.
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
fetch_script 0-ETL/Step2Parquet2Rdata.R
sbatch --job-name=impc_stats_pipeline_job --time=01:00:00 --mem=1G -o ../compressed_logs/step2_job_id.txt --wrap="bash jobs_step2_Parquet2Rdata.bch"
waitTillCommandFinish
rm Step2Parquet2Rdata.R
find ../ -type f -name '*.log' -exec zip -q -m ../compressed_logs/step2_logs.zip {} +
find ../ -type f -name '*.err' -exec zip -q -m ../compressed_logs/step2_logs.zip {} +

message0 "Step 3. Merging pseudo Rdata files into single file for each procedure - jobs creator"
dirs=$(find "${sp_results}/ProcedureScatterRdata" -maxdepth 1 -type d)
for dir in $dirs; do
  echo "sbatch --job-name=impc_stats_pipeline_job --mem=50G --time=01:30:00 -e ${dir}/step4_merge_rdatas.err -o ${dir}/step4_merge_rdatas.log --wrap='Rscript Step4MergingRdataFiles.R ${dir}'" >> jobs_step4_MergeRdatas.bch
done

message0 "Step 4. Merging pseudo Rdata files into single files per procedure"
fetch_script 0-ETL/Step4MergingRdataFiles.R
sbatch --job-name=impc_stats_pipeline_job --time=01:00:00 --mem=1G -o ../compressed_logs/step4_job_id.txt --wrap="bash jobs_step4_MergeRdatas.bch"
waitTillCommandFinish
rm Step4MergingRdataFiles.R
find . -type f -name '*.log' -exec zip -q -m ../compressed_logs/step4_logs.zip {} +
find . -type f -name '*.err' -exec zip -q -m ../compressed_logs/step4_logs.zip {} +

message0 "Phase I. Compressing the log files and house cleaning..."
zip -q -rm ../compressed_logs/phase1_jobs.zip *.bch
rm -rf ProcedureScatterRdata

message0 "Starting Phase II, packaging the big data into small packages ..."
mkdir DataGeneratingLog
for file in $(find Rdata -type f -exec realpath {} \;); do
  file_basename=$(basename $file .Rdata)
  echo "sbatch --job-name=impc_stats_pipeline_job --mem=45G --time=6-00 -e DataGeneratingLog/${file_basename}_errorlog.log -o DataGeneratingLog/${file_basename}_outputlog.log --wrap='Rscript InputDataGenerator.R ${file} ${file_basename}'" >> DataGenerationJobList.bch
done
fetch_script jobs/InputDataGenerator.R
sbatch --job-name=impc_stats_pipeline_job --time=01:00:00 --mem=1G -o ../compressed_logs/phase2_job_id.txt --wrap="bash DataGenerationJobList.bch"
waitTillCommandFinish
rm InputDataGenerator.R

message0 "End of packaging data."
message0 "Phase II. Compressing the log files and house cleaning..."
mv *.bch  DataGeneratingLog/
zip -q -rm phase2_logs.zip DataGeneratingLog/
mv phase2_logs.zip ../compressed_logs/

message0 "Appending all procedure based jobs into one single file..."
mkdir jobs
find ./*/*_RawData/*.bch -type f | xargs  cat >> jobs/AllJobs.bch

message0 "Phase III. Initialising the statistical analysis..."
cd jobs
message0 "Updating the dynamic contents from the IMPReSS..."
R --quiet -e \
"DRrequiredAgeing:::updateImpress( \
  updateImpressFileInThePackage = TRUE, \
  updateOptionalParametersList = TRUE, \
  updateTheSkipList = TRUE, \
  saveRdata = FALSE \
)"

message0 "Running the IMPC statistical pipeline by submitting jobs..."
if [ "${WINDOWING_PIPELINE}" = true ]; then
  fetch_script jobs/function_windowed.R
  mv function_windowed.R function.R
else
  fetch_script jobs/function.R
fi

R --quiet -e \
"DRrequiredAgeing:::ReplaceWordInFile( \
  '$(realpath function.R)', \
  'DRversionNotSpecified', \
  ${VERSION} \
)"
chmod 775 AllJobs.bch
