#!/bin/bash
set -e

# Check if the correct number of arguments are provided
if [ "$#" -ne 6 ] && [ "$#" -ne 7 ]; then
    echo "Usage: $0 <VERSION> <REMOTE>"
    exit 1
fi

# Assign arguments to variables.
VERSION="$1"
REMOTE="$2"
BRANCH="$3"
KOMP_PATH="$4"
PARQUET_FOLDER="$5"
MP_CHOOSER_FOLDER="$6"
WINDOWING_PIPELINE=${7:-"true"}

# Redirect all output and errors to the log file
LOGFILE=${KOMP_PATH}/impc_statistical_pipeline/IMPC_DRs/stats_pipeline_logs/orchestration_${VERSION}.log
exec > >(tee -a "$LOGFILE") 2>&1

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
    folder_path=${2:-"Late%20adults%20stats%20pipeline/DRrequiredAgeing/DRrequiredAgeingPackage/inst/extdata/StatsPipeline"}
    wget -O ${file_name} --quiet "https://github.com/${REMOTE}/impc_stats_pipeline/raw/${BRANCH}/${folder_path}/$1"
}

# Function to submit limited number of jobs.
function submit_limit_jobs() {
    local bch_file=$1
    local job_id_logfile=$2
    local max_jobs=${3:-800}

    echo "Start submit_limit_jobs"
    echo > "$job_id_logfile"

    while IFS= read -r command; do
        while true; do
            num_running=$(squeue | wc -l)

            if [ "$num_running" -le "$max_jobs" ]; then
                break
            fi
            sleep 0.5
        done
        job_id=$(eval "$command")
        echo "$job_id" >> "$job_id_logfile"
    done < "$bch_file"

    echo "End submit_limit_jobs"
}

# Preparation.
mkdir --mode=775 ${KOMP_PATH}/impc_statistical_pipeline/IMPC_DRs/stats_pipeline_input_dr${VERSION}
cd ${KOMP_PATH}/impc_statistical_pipeline/IMPC_DRs/stats_pipeline_input_dr${VERSION}
cp ${PARQUET_FOLDER}/*.parquet ./
cp ${MP_CHOOSER_FOLDER}/part*.txt ./mp_chooser.json

message0 "Update packages to the latest version"
fetch_script UpdatePackagesFromGithub.R
Rscript UpdatePackagesFromGithub.R ${REMOTE} ${BRANCH} FALSE
rm UpdatePackagesFromGithub.R
message0 "Update completed"

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
submit_limit_jobs AllJobs.bch ../../compressed_logs/phase3_job_id.txt
waitTillCommandFinish

message0 "Postprocessing the IMPC statistical analysis results..."
cd ${sp_results}

message0 "Creating backups from the main R packages..."
mkdir RPackage_backup
R --quiet -e \
"DRrequiredAgeing:::packageBackup( \
  package='DRrequiredAgeing', \
  storepath='$(realpath RPackage_backup)' \
)"
R --quiet -e \
"DRrequiredAgeing:::packageBackup( \
  package='OpenStats', \
  storepath='$(realpath RPackage_backup)' \
)"

message0 "Compress phase III log files"
find . -type f -name '*.ClusterOut' -exec zip -m ../compressed_logs/phase3_logs.zip {} +
message0 "Compress phase III error files"
find . -type f -name '*.ClusterErr' -exec zip -m ../compressed_logs/phase3_errs.zip {} +

message0 "This is the last step. If you see no file in the list below, the SP is successfully completed."

# Annotation pipeline.
message0 "Starting the IMPC annotation pipeline..."
cd jobs/Results_IMPC_SP_Windowed/
message0 "Step 1: Clean ups and creating the global index for the results."
message0 "Indexing the results..."
for dir in $(find . -mindepth 2 -maxdepth 2 -type d); do
  base_dir=$(basename "$dir")
  output_file="FileIndex_${base_dir}_$(printf "%.6f" $(echo $RANDOM/32767 | bc -l)).Ind"
  echo "sbatch --job-name=impc_stats_pipeline_job --mem=1G --time=2-00 \
-e ${base_dir}_error.err -o ${base_dir}_output.log --wrap=\"find $dir -type f -name '*.tsv' > $output_file\"" >> minijobs.bch
done
chmod 775 minijobs.bch
submit_limit_jobs minijobs.bch ../../../compressed_logs/minijobs_job_id.txt
waitTillCommandFinish
mv minijobs.bch ../../../compressed_logs

find . -type f -name '*_output.log' -exec zip -m ../../../compressed_logs/minijobs_logs.zip {} +
find . -type f -name '*_error.err' -exec zip -m ../../../compressed_logs/minijobs_logs.zip {} +
message0 "Moving single indeces into a separate directory called AnnotationExtractorAndHadoopLoader..."
mkdir AnnotationExtractorAndHadoopLoader
chmod 775 AnnotationExtractorAndHadoopLoader
mv *.Ind AnnotationExtractorAndHadoopLoader
cd AnnotationExtractorAndHadoopLoader

message0 "Concatenating single index files to create a global index for the results..."
cat *.Ind >> AllResultsIndeces.txt
message0 "Zipping the single indeces..."
zip -rm allsingleindeces.zip *.Ind
split -50 AllResultsIndeces.txt split_index_

message0 "Convert the mp_chooser JSON file to Rdata..."
R --quiet -e "a = jsonlite::fromJSON('mp_chooser.json');save(a,file='mp_chooser.json.Rdata')"
export MP_CHOOSER_FILE=$(echo -n '"'; realpath mp_chooser.json.Rdata | tr -d '\n'; echo -n '"')

if [[ -z "${MP_CHOOSER_FILE}" || ! -f "${MP_CHOOSER_FILE}" ]]; then
    echo -e "ERROR: mp_chooser not found at location\n\t${MP_CHOOSER_FILE}"
    exit 1
fi

mkdir err log out
for file in $(find . -maxdepth 1 -type f -name "split_index*"); do
  echo "sbatch --job-name=impc_stats_pipeline_job --mem=5G --time=2-00 \
  -e err/$(basename "$file").err -o out/$(basename "$file").out --wrap='Rscript loader.R $(basename "$file")'" >> annotation_jobs.bch
done
chmod 775 annotation_jobs.bch

message0 "Downloading the action script..."
fetch_script loaderHadoop.R IMPC%20annotation%20pipeline
submit_limit_jobs annotation_jobs.bch ../../../../compressed_logs/annotation_job_id.txt
waitTillCommandFinish

message0 "Zipping logs..."
mv annotation_jobs.bch ../../../../compressed_logs
zip -rm ../../../../compressed_logs/annotation_logs.zip log/* err/* out/*
zip -rm splits.zip split_index_*
message0 "Job done."
