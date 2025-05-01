# IMPC Statistical Pipeline
This is the main R source code repository for IMPC statistical pipeline.

The IMPC statistical pipeline requires 2 steps to complete:
1. Pre-processing the data and run the statistical analysis.
2. Run the annotation pipeline.

```mermaid
%%{
    init: {
        "theme": "default",
        "themeVariables": {
            "fontSize": "15px"
        },
        "sequence": {
            "useMaxWidth": false
        }
    }
}%%
flowchart TB
    subgraph container[ ]
    style container fill:#ffffff
    direction TB
        subgraph stats_pipeline ["Step 1. Analysis ±2 weeks"]
            style stats_pipeline fill:#E6E6FA
            style main_ageing fill:#E6FAE6,stroke:#6FC56D
            style main_ageing_phase3 fill:#E6FAE6,stroke:#6FC56D
            direction LR
            subgraph phase1 ["Phase I. Preparing parquet files ±36 min"]
                direction TB
                inputStatsPipeline[StatsPipeline]-->|DRversion=20.2| step1[far:fa-file Step1MakePar2RdataJobs.R]
                step1 --> |Generate file with a list of jobs| step2_parquet2rdata{{jobs_step2_Parquet2Rdata.bch}}
                step2_parquet2rdata --> step2[far:fa-file Step2Parquet2Rdata.R] 
                step2_parquet2rdata --> |Run all jobs in .bch and \nwait until it's finished| step3[far:fa-file Step3MergeRdataFilesJobs.R]
                step2 --> step3
                step3 --> |Generate file with a list of jobs| step4_merge_rdatas{{jobs_step4_MergeRdatas.bch}}
                step4_merge_rdatas --> step4[far:fa-file Step4MergingRdataFiles.R] 
                step4_merge_rdatas --> |Run all jobs in .bch and \nwait until it's finished| compress_cleaning[Compress log files and clean up]
                step4 --> compress_cleaning
                compress_cleaning --> |zip -rm| parquet_to_rdata_jobs{{far:fa-folder Parquet2RdataJobs.zip}}
                compress_cleaning --> |zip -rm| parquet_to_rdata_logs{{far:fa-folder Parquet2RdataLogs.zip}}
                compress_cleaning --> |rm -rf| procedure_scatter_data{{far:fa-folder ProcedureScatterRdata}}
            end
            subgraph phase2 ["Phase II. Reprocessing the data ±5 days 14 hours"]
                direction TB
                job_creator[jobCreator from\nsideFunctions.R] --> |Generate file with jobs| data_generation_job_list{{DataGenerationJobList.bch}}
                data_generation_job_list --> input_data_generator[far:fa-file InputDataGenerator.R]
                data_generation_job_list --> |Run all jobs in .bch and \nwait until it's finished| compress_logs[Compress logs]
                input_data_generator --> generate_data[GenerateData from\nInputDataGenerator.R] 
                generate_data --> |GenerateData run\nmainAgeing function| main_ageing[mainAgeing from\nDRrequiredAgeing]
                main_ageing --> |BatchProducer = TRUE| compress_logs
                compress_logs --> remove_logs[Remove logs]
            end
            subgraph phase3 ["Phase III. Initialising the statistical analysis... ±6 days 22 hours"]
                direction TB
                update_impress[updateImpress from\nsideFunctions.R] --> windowing_pipeline{Is\nwindowingPipeline\nTrue?}
                windowing_pipeline --> |"True — default"| window_true[Copy function_windowed.R\n and rename to function.R]
                windowing_pipeline --> |Else| window_else[Copy function.R]
                window_true --> replace_word[ReplaceWordInFile from\nsideFunctions.R]
                window_else --> replace_word
                replace_word --> |ReplaceWordInFile use function.R| main_ageing_phase3[mainAgeing from\nDRrequiredAgeing]
                main_ageing_phase3 --> |BatchProducer = FALSE\nWait until completion| package_backup[packageBackup from\nsideFunctions.R]
            end
        end
        subgraph further_steps[ ]
            direction LR
            annotation["Step 2.Annotation\nand transfer pipeline\n±1 Day"] --> report["Step 3. Report\ngenerating pipeline\n±½ day"]
            report --> risky["Step 4. Extraction\nof risky genes pipeline\n±30 minutes"]
        end
        input[/ETL Parquet Files\] --> stats_pipeline --> further_steps
        mp_chooser[/mp_chooser\] --> stats_pipeline
        phase1 --> phase2
        phase2 --> phase3
    end

    classDef title font-size:30px
    class stats_pipeline title
```
# How to Run IMPC Statistical and Annotation Pipeline
These instructions are tailored for Release 21.0. To know more about input files for statistical pipeline refer to the [Observations Output Schema](https://github.com/mpi2/impc-etl/wiki/Observations-Output-Schema). In the current dataset, some fields that should be arrays are presented as comma-separated lists.

## Step 1. Data Preprocessing and Analysis
### 0. Switch to the mi_stats virtual user, start screen and activate working environment called R2D2
```console
become mi_stats
screen -S stats-pipeline
conda deactivate
conda activate R2D2
```

### 1. Set necessary variables
```console
export VERSION="21.0"
export REMOTE="mpi2"
export BRANCH="master"
export KOMP_PATH="<absolute_path_to_directory>"
```

### 2. Download script `orchestration.sh` that run both statistical and annotation pipeline on SLURM and add execute permission to a file
```console
cd ${KOMP_PATH}/impc_statistical_pipeline/IMPC_DRs/orchestration_scripts
wget https://raw.githubusercontent.com/${REMOTE}/impc_stats_pipeline/${BRANCH}/orchestration/orchestration.sh -O ${VERSION}_orchestration.sh
chmod +x ${VERSION}_orchestration.sh
```

### 3. Execute `orchestration.sh` script
```console
bash ${VERSION}_orchestration.sh ${VERSION} ${REMOTE} ${BRANCH} ${KOMP_PATH} ${KOMP_PATH}/data-releases/latest-input/dr${VERSION}/output/flatten_observations_parquet/ ${KOMP_PATH}/data-releases/latest-input/dr${VERSION}/output/mp_chooser_json/
```
- To leave screen, press combination `Ctrl + A + D`. Save screen session name, for example `3773511.stats-pipeline`. You will need it to reattach to the screen.
- You can track progress in the `${KOMP_PATH}/impc_statistical_pipeline/IMPC_DRs/stats_pipeline_logs/orchestration_${VERSION}.log` file or reattach to the screen with following command.
```console
screen -r 3773511.stats-pipeline
```

**Note:** Be cautious, the location of the input files may vary.
To execute `orchestration.sh` we need to pass six parameters:

  1. Version of the data release.
  2. Remote name.
  3. Branch name.
  4. Path to the initial directory.
  5. Path to the input parquet files.
  6. Path to the MP chooser file.

Seventh parameter is optional and by default is true. It indicates whether to use windowing or not. 

### Monitor progress using the following commands
- Use `squeue` to check list of running jobs.
- Use `jobinfo -v <job_id>` to check the job status.
- Review the log files:
```console
less ${KOMP_PATH}/impc_statistical_pipeline/IMPC_DRs/stats_pipeline_logs/stats_pipeline_${VERSION}.log
less ${KOMP_PATH}/impc_statistical_pipeline/IMPC_DRs/stats_pipeline_logs/stats_pipeline_${VERSION}.err
```
