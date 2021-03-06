## Copyright © 2012-2014 EMBL - European Bioinformatics Institute
##
## Licensed under the Apache License, Version 2.0 (the "License");
## you may not use this file except in compliance with the License.
## You may obtain a copy of the License at
##
##     http://www.apache.org/licenses/LICENSE-2.0
##
## Unless required by applicable law or agreed to in writing, software
## distributed under the License is distributed on an "AS IS" BASIS,
## WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
## See the License for the specific language governing permissions and
## limitations under the License.
##------------------------------------------------------------------------------
## PhenList.R contains PhenList and checkDataset functions to construct
## new PhenList object from components and to check dataset integrity
##------------------------------------------------------------------------------
## Construct PhenList object from components with data quality checks
PhenList <-
  function(dataset,
           testGenotype,
           refGenotype = '+/+',
           hemiGenotype = NULL,
           outputMessages = TRUE,
           dataset.clean = TRUE,
           dataset.colname.batch = NULL,
           dataset.colname.genotype = NULL,
           dataset.colname.sex = NULL,
           dataset.colname.weight = NULL,
           dataset.values.missingValue = " ",
           dataset.values.male = NULL,
           dataset.values.female = NULL)
  {
    if (class(dataset) == "data.frame") {
      dataset <- dataset[, order(names(dataset))]

      ## Rename columns if needed
      if (dataset.clean) {
        if (!is.null(dataset.colname.batch) &&
            (dataset.colname.batch != 'Batch')) {
          # check here for the existing column named 'Batch'
          if ('Batch' %in% colnames(dataset)) {
            colnames(dataset)[colnames(dataset) == 'Batch'] <-
              'Original.Batch'
          }
          colnames(dataset)[colnames(dataset) == dataset.colname.batch] <-
            'Batch'
        }
        else {
          if ('Assay.Date' %in% colnames(dataset)) {
            if ('Batch' %in% colnames(dataset)) {
              colnames(dataset)[colnames(dataset) == 'Batch'] <-
                'Original.Batch'
            }
            colnames(dataset)[colnames(dataset) == 'Assay.Date'] <-
              'Batch'
            if (outputMessages)
              message(
                paste(
                  "Warning:\nDataset's column 'Assay.Date' has been ",
                  "renamed to 'Batch' and will be used for the batch effect modelling.\n",
                  sep = ""
                )
              )
          }
          else
            if ('AssayDate' %in% colnames(dataset)) {
              # check for the existing column named 'Batch'
              if ('Batch' %in% colnames(dataset)) {
                colnames(dataset)[colnames(dataset) == 'Batch'] <-
                  'Original.Batch'
              }
              colnames(dataset)[colnames(dataset) == 'AssayDate'] <-
                'Batch'
              if (outputMessages)
                message(
                  paste(
                    "Warning:\nDataset's column 'AssayDate' has been ",
                    "renamed to 'Batch' and will be used for the batch effect modelling.\n",
                    sep = ""
                  )
                )
            }
          else
            if (length(colnames(dataset)[grep("batch",
                                              tolower(colnames(dataset)))]) >
                0 && outputMessages) {
              batch_potential_columns <- paste(colnames(dataset)[grep
                                                                 ("batch", tolower(colnames(dataset)))], collapse =
                                                 "', '")

              message(
                paste(
                  "Warning:\nDataset contains columns that might ",
                  "be used for Batch effect modeling, for instance '",
                  batch_potential_columns,
                  "'.\n",
                  sep = ""
                )
              )
            }

        }
        if (!is.null(dataset.colname.genotype) &&
            (dataset.colname.genotype != 'Genotype')) {
          # check for the existing column named 'Genotype'
          if ('Genotype' %in% colnames(dataset)) {
            colnames(dataset)[colnames(dataset) == 'Genotype'] <-
              'Original.Genotype'
          }
          colnames(dataset)[colnames(dataset) == dataset.colname.genotype] <-
            'Genotype'
        }
        if (!is.null(dataset.colname.sex) &&
            (dataset.colname.sex != 'Sex') &&
            (dataset.colname.sex %in% colnames(dataset))) {
          # check for the existing column named 'Sex'
          if ('Sex' %in% colnames(dataset)) {
            colnames(dataset)[colnames(dataset) == 'Sex'] <- 'Original.Sex'
          }
          colnames(dataset)[colnames(dataset) == dataset.colname.sex] <-
            'Sex'
        }
        else {
          if ('Gender' %in% colnames(dataset)) {
            # check for the existing column named 'Sex'
            if ('Sex' %in% colnames(dataset)) {
              colnames(dataset)[colnames(dataset) == 'Sex'] <- 'Original.Sex'
            }
            colnames(dataset)[colnames(dataset) == 'Gender'] <-
              'Sex'
            if (outputMessages)
              message(
                paste(
                  "Warning:\nDataset's column 'Gender' has been ",
                  "renamed to 'Sex'.\n",
                  sep = ""
                )
              )
          }
        }

        if (!is.null(dataset.colname.weight) &&
            (dataset.colname.weight != 'Weight')) {
          # check for the existing column named 'Weight'
          if ('Weight' %in% colnames(dataset)) {
            colnames(dataset)[colnames(dataset) == 'Weght'] <-
              'Original.Weight'
          }
          colnames(dataset)[colnames(dataset) == dataset.colname.weight] <-
            'Weight'
        }
        ## Replace missing values specified in the user format with NA
        if (!is.null(dataset.values.missingValue))
          dataset[dataset == dataset.values.missingValue] <- NA
        ## Replace empty strings with NA
        dataset[dataset == ""] <- NA


        # make Weight column numeric if possible (if there are no strings)
        if ('Weight' %in% colnames(dataset)) {
          #for (i in 1:length(colnames(dataset))){
          columnName <- "Weight"
          checkNA_transformed <-
            sum(is.na(suppressWarnings(as.numeric(
              as.character(dataset[, c(columnName)])
            ))))
          checkNA_initial <- sum(is.na(dataset[, c(columnName)]))
          if (checkNA_transformed == checkNA_initial) {
            dataset[, c(columnName)] <-
              as.numeric(as.character(dataset[, c(columnName)]))
          }

          #}
        }



        #        if ('Weight' %in% colnames(dataset)){
        #            if (is.numeric(dataset$Weight) && (!all(sapply(dataset$Weight,is.na)))) {
        #                dataset$Weight<-as.numeric(dataset$Weight)
        #            }
        #            else {
        #                colnames(dataset)[colnames(dataset) == 'Weight'] <-'Weight_labels'
        #                if (outputMessages)
        #                message(paste("Warning:\nWeight column values are not numeric or NA. ",
        #                              "In order to avoid erroneous execution of statistical ",
        #                              "functions column is renamed to 'Weight_labels'.\n",sep=""))
        #
        #            }
        #        }

        ## Renew levels
        if ('Sex' %in% colnames(dataset)) {
          dataset$Sex <- factor(dataset$Sex)
          # # Replace values for sexes with 'Male','Female'
          levels(dataset$Sex)[levels(dataset$Sex) == "female"] <-
            "Female"
          levels(dataset$Sex)[levels(dataset$Sex) == "male"] <-
            "Male"

          if (!is.null(dataset.values.female))
            levels(dataset$Sex)[levels(dataset$Sex) == dataset.values.female] <-
            "Female"

          if (!is.null(dataset.values.male))
            levels(dataset$Sex)[levels(dataset$Sex) == dataset.values.male] <-
            "Male"
        }
        if ('Genotype' %in% colnames(dataset))
          dataset$Genotype <- factor(dataset$Genotype)
        if ('Batch' %in% colnames(dataset))
          dataset$Batch <- factor(dataset$Batch)




        ## Hemi to test genotype replacement
        if (!is.null(hemiGenotype)) {
          if (length(rownames(dataset[dataset$Genotype == hemiGenotype,])) > 0) {
            levels(dataset$Genotype)[levels(dataset$Genotype) == hemiGenotype] <-
              testGenotype
            if (outputMessages)
              message(
                paste(
                  "Warning:\nHemizygotes '",
                  hemiGenotype,
                  "' have been relabelled to test genotype '",
                  testGenotype,
                  "'.\nIf you don't want this behaviour then don't define ",
                  "'hemiGenotype' argument.\n",
                  sep = ""
                )
              )
          }
        }

        ## Clean genotypes
        if (length(setdiff(rownames(dataset),
                           rownames(dataset[dataset$Genotype %in% c(testGenotype, refGenotype),]))) >
            0) {
          dataset <-
            dataset[dataset$Genotype %in% c(testGenotype, refGenotype),]
          if (outputMessages)
            message(
              paste(
                "Warning:\nDataset has been cleaned by ",
                "filtering out records with genotype value other than test ",
                "genotype '",
                testGenotype,
                "' or reference genotype '",
                refGenotype,
                "'.\n",
                sep = ""
              )
            )

        }
      }


      ## Clean the empty records  NB - after renaming/cleaning !
      if ('Sex' %in% colnames(dataset)) {
        dataset <- dataset[dataset$Sex != "",]
        dataset <- dataset[!is.na(dataset$Sex),]
      }
      if ('Genotype' %in% colnames(dataset)) {
        dataset <- dataset[dataset$Genotype != "",]
        dataset <- dataset[!is.na(dataset$Genotype),]
      }
      if ('Batch' %in% colnames(dataset)) {
        dataset <- dataset[dataset$Batch != "",]
        dataset <- dataset[!is.na(dataset$Batch),]
      }


      ## Renew levels
      if ('Sex' %in% colnames(dataset))
        dataset$Sex <- factor(dataset$Sex)
      if ('Genotype' %in% colnames(dataset))
        dataset$Genotype <- factor(dataset$Genotype)
      if ('Batch' %in% colnames(dataset))
        dataset$Batch <- factor(dataset$Batch)

      ## CHECKS
      dataset_unfiltered <-
        dataset  # RR special case - unfiltered dataset
      dataset <- checkDataset(dataset,
                              testGenotype,
                              refGenotype,
                              outputMessages,
                              dataset.clean)

      checkWeight <- columnChecks(dataset, "Weight", 2)

      if (!checkWeight[1]) {
        if (outputMessages)
          message("Warning:\nWeight column is not present in the database.\n")
      }
      else {
        if (!checkWeight[2]) {
          if (outputMessages)
            message(
              paste(
                "Warning:\nWeight column values are not numeric.\n",
                "In order to avoid erroneous execution of statistical ",
                "functions column is renamed to 'Weight_labels'.\n",
                sep = ""
              )
            )
          colnames(dataset)[colnames(dataset) == 'Weight'] <-
            'Weight_labels'
        }
        if (!checkWeight[3]) {
          if (outputMessages)
            message(
              paste(
                "Warning:\nWeight column does not have enough data points ",
                "for genotype/sex combinations.\n",
                "In order to avoid erroneous execution of statistical ",
                "functions column is renamed to 'Weight_labels'.\n",
                sep = ""
              )
            )
          colnames(dataset)[colnames(dataset) == 'Weight'] <-
            'Weight_labels'
        }
      }

      Genotype_levels <- levels(dataset$Genotype)
      Sex_levels <- levels(dataset$Sex)

      new(
        "PhenList",
        datasetPL = as.data.frame(dataset),
        refGenotype = refGenotype,
        testGenotype = testGenotype,
        hemiGenotype = ifelse(is.null(hemiGenotype), character(0), hemiGenotype),
        dataset.colname.batch = ifelse(
          is.null(dataset.colname.batch),
          character(0),
          dataset.colname.batch
        ),
        dataset.colname.genotype = ifelse(
          is.null(dataset.colname.genotype),
          character(0),
          dataset.colname.genotype
        ),
        dataset.colname.sex = ifelse(
          is.null(dataset.colname.sex),
          character(0),
          dataset.colname.sex
        ),
        dataset.colname.weight = ifelse(
          is.null(dataset.colname.weight),
          character(0),
          dataset.colname.weight
        ),
        dataset.values.missingValue = ifelse(
          is.null(dataset.values.missingValue),
          character(0),
          dataset.values.missingValue
        ),
        dataset.values.male = ifelse(
          is.null(dataset.values.male),
          character(0),
          as.character(dataset.values.male)
        ),
        dataset.values.female = ifelse(
          is.null(dataset.values.female),
          character(0),
          as.character(dataset.values.female)
        ),
        dataset.clean = dataset.clean,
        datasetUNF = as.data.frame(dataset_unfiltered)
      )
      # RR special case - unfiltered dataset

    }
    else {
      message <-
        "Error: PhenList function's first argument should be data frame.\n"
      if (outputMessages) {
        message(message)
        opt <- options(show.error.messages = FALSE)
        on.exit(options(opt))
        stop()
      }
      else {
        stop(message)
      }
    }
  }
#-------------------------------------------------------------------------------
## Check dataset for the minimum required info and additional cleaning steps
checkDataset <- function(dataset,
                         testGenotype,
                         refGenotype = "+/+",
                         outputMessages = TRUE,
                         dataset.clean = TRUE)
{
  message <- ""
  message_dp <- ""
  pass <- TRUE

  nvar <- ncol(dataset)
  ntags <- nrow(dataset)

  ## Column names should be given
  if (nvar > 0 && is.null(colnames(dataset))) {
    pass <- FALSE
    message <- paste(message,
                     "\nCheck failed:\nDataset with no column names.\n",
                     sep = "")
  }

  ## Check for mandatory columns: Genotype and Sex
  if (!('Genotype' %in% colnames(dataset))) {
    pass <- FALSE
    message <- paste(message,
                     "\nCheck failed:\nDataset's 'Genotype' column is missed.\n",
                     sep = "")
  }


  if (!('Sex' %in% colnames(dataset))) {
    pass <- FALSE
    message <- paste(message,
                     "\nCheck failed:\nDataset's 'Sex' column is missed.\n",
                     sep = "")
  }

  ## Check for other columns: Weight and Batch
  if (!('Weight' %in% colnames(dataset)) && outputMessages) {
    message(
      paste(
        "Warning:\nDataset's 'Weight' column is missed.\n",
        "You can define 'dataset.colname.weight' argument to specify column ",
        "for the weight effect modeling.\n",
        sep = ""
      )
    )
  }

  if (!('Batch' %in% colnames(dataset))) {
    message(
      paste(
        "Warning:\nDataset's 'Batch' column is missed.\n",
        "You can define 'dataset.colname.batch' argument to specify column ",
        "for the batch effect modeling.\n",
        sep = ""
      )
    )
  }

  if (('Sex' %in% colnames(dataset)) &&
      ('Genotype' %in% colnames(dataset))) {
    Genotype_levels <- levels(dataset$Genotype)
    Sex_levels <- levels(dataset$Sex)

    ## Genotype/Sex combinations with less than two data points
    combinations_list <-
      "" # String with Genotype/Sex (count) in initial dataset

    for (i in 1:length(Genotype_levels)) {
      GenotypeSubset <-
        subset(dataset, dataset$Genotype == Genotype_levels[i])
      for (j in 1:length(Sex_levels)) {
        nr <-
          sum(is.finite(GenotypeSubset[GenotypeSubset$Sex == Sex_levels[j],][, "Sex"]))
        combinations_list <-
          paste(
            combinations_list,
            paste(
              "'",
              Genotype_levels[i],
              "'/'",
              Sex_levels[j],
              "' (",
              nr,
              "), ",
              sep = ""
            ),
            sep = ""
          )
      }

    }

    ## String with Genotype/Sex (count), where Sex would be filtered out
    filtered_list_combinations <- ""
    dataset_filtered <- dataset
    for (i in 1:length(Genotype_levels)) {
      GenotypeSubset <-
        subset(dataset, dataset$Genotype == Genotype_levels[i])
      for (j in 1:length(Sex_levels)) {
        nr <-
          sum(is.finite(GenotypeSubset[GenotypeSubset$Sex == Sex_levels[j],][, "Sex"]))

        ## There are combinations with less than two data points
        if (nr < 2) {
          filtered_list_combinations <- paste(
            filtered_list_combinations,
            paste(
              "'",
              Genotype_levels[i],
              "'/'",
              Sex_levels[j],
              "' (",
              nr,
              "), ",
              sep = ""
            ),
            sep = ""
          )

          if (dataset.clean) {
            ## If you have data in one genotype for both sexes
            ## but not in the other then you have to revert to a
            ## one sex analysis.
            subset_to_filter <- subset(dataset_filtered,
                                       (dataset_filtered$Sex == Sex_levels[j]))

            dataset_filtered <-
              dataset_filtered[setdiff(rownames(dataset_filtered),
                                       rownames(subset_to_filter)),]

          }
          else{
            pass <- FALSE
          }
        }
      }

    }


    if (nchar(filtered_list_combinations) > 2) {
      combinations_list <- substr(combinations_list, 1,
                                  nchar(combinations_list) - 2)

      filtered_list_combinations <-
        substr(filtered_list_combinations,
               1,
               nchar(filtered_list_combinations) - 2)

      if (dataset.clean) {
        if (outputMessages) {
          message(
            paste(
              "Warning:\nSince dataset has to have at least ",
              "two data points for each genotype/sex combination and ",
              "there are not enough records for the combination(s): ",
              filtered_list_combinations,
              ", appropriate sex records ",
              "have been filtered out from the dataset.",
              sep = ""
            )
          )
          message(
            paste(
              "N.B. for RR and FE methods:\nBy default RR and FE will use filtered dataset. ",
              "However, you can use unfiltered dataset which includes all data points by setting 'testDataset' function's argument 'useUnfiltered' to TRUE .\n",
              sep = ""
            )
          )
        }
      }
      else
        message_dp <-
        paste(
          "\nCheck failed:\nDataset should have at ",
          "least two data points for each genotype/sex combination.\n",
          "At the moment there are no enough data points for the following ",
          "combination(s): ",
          filtered_list_combinations,
          ".\n",
          sep = ""
        )
    }

    dataset <- dataset_filtered

    ## Renew levels
    dataset$Sex <- factor(dataset$Sex)
    dataset$Genotype <- factor(dataset$Genotype)

    if ('Batch' %in% colnames(dataset))
      dataset$Batch <- factor(dataset$Batch)


    Genotype_levels <- levels(dataset$Genotype)
    Sex_levels <- levels(dataset$Sex)

    ## INFO about genotype and sex levels
    genotype_values <- paste(Genotype_levels, collapse = "', '")
    sex_values <- paste(Sex_levels, collapse = "', '")
    if (outputMessages) {
      message(
        paste(
          "Information:\nDataset's 'Genotype' ",
          "column has following values: '",
          genotype_values,
          "'\n",
          sep = ""
        )
      )
      if (nchar(filtered_list_combinations) > 2) {
        message(
          paste(
            "Information:\nFiltered dataset's 'Sex' ",
            "column has following value(s): '",
            sex_values,
            "'\n",
            sep = ""
          )
        )
      }
      else {
        message(
          paste(
            "Information:\nDataset's 'Sex' ",
            "column has following value(s): '",
            sex_values,
            "'\n",
            sep = ""
          )
        )
      }
    }

    ## Check of genotype and sex levels after cleaning
    if (length(Genotype_levels) != 2)  {
      pass <- FALSE
      message <-
        paste(
          message,
          "\nCheck failed:\nDataset's 'Genotype' ",
          "column has to have two values.\nYou can define 'testGenotype' and ",
          "'refGenotype' arguments to automatically filter out records with ",
          "genotype values other than specified.\nAlternatively you can define ",
          "'hemiGenotype' and 'testGenotype' arguments to relabel hemizygotes ",
          "to homozygotes.\n",
          sep = ""
        )
    }


    if (!(length(Sex_levels) %in% c(1, 2))) {
      pass <- FALSE
      message <-
        paste(
          message,
          "\nCheck failed:\nDataset's 'Sex' ",
          "column has to have one or two values and currently the data has ",
          "more than two.\n",
          sep = ""
        )

    }

    ## Check for sex levels - we want to have 'Female' and/or 'Male' only
    wrong_sex_levels <- setdiff(Sex_levels, c("Female", "Male"))
    wrong_sex_values <- paste(wrong_sex_levels, collapse = "', '")

    if (!length(wrong_sex_levels) == 0) {
      pass <- FALSE
      if (length(Sex_levels) <= 2)
        message <- paste(
          message,
          "\nCheck failed:\nDataset's ",
          "'Sex' column has '",
          sex_values,
          "' values instead of ",
          "'Female' and/or 'Male' values. You can define ",
          "'dataset.values.male' and 'dataset.values.female' ",
          "arguments to replace those values automatically.\n",
          sep = ""
        )
      else
        message <- paste(
          message,
          "\nCheck failed:\nDataset's ",
          "'Sex' column has '",
          sex_values,
          "' values instead of ",
          "'Female' and/or 'Male' values only. ",
          "Please delete records with sex(es) '",
          wrong_sex_values,
          "' from the dataset.\n",
          sep = ""
        )
    }


    message <- paste(message, message_dp, sep = "")

    ## Check for reference genotype records
    #if (sum(grepl(refGenotype, Genotype_levels, fixed=TRUE))==1)
    if (refGenotype %in% Genotype_levels)
      dataset$Genotype = relevel(dataset$Genotype, ref = refGenotype)
    else {
      pass <- FALSE
      message <-
        paste(
          message,
          paste(
            "\nCheck failed:\nDataset with not ",
            "enough records for statistical analysis with reference genotype '",
            refGenotype,
            "'.\n",
            sep = ""
          )
        )
    }

    ## Check for test genotype records
    #if (!(sum(grepl(testGenotype, Genotype_levels, fixed=TRUE))==1)){
    if (!(testGenotype %in% Genotype_levels)) {
      pass <- FALSE
      message <-
        paste(
          message,
          "\nCheck failed:\nDataset ",
          "with not enough records for statistical analysis with test ",
          "genotype '",
          testGenotype,
          "'.\n",
          sep = ""
        )
    }
  }

  if (!pass) {
    if (outputMessages) {
      message(paste("********* Errors start *********\n", message, sep = ""))
      message("********* Errors end ***********")
      opt <- options(show.error.messages = FALSE)
      on.exit(options(opt))
      stop()
    }
    else {
      stop(message)
    }

  }

  return(dataset)
}
##------------------------------------------------------------------------------
