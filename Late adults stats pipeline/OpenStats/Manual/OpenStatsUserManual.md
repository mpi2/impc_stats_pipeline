**OpenStats**
================
Hamed Haseli Mashhadi (<hamedhm@ebi.ac.uk>),
03 January 2020

## Introduction

Coping with a large volumes of categorical and continuous data in the
high-throughput phenotyping pipelines requires a robust automated
statistical pipeline that minimal alleviates manual intervention is
required.

**OpenStats** provides a variety of the statistical work-flows for the
identification of abnormal phenotypes with an emphasize on
high-throughput. The embeded set of checks and on-fly fixes assures an
accurate and robust analysis for high-throughput data and a
comprehensive output allows transfering the results to the downstream
steges.

## Package architecture

A schematic illustration of three layer internal workflows of OpenStats
is given below. More details can be found in the appropriate sections of
the User’s Guide.

![alt text here](workflow.png)

Three first lawer allows a dataset as the input and a proper input
model. In the second layer, the input data as well as the model are
checked for feasibility and correctness and reports any unusuality in
the data or the input model. The successful output of this stage is a
OpenStatsList object that is the input of the other downstream layer.

The statistical analysis layer encompasses three analysis framework
precisely, Linear Mixed Model and Reference range plus for the
continuous and the Fisher’s exact test for the categorical data.The are
extra checks and on-fly fixes in this stage to assure a successful,
correct and comprehensive analysis will be applied to the data. More
details are provided in the corresponding section.

The final layer provides feed for the downstream processes and consists
of summaries, plots list an JSON objects.

## Input data and model

OpenStats allows a data frame for the input. The input model must be in
the standard form of `y~x+z+t+...` no function is allowed in the
formula. For instance, `log(y)~x+z+t+..` or `y~x^2+z+t+...` are invalid.

## Data preparation with `OpenStatsList` function

`OpenStatsList` function performs data processing and creates a
OpenStats object. This function allows a data frame as the input so that
the rows and columns represent the samples and features respectively. In
addition to dependent variable column (the variable of interest)
mandatory column “Genotype” or the corresponded proxy is required. Note
that any other name for “Genotype” is accepted however the OpenStatsList
function rename that to “Genotype”. To align with **PhenStat** the
function allows optional *Sex*, *Batch*, *Weight* (bodyweight) and an
extra *LifeStage* in the input parameters.

###### The main tasks performed by the OpenStats package’s function OpenStatsList are:

  - Preparing a workind dataset for the downstream operations
  - terminology normalisation
  - showing a general view of the dataset
  - filtering out undesirable records with respect to the input model
  - checking whether the dataset can be used for the statistical
    analysis

Allchecks are accompanied with the informative messages, warnings and
errors. One example of the `OpenStatsList` is shown below:

``` r
  ####################################################################
  df = read.csv(system.file("extdata", "test_continuous.csv", package = "OpenStats"))
  ####################################################################
  # OpenStatsList object
  ####################################################################
  OpenStatsList   = OpenStatsList(
    dataset       = df,
    testGenotype  = 'experimental',
    refGenotype   = 'control',
    dataset.colname.batch = 'date_of_experiment',
    dataset.colname.genotype = 'biological_sample_group',
    dataset.colname.sex = 'sex',
    dataset.colname.weight = 'weight'
  )
```

    ## 2020-01-03 10:15:46. Input data of the dimensions, rows = 410, columns = 75

    ## 2020-01-03 10:15:46. Checking the input data in progress ...

    ## 2020-01-03 10:15:46. Checking the specified missing values [x2] (` `, ``) ...

    ## 2020-01-03 10:15:46.     1/2. Checking (` `) ...

    ## 2020-01-03 10:15:47.     2/2. Checking (``) ...

    ## 2020-01-03 10:15:47. Checking whether variable `biological_sample_group` exists in the data ... 
    ## 2020-01-03 10:15:47.     Result = TRUE

    ## 2020-01-03 10:15:47.      Levels (Total levels = 2, missings = 0%): 
    ## 2020-01-03 10:15:47.       control
    ## 2020-01-03 10:15:47.       experimental

    ## 2020-01-03 10:15:47. Checking whether variable `sex` exists in the data ... 
    ## 2020-01-03 10:15:47.     Result = TRUE

    ## 2020-01-03 10:15:47.      Levels (Total levels = 2, missings = 0%): 
    ## 2020-01-03 10:15:47.       female
    ## 2020-01-03 10:15:47.       male

    ## 2020-01-03 10:15:47. Checking whether variable `date_of_experiment` exists in the data ... 
    ## 2020-01-03 10:15:47.     Result = TRUE

    ## 2020-01-03 10:15:47.      Levels (Total levels = 43, missings = 0%): 
    ## 2020-01-03 10:15:47.       2012-07-23T00:00:00Z
    ## 2020-01-03 10:15:47.       2012-07-30T00:00:00Z
    ## 2020-01-03 10:15:47.       2012-08-06T00:00:00Z
    ## 2020-01-03 10:15:47.       2012-08-13T00:00:00Z
    ## 2020-01-03 10:15:47.       2012-08-20T00:00:00Z
    ## 2020-01-03 10:15:47.       2012-11-26T00:00:00Z
    ## 2020-01-03 10:15:47.       2012-12-24T00:00:00Z
    ## 2020-01-03 10:15:47.       2013-01-02T00:00:00Z
    ## 2020-01-03 10:15:47.       2013-01-15T00:00:00Z
    ## 2020-01-03 10:15:47.       2013-01-21T00:00:00Z
    ## 2020-01-03 10:15:47.       2013-01-28 ...

    ## 2020-01-03 10:15:47. Checking whether variable `LifeStage` exists in the data ... 
    ## 2020-01-03 10:15:47.     Result = FALSE

    ## 2020-01-03 10:15:47. Checking whether variable `weight` exists in the data ... 
    ## 2020-01-03 10:15:47.     Result = TRUE

    ## 2020-01-03 10:15:47.      Summary:
    ## 2020-01-03 10:15:47.       mean      = 20.0036585365854
    ## 2020-01-03 10:15:47.       sd        = 2.63972182732584
    ## 2020-01-03 10:15:47.       Missings  = 0%

    ## 2020-01-03 10:15:47. Variable `biological_sample_group` renamed to `Genotype`

    ## 2020-01-03 10:15:47. Variable `sex` renamed to `Sex`

    ## 2020-01-03 10:15:47. Variable `date_of_experiment` renamed to `Batch`

    ## 2020-01-03 10:15:47. Variable `weight` renamed to `Weight`

    ## 2020-01-03 10:15:47. Total samples in Genotype:Sex interactions: 
    ## 2020-01-03 10:15:47.     Level(frequency): 
    ## 2020-01-03 10:15:47.      control.Female(201)
    ## 2020-01-03 10:15:47.      experimental.Female(0)
    ## 2020-01-03 10:15:47.      control.Male(201)
    ## 2020-01-03 10:15:47.      experimental.Male(8)

    ## 2020-01-03 10:15:47. No observations detected in Genotype:Sex interactions for:
    ## 2020-01-03 10:15:47.     experimental.Female

    ## 2020-01-03 10:15:47. Total `Weight` data points for Genotype:Sex interactions:
    ## 2020-01-03 10:15:47.      Level(frequency): 
    ## 2020-01-03 10:15:47.       control.Female(201),
    ## 2020-01-03 10:15:47.       experimental.Female(0),
    ## 2020-01-03 10:15:47.       control.Male(201),
    ## 2020-01-03 10:15:47.       experimental.Male(8)

    ## 2020-01-03 10:15:47. `Weight` column has (<2) data points for at least one level of Genotype:Sex interactions.
    ## 2020-01-03 10:15:47.      The `Weight` column renamed to `Weight_labels`

    ## 2020-01-03 10:15:47. Successfully performed checks in 0.1 second(s).

## Including Plots

You can also embed plots, for example:

    ## 12345678910
    ## 12345678910

Note that the `echo = FALSE` parameter was added to the code chunk to
prevent printing of the R code that generated the plot.
