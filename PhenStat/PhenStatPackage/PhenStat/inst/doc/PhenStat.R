### R code from vignette source 'PhenStat.Rnw'

###################################################
### code chunk number 1: R.hide
###################################################
library(PhenStat)

dataset1 <- system.file("extdata", "test1.csv", package="PhenStat")

dataset2 <- system.file("extdata", "test1.txt", package="PhenStat")


###################################################
### code chunk number 2: R.hide
###################################################
# Default behaviour with messages
library(PhenStat)
dataset1 <- system.file("extdata", "test1.csv", package="PhenStat")
test <- PhenList(dataset=read.csv(dataset1,na.strings = '-'),
        testGenotype="Sparc/Sparc")

# Out-messages are switched off
test <- PhenList(dataset=read.csv(dataset1,na.strings = '-'),
        testGenotype="Sparc/Sparc",
        outputMessages=FALSE)


###################################################
### code chunk number 3: R.hide
###################################################
library(PhenStat)
dataset1 <- system.file("extdata", "test3.csv", package="PhenStat")

test <- PhenList(dataset=read.csv(dataset1,na.strings = '-'),
        dataset.clean=TRUE,
        dataset.values.female=1,
        dataset.values.male=2,
        testGenotype="Mysm1/+")


###################################################
### code chunk number 4: R.hide
###################################################
library(PhenStat)
dataset1 <- system.file("extdata", "test3.csv", package="PhenStat")

test <- PhenList(dataset=read.csv(dataset1,na.strings = '-'),
                dataset.clean=TRUE,
                dataset.values.female=1,
                dataset.values.male=2,
                testGenotype="Mysm1/+")

PhenStat:::getDataset(test)
test


###################################################
### code chunk number 5: R.hide
###################################################
library(PhenStat)
dataset2 <- system.file("extdata", "test2.csv", package="PhenStat")

test2 <- PhenList(dataset=read.csv(dataset2,na.strings = '-'),
                testGenotype="Arid4a/Arid4a",
                dataset.colname.weight="Weight.Value")

PhenStat:::testGenotype(test2)
PhenStat:::refGenotype(test2)


###################################################
### code chunk number 6: R.hide
###################################################
library(PhenStat)
dataset1 <- system.file("extdata", "test1.csv", package="PhenStat")

test <- PhenList(dataset=read.csv(dataset1,na.strings = '-'),
        testGenotype="Sparc/Sparc",
        outputMessages=FALSE)

# Default behaviour
result <- testDataset(test,
        depVariable="Bone.Area",
        equation="withoutWeight")

# Perform each step of the MM framework separatly
result <- testDataset(test,
        depVariable="Bone.Area",
        equation="withoutWeight",callAll=FALSE)

# Estimated model effects
linearRegressionResults <- PhenStat:::analysisResults(result)
linearRegressionResults$model.effect.batch

linearRegressionResults$model.effect.variance

linearRegressionResults$model.effect.weight

linearRegressionResults$model.effect.sex

linearRegressionResults$model.effect.interaction

# Change the effect values: interaction effect will stay in the model
result <- testDataset(test,
        depVariable="Bone.Area",
        equation="withoutWeight",
        keepList=c(TRUE,TRUE,FALSE,TRUE,TRUE),
        callAll=FALSE)

result <- PhenStat:::finalModel(result)

PhenStat:::summaryOutput(result)


###################################################
### code chunk number 7: R.hide
###################################################
PhenStat:::testFinalModel(result)

PhenStat:::classificationTag(result)


###################################################
### code chunk number 8: R.hide
###################################################
file <- system.file("extdata", "test7_TFE.csv", package="PhenStat")
test <- PhenList(dataset=read.csv(file,na.strings = '-'),
        testGenotype="het",
        refGenotype = "WT",
        dataset.colname.sex="sex",
        dataset.colname.genotype="Genotype",
        dataset.values.female="f",
        dataset.values.male= "m",
        dataset.colname.weight="body.weight",
        dataset.colname.batch="Date_of_procedure_start")

# TFDataset function creates cleaned dataset - concurrent controls dataset
test_TF <- PhenStat:::TFDataset(test,depVariable="Cholesterol")

# TF method is called
result  <- testDataset(test_TF,
        depVariable="Cholesterol",
        method="TF")
PhenStat:::summaryOutput(result)


###################################################
### code chunk number 9: R.hide
###################################################
library(PhenStat)
file <- system.file("extdata", "test1.csv", package="PhenStat")
test <- PhenList(dataset=read.csv(file,na.strings = '-'),
        testGenotype="Sparc/Sparc")

# RR method is called
result <- testDataset(test,
        depVariable="Lean.Mass",
        method="RR")
PhenStat:::summaryOutput(result)


###################################################
### code chunk number 10: R.hide
###################################################
library(PhenStat)
dataset_cat <- system.file("extdata", "test_categorical.csv", package="PhenStat")

test_cat <- PhenList(read.csv(dataset_cat,na.strings = '-'),testGenotype="Aff3/Aff3")

result_cat <- testDataset(test_cat,
        depVariable="Thoracic.Processes",
        method="FE")

PhenStat:::getVariable(result_cat)

PhenStat:::method(result_cat)

PhenStat:::summaryOutput(result_cat)


###################################################
### code chunk number 11: R.hide
###################################################
library(PhenStat)
dataset1 <- system.file("extdata", "test1.csv", package="PhenStat")

# MM framework
test <- PhenList(dataset=read.csv(dataset1,na.strings = '-'),
        testGenotype="Sparc/Sparc",outputMessages=FALSE)

result <- testDataset(test,
        depVariable="Lean.Mass",
        outputMessages=FALSE)

PhenStat:::summaryOutput(result)


###################################################
### code chunk number 12: R.hide
###################################################
library(PhenStat)
dataset_cat <- system.file("extdata", "test_categorical.csv", package="PhenStat")

test2 <- PhenList(dataset=read.csv(dataset_cat,na.strings = '-'),
        testGenotype="Aff3/Aff3",outputMessages=FALSE)

result2 <- testDataset(test2,
        depVariable="Thoracic.Processes",
        method="FE",outputMessages=FALSE)

PhenStat:::summaryOutput(result2)


###################################################
### code chunk number 13: R.hide
###################################################
library(PhenStat)
dataset_cat <- system.file("extdata", "test_categorical.csv", package="PhenStat")

test_cat <- PhenList(dataset=read.csv(dataset_cat,na.strings = '-'),
        testGenotype="Aff3/Aff3",outputMessages=FALSE)

result_cat <- testDataset(test_cat,
        depVariable="Thoracic.Processes",
        method="FE",outputMessages=FALSE)

PhenStat:::vectorOutput(result_cat)


###################################################
### code chunk number 14: R.hide
###################################################
library(PhenStat)
dataset_cat <- system.file("extdata", "test_categorical.csv", package="PhenStat")

test_cat <- PhenList(dataset=read.csv(dataset_cat,na.strings = '-'),
        testGenotype="Aff3/Aff3",outputMessages=FALSE)

result_cat <- testDataset(test_cat,
        depVariable="Thoracic.Processes",
        method="FE",outputMessages=FALSE)

#vectorOutputMatrices(result_cat)


###################################################
### code chunk number 15: R.hide
###################################################
library(PhenStat)
dataset_cat <- system.file("extdata", "test_categorical.csv", package="PhenStat")

test_cat <- PhenList(dataset=read.csv(dataset_cat,na.strings = '-'),
        testGenotype="Aff3/Aff3",outputMessages=FALSE)

result_cat <- testDataset(test_cat,
        depVariable="Thoracic.Processes",
        method="FE",outputMessages=FALSE)

PhenStat:::categoricalBarplot(result_cat)


###################################################
### code chunk number 16: R.hide
###################################################
library(PhenStat)
dataset1 <- system.file("extdata", "test1.csv", package="PhenStat")

# MM framework
test <- PhenList(dataset=read.csv(dataset1,na.strings = '-'),
        testGenotype="Sparc/Sparc",outputMessages=FALSE)

result <- testDataset(test,
        depVariable="Lean.Mass",
        outputMessages=FALSE)

PhenStat:::boxplotSexGenotype(test,
        depVariable="Lean.Mass",
        graphingName="Lean Mass")

PhenStat:::scatterplotSexGenotypeBatch(test,
        depVariable="Lean.Mass",
        graphingName="Lean Mass")

PhenStat:::scatterplotGenotypeWeight(test,
        depVariable="Bone.Mineral.Content",
        graphingName="BMC")


###################################################
### code chunk number 17: R.hide
###################################################
library(PhenStat)
dataset1 <- system.file("extdata", "test1.csv", package="PhenStat")

# MM framework
test <- PhenList(dataset=read.csv(dataset1,na.strings = '-'),
        testGenotype="Sparc/Sparc",outputMessages=FALSE)

result <- testDataset(test,
        depVariable="Lean.Mass",
        outputMessages=FALSE)

PhenStat:::qqplotGenotype(result)

PhenStat:::qqplotRandomEffects(result)

PhenStat:::qqplotRotatedResiduals(result)

PhenStat:::plotResidualPredicted(result)

PhenStat:::boxplotResidualBatch(result)


