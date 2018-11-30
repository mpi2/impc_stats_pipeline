library(nlme)
library(PhenStatAgeing)
library(PhenStat)
library(MASS)
# test = PhenList(
# 	dataset = read.csv(file = 'https://www.mousephenotype.org/data/exportraw?phenotyping_center=UC%20Davis&parameter_stable_id=IMPC_CSD_045_002&allele_accession_id=MGI:5548568&strain=MGI:2683688&pipeline_stable_id=UCD_001&&zygosity=homozygote&', na.strings = '-'),
# 	testGenotype = 'BL3208'
# )
# test2 = PhenListAgeing(x = test)


test2 = test
test2@datasetPL = read.csv('d:/ABR.test.csv')
#test2 = PhenListAgeing(x = test2)
cat('\014')
###
outputmodel = testDatasetAgeing(
	procedure = 'TYPICAL',
	phenListAgeing = (test2),
	depVariable = 'Value'   ,
	categorical = 0         ,
	withWeight = 0,
	LifeStage = 1
)

a = summary(outputmodel)
write(a$JSON,'d:hamed.json')
#write(stargazer(a$object$initial.model,a$object$final.model,type = 'html'),file = 'd:h.html')

#plot(a$object$final.model)
# phenListAgeing = (test)
# depVariable = 'Value'
# withWeight = TRUE


####################
# Mixed Model framework
file <- system.file("extdata", "test1.csv", package = "PhenStat")
test <-
	PhenStat:::PhenList(dataset = read.csv(file, na.strings = '-'),
											testGenotype = "Sparc/Sparc")
result <- PhenStat:::testDataset(test,
																 depVariable = "Lean.Mass")
tt = PhenListAgeing(
	x = test,
	DOB = 'Birth.Date',
	DOE = 'Batch',
	d.threshold = 95
)
a = testDatasetAgeing(
	phenListAgeing = tt,
	depVariable = 'Lean.Mass',
	LifeStage = 0,
	withWeight = 1
)
b=summary(a)
write(b$JSON,'d:hamed.json')
summary(result)
###################
# Fisher Exact Test framework
file =
	system.file("extdata", "test_categorical.csv", package = "PhenStat")
test2 =
	PhenStat:::PhenList(dataset = read.csv(file, na.strings = '-'),
											testGenotype = "Aff3/Aff3")

result2 = PhenStat:::testDataset(test2,
																	depVariable = "Thoracic.Processes",
																	method = "FE")



h = PhenListAgeing(x = test2, DOB = 'Birth.Date')
a=testDatasetAgeing(
	phenListAgeing = h,
	LifeStage = 0,
	depVariable = "Thoracic.Processes"
)

summary(result2)
b=summary(a)
write(b$JSON,'d:hamed.json')

