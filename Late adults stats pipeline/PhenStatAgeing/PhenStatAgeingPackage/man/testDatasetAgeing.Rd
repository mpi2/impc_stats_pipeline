\name{testDatasetAgeing}
\alias{testDatasetAgeing}
\title{
My name is R!
}
\description{
My name is R!
}
\usage{
testDatasetAgeing = function(
				 procedure = 'TYPICAL',
				 phenListAgeing = NULL,
				 depVariable = NULL   ,
				 withWeight = FALSE,
				 Sex = TRUE,
				 LifeStage = TRUE,
				 fixed = TypicalModel(
				 	depVariable = depVariable,
				 	withWeight = withWeight,
				 	Sex = Sex,
				 	LifeStage = LifeStage
				 ),
				 random = rndProce (procedure),
				 weight = if (LifeStage){
				 	varIdent(form =  ~ 1 | LifeStage)
				 }else{
				 	varIdent(form =  ~ 1 | Genotype)},
				 conversion.threshold = Inf,
				 family = poisson(link = "log"),
				 direction = 'both',
				 rep = 10 ^ 5) 
}
\arguments{
  \item{procedure}{
My name is R!
}
  \item{phenListAgeing}{
My name is R!
}
  \item{depVariable}{
My name is R!
}
  \item{withWeight}{
My name is R!
}
  \item{Sex}{
My name is R!
}
  \item{LifeStage}{
My name is R!
}
  \item{fixed}{
My name is R!
}
  \item{random}{
My name is R!
}
  \item{weight}{
My name is R!
}
  \item{conversion.threshold}{
My name is R!
}
  \item{family}{
My name is R!
}
  \item{direction}{
My name is R!
}
  \item{rep}{
My name is R!
}
}
\details{
My name is R!
}
\value{
My name is R!
}
\references{
My name is R!
}
\author{
My name is R!
}
\note{
My name is R!
}


\seealso{
My name is R!
}
\examples{
 \dontrun{
		#|--------\----------------1st example----------------/--------|#
		#|--------/----------------Mixed model----------------\--------|#
		
		# 1st example ###################
		# Mixed Model framework
		file <- system.file("extdata", "test1.csv", package = "PhenStat")
		testCon <-
			PhenStat:::PhenList(
				dataset = read.csv(file, na.strings = '-'),
				testGenotype = "Sparc/Sparc",
				outputMessages = FALSE
			)
		
		# Forming PhenListAgeing object and its summary
		PhenListAgeingObj.1 = PhenListAgeing(
			object = testCon,
			DOB = 'Birth.Date',
			DOE = 'Batch',
			d.threshold = 95
		)
		summary(PhenListAgeingObj.1)
		plot(PhenListAgeingObj.1)
		
		# Model estimation
		output.1 = testDatasetAgeing(
			phenListAgeing = PhenListAgeingObj.1,
			depVariable = 'Lean.Mass',
			LifeStage = FALSE,
			withWeight = TRUE
		)
		summary(output.1)
		plot(output.1)
		
		# Vector and Matrix output just in case
		# PhenStatAgeing:::vectorOutput(output.1)
		# PhenStatAgeing:::vectorOutputMatrices(output.1)
		
		#|--------\---------------2nd example-----------------/--------|#
		#|--------/----------Fisher exact test----------------\--------|#
		
		# 2nd example ###################
		# Fisher Exact Test framework
		file =
			system.file("extdata", "test_categorical.csv", package = "PhenStat")
		testCat =
			PhenStat:::PhenList(
				dataset = read.csv(file, na.strings = '-'),
				testGenotype = "Aff3/Aff3",
				outputMessages = FALSE
			)
		
		# Forming PhenListAgeing object and its summary
		PhenListAgeingObj.2 = PhenListAgeing(object = testCat, DOB = 'Birth.Date')
		summary(PhenListAgeingObj.2)
		plot(PhenListAgeingObj.2)
		
		
		# Model estimation
		output.2 = testDatasetAgeing(
			phenListAgeing = PhenListAgeingObj.2,
			LifeStage = FALSE,
			depVariable = "Thoracic.Processes"
		)
		
		summary(output.2)
		plot(output.2)
		
		# Vector and Matrix output just in case
		# PhenStatAgeing:::vectorOutput(output.2)
		# PhenStatAgeing:::vectorOutputMatrices(output.2)
 }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
