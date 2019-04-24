\name{testDatasetAgeing}
\alias{testDatasetAgeing}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Method "testDatasetAgeing"
}
\description{
The driver function in the PhenStatAgeing package for running statisical analysis on phenotypic data. \cr

It performs several checks on the data and the input model before running the analysis. This function supports three main analysis frameworks precisely, Linear Mix model (MM), Fisher Exact test (FE) and Reference Range plus (RR).\cr 

It further monitors the process for failures and errors and applies some runtime patches/fixes.\cr 

The function parameters are designed to be human-friendly by initialising the inputs by the model that will be applied to the data.

}
\usage{
testDatasetAgeing(
	phenListAgeing = NULL,
	method         = NULL,
	MM_fixed       = TypicalModel(
		depVariable = "data_point",
		withWeight  = MM_BodyWeightIncluded,
		Sex         = TRUE,
		LifeStage   = TRUE,
		data   = phenListAgeing@datasetPL,
		others = NULL,
		debug  = debug
	),
	MM_random             = rndProce("TYPICAL"),
	MM_BodyWeightIncluded = TRUE,
	MM_lower              = ~ Genotype + 1,
	MM_weight = if (TermInModelAndnLevels(
		model = MM_fixed, 
		data = phenListAgeing@datasetPL))
		varIdent(form = ~ 1 |
						 	LifeStage)
	else
		varIdent(form = ~ 1 |
						 	Genotype),
	MM_direction = "both",
	MM_checks    = c(TRUE, TRUE, TRUE),
	MM_optimise  = c(TRUE, TRUE, TRUE),
	FE_formula   = category ~ Genotype + Sex + LifeStage,
	RR_formula   = data_point ~ Genotype + Sex + LifeStage,
	RR_prop      = 0.95,
	FERR_rep     = 1500,
	debug        = TRUE,
	...
)
}
\arguments{
  \item{phenListAgeing}{
	 An instance of the PhenListAgeing or PhenList class; mandatory argument.
}
  \item{method}{
		A character string ("MM", "FE" or "RR") defining the method to use for model building.
}
  \item{MM_fixed}{
  Only applies to the "MM" framework. A formula that specifies the fixed effect in the linear mix model. The default is data_point~Genotype+Sex+LifeStage +/- BodyWeight. Note that the algorithm checks and removes the formula terms that do not exist in data. 
}
  \item{MM_random}{
  Only applies to the "MM" framework. The random effect in the linear mixed model. See the lme() function for more details about the formula. The default is ~1|Batch.
}
  \item{MM_BodyWeightIncluded}{
  Only applies to the default MM_fixed in the "MM" framework. If TRUE then the default model includes the body weight and the model would be: data_point~Genotype+Sex+LifeStage + Weight.
}
  \item{MM_lower}{
  Only applies to the "MM" framework. A right-sided formula, for example ~Genotype+1 or ~Sex+Genotype+1 or ~Sex+Genotype+Sex:Genotype. The lowest model that must not be included in the model optimisation. In other words, the terms in this model won't be removed during the optimisation process. The default is ~ Genotype + 1 that is the genotype effect and the intercept will be kept in the model during the optimisation process.
}
  \item{MM_weight}{
  Only applies to the "MM" framework. From weight in the lme() manual: \cr
  "an optional varFunc object or one-sided formula describing the within-group heteroscedasticity structure. If given as a formula, it is used as the argument to varFixed, corresponding to fixed variance weights. See the documentation on varClasses for a description of the available varFunc classes. Defaults to NULL, corresponding to homoscedastic within-group errors". \cr\cr
  The default is varIdent(form = ~ 1 |	LifeStage) if the LifeStage includes in the data. Otherwise	varIdent(form = ~ 1 |	Genotype).
}
  \item{MM_direction}{
	Only applies to the "MM" framework. Select from  "both" (for stepwise optimisation), "backward" (for backward elimination) or "forward" (for forward selection) for the optimisation algorithm. The default is "both".
	}
  \item{MM_checks}{
	Only applies to the "MM" framework. A vector of three 1/0 or TRUE/FALSE values such as c(TRUE, TRUE, TRUE). The first element of the vector activates checking the model terms (See MM_fixed) to be existed in data. The second term removes any single level factor from the model (in MM_fixed). The third element checks the interaction term to make sure all interactions have some data attached to them. Caution is needed for the third element of the check as it may take longer than usual if the formula in MM_fixed contains many factors. The default is c(TRUE, TRUE, TRUE).
}
	\item{MM_optimise}{
	Only applies to the "MM" framework. A vector of three binary values such as c(1,1,1) or c(TRUE,TRUE,TRUE), default. The first element of the vector activates the fixed effect optimisation. The algorithm uses AICc to optimise the fixed effects (Check AICcmodavg package for more details about AICc). The second and third elements of the vector activate optimisation on weight and random effects respectively. The optimisation of weight and random effects refers to comparing the AICc between a model with and without those effects.
	}
  \item{FE_formula}{
  Only applies to the "FE" framework. The model for analysing the categorical data. The default is: category ~ Genotype + Sex + LifeStage. \cr
  Note that similar to MM_fixed, the terms that do not exist in data or the interaction terms will be dismissed from the model. 
}
  \item{RR_formula}{
	Only applies to the "RR" framework. The model for analysing the RR+ compatible data. The default is: data_point ~ Genotype + Sex + LifeStage. \cr
  Note that similar to MM_fixed,  the terms that do not exist in data or the interaction terms will be dismissed from the model. 
}
  \item{RR_prop}{
	Only applies to the "RR" framework. A single value between (0.5,1) not including the boundaries. The threshold for the variation ranges in the RR framework. The default value is 0.95. 
}
  \item{FERR_rep}{
  Only applies to the "RR" or "FE" framework. The number of repetition for the Bayesian fisher exact test. See "B" parameter in fisher.test() function. Set to 0 for non-bayesian results (not recommended). The default is 1500.
}
  \item{debug}{
	A logical flag. Set to TRUE to see more details about the progress of the function. Default TRUE
}
  \item{\dots}{
  Other parameters that can be passed to:\cr
  If the model is set to Linear Mixed Model (MM) then the parameters that can be passed to lme function. See ?lme() manual page\cr
  If the model is either Fisher Exact test (FE) or Reference Range + (RR) then the parameters that can be passed to the fisher.test function. See ?fisher.test() manual page\cr
}
}

\value{
A successful execution of the function will return a list of three elements:
 \item{input }{This contains the list of inputs}
 \item{output }{A list of outputs}
 \item{extra }{A placeholder for extra information if exists}
 \item{* messages}{A placeholder for the errors/warnings in the case of failure}
}
\author{
Hamed Haseli Mashhadi <hamedhm@ebi.ac.uk>
}
\seealso{
\code{\link{PhenListAgeing}}, \code{\link{plot}}, \code{\link{summary}}
}
\examples{
		library(PhenStat)
		library(nlme)
		####################################################################
		# 1 Data preparation
		####################################################################
		
		#################
		# 1.1 Continuous data
		file = system.file("extdata", "test1.csv", package = "PhenStat")
		test_Cont = PhenList(dataset = read.csv(file, na.strings = '-'), 
		 testGenotype = "Sparc/Sparc")
		
		#################
		# 1.2 Categorical data
		file  = system.file("extdata", "test_categorical.csv", package = "PhenStat")
		test_Cat = PhenList(dataset = read.csv(file, na.strings = '-'), 
		 testGenotype = "Aff3/Aff3")
		
		####################################################################
		# 2 Testing frameworks
		####################################################################
		
		#################
		# 2.1 Optimised Mixed model (MM) framework
		#################
		MM1_result = testDatasetAgeing(
			phenListAgeing = test_Cont,
			method = 'MM',
			MM_fixed = Lean.Mass ~ Genotype + Sex + Weight
		)
		VO_MM1 = vectorOutputAgeing(MM1_result)
		plot(MM1_result, col = 2,main = 'Optimised model')
		summary(MM1_result)
		
		
		#################
		# 2.2 Mixed model (MM) with NO optimisation 
		# for the fixed effects but random/weight effects
		#################
		MM2_result = testDatasetAgeing(
			phenListAgeing = test_Cont,
			method = 'MM',
			MM_fixed = Lean.Mass ~ Genotype + Sex + Weight,
			MM_lower =           ~ Genotype + Sex + Weight + 1
			# Or simply MM_optimised = c(0, 1, 1)
		)
		VO_MM2 = vectorOutputAgeing(MM2_result)
		plot(MM2_result, col = 4,main = 'No optimisation on the fixed effects')
		summary(MM2_result)
		
		#################
		# 2.3 Mixed model (MM) with NO optimisation at all
		#################
		MM3_result = testDatasetAgeing(
			phenListAgeing = test_Cont,
			method = 'MM',
			MM_fixed = Lean.Mass ~ Genotype + Sex + Weight,
			MM_optimise = c(0,0,0)
		)
		VO_MM3 = vectorOutputAgeing(MM3_result)
		plot(MM3_result, col = 5,main = 'Not optimised model')
		summary(MM3_result)
		
		
		#################
		# 2.3 Reference range model
		#################
		RR_result = testDatasetAgeing(
			phenListAgeing = test_Cont,
			method = 'RR',
			RR_formula  = Lean.Mass ~ Genotype + Sex
		)
		VO_RR = vectorOutputAgeing(RR_result)
		plot(RR_result, col = 3:4)
		summary(RR_result)
		
		
		#################
		# 2.4 Fisher exact model
		#################
		FE_result = testDatasetAgeing(
			phenListAgeing = test_Cat,
			method = "FE",
			FE_formula =  Thoracic.Processes ~ Genotype + Sex + Weight
		)
		VO_FE = vectorOutputAgeing(FE_result)
		plot(FE_result, col = 1:2)
		summary(FE_result)
		
}
\keyword{ ~PhenStatAgeing }
