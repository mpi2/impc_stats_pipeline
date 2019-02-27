testDatasetAgeing = function(procedure = 'TYPICAL',
														 phenListAgeing = NULL,
														 depVariable = NULL   ,
														 withWeight = FALSE,
														 Sex = TRUE,
														 LifeStage = TRUE,
														 fixed = TypicalModel(depVariable,
														 										 withWeight  ,
														 										 Sex,
														 										 LifeStage,
														 										 data = phenListAgeing@datasetPL),
														 random = rndProce (procedure),
														 weight = if (LifeStage)
														 	varIdent(form =  ~ 1 |
														 					 	LifeStage)
														 else
														 	varIdent(form =  ~ 1 |	Genotype),
														 conversion.threshold = Inf,
														 family = poisson(link = "log"),
														 direction = 'both',
														 rep = 10 ^ 5) {
	r = tryCatch(
		expr = {
			testDatasetAgeing0(
				procedure = procedure,
				phenListAgeing = phenListAgeing,
				depVariable = depVariable   ,
				withWeight = withWeight,
				Sex = Sex,
				LifeStage = LifeStage,
				fixed = fixed,
				random = random,
				conversion.threshold = conversion.threshold,
				family = family,
				direction = direction,
				rep = rep,
				weight = weight
			)
		},
		warning = function(war) {
			warning('\n ~> Warning : ', war, '\n')
			return(NULL)
		},
		error = function(err) {
			warning('\n ~> error : ', err, '\n')
			return(NULL)
		}
	)
	return(r)
}

# Core
testDatasetAgeing0 = function(procedure = 'TYPICAL',
															phenListAgeing = NULL,
															depVariable = NULL   ,
															withWeight = FALSE,
															Sex = TRUE,
															LifeStage = TRUE,
															fixed ,
															random ,
															weight ,
															direction = 'both',
															conversion.threshold = Inf,
															family = poisson(link = "log"),
															rep = 10 ^ 5) {
	if (!is(phenListAgeing, 'PhenListAgeing'))
		stop('\n ~> function requires a "PhenListAgeing" object \n')
	if (noVariation(data = phenListAgeing@datasetPL))
		stop('\n ~> There is no variation on Genotype.\n')
	if (LifeStage &&
			noVariation(data = phenListAgeing@datasetPL, f = '~ Genotype+LifeStage'))
		stop('\n ~> There is no variation on a cell in the Genotype*LifeStage table\n')
	
	chklist = columnChecks0(
		dataset = phenListAgeing@datasetPL,
		columnName = depVariable,
		dataPointsThreshold = 4
	)
	categorical = IsCategorical(phenListAgeing@datasetPL[, depVariable], max.consider = conversion.threshold) # Automatic type detection
	if (all(chklist[1:2])) {
		# Important to remove missings
		Obj =  CheckMissing(phenListAgeing@datasetPL, fixed)
		data = Obj$new.data
		# All combination of the variables
		AllCombinations = AllTables(dframe = data[, all.vars(fixed)],
																response.name = depVariable,
																shrink = FALSE)
		o.Model = M.opt(
			fixed  = fixed,
			random = random,
			data   = data,
			direction = direction,
			family = family,
			categorical = categorical$cat,
			LifeStage = LifeStage,
			weight = weight,
			trace = FALSE
		)
		if (LifeStage) {
			effA = eff.size(
				object = o.Model$Final.Model,
				depVariable = depVariable,
				effOfInd = 'LifeStage'
			)
			effE = eff.size(
				object = o.Model$Final.Model,
				data = subset(data, LifeStage == 'Early'),
				depVariable = depVariable
			)
			effL = eff.size(
				object = o.Model$Final.Model,
				data = subset(data, LifeStage == 'Late'),
				depVariable = depVariable
			)
			effs = list(Early = effE,
									Late = effL,
									LifeStage = effA)
		} else{
			effs = eff.size(object = o.Model$Final.Model,
											depVariable = depVariable)
		}
		####
		output =    list(
			initial.model   = o.Model$Initial.Model,
			final.model     = o.Model$Final.Model  ,
			split.models    = o.Model$SplitModels  ,
			fixed           = fixed                ,
			random          = random               ,
			fin.fixed       = formula(o.Model$Final.Model),
			fin.random      = if (is.null(o.Model$Final.Model$call$random)) {
				NULL
			} else{
				o.Model$Final.Model$call$random
			},
			VarHomo         = o.Model$VarHomo,
			weight          = weight,
			effect          = effs,
			data            = phenListAgeing,
			OriginalDataSet = Obj$org.data,
			depVariable     = depVariable,
			procedure       = procedure,
			Sex             = Sex,
			LifeStage       = LifeStage,
			categorical     = categorical$cat,
			method          = categorical$method,
			direction       = direction,
			missings        = Obj$missings,
			AllCombinations   = AllCombinations
		)
		class(output) = 'PhenListAgeingModel'
		
	} else if (chklist[1] && !chklist[2]) {
		output = crunner(
			object = phenListAgeing,
			Sex = Sex,
			LifeStage = LifeStage,
			depVariable = depVariable,
			rep = rep,
			categorical = categorical$cat,
			method      = categorical$method
		)
		class(output) = 'PhenlistCategoricalAgeingModel'
	} else{
		stop(
			'\n ~> dependent variable must be (1) numeric/nominal (2) exist in the data frame (3) at least 4 datapoints length \n'
		)
	}
	return(output)
}
