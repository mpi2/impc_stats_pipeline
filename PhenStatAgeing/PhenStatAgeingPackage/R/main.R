testDatasetAgeing = function(phenListAgeing = NULL                        ,
														 method         = NULL                        ,
														 ### MM
														 MM_fixed = TypicalModel(
														 	depVariable = 'data_point'           ,
														 	withWeight = MM_BodyWeightIncluded   ,
														 	Sex = TRUE                           ,
														 	LifeStage = TRUE                     ,
														 	data = phenListAgeing@datasetPL      ,
														 	others = NULL                        ,
														 	debug = debug
														 )                                     ,
														 MM_random = rndProce ('TYPICAL')      ,
														 MM_BodyWeightIncluded = TRUE          ,
														 MM_lower  =  ~ Genotype + 1           ,
														 MM_weight = if (TermInModelAndnLevels(model = MM_fixed,
														 																			data = phenListAgeing@datasetPL))
														 	varIdent(form =  ~ 1 |
														 					 	LifeStage)
														 else
														 	varIdent(form =  ~ 1 |	Genotype),
														 MM_direction = 'both',
														 MM_checks    = c(1, 1, 1),
														 ### FE or RR
														 FE_formula = category ~  Genotype + Sex + LifeStage,
														 RR_formula = data_point ~ Genotype + Sex + LifeStage,
														 RR_prop    = 0.95,
														 FERR_rep   = 1500,
														 ##### Others
														 debug      = TRUE,
														 ...) {
	r = NULL
	s = tryCatch(
		expr = {
			suppressMessagesANDErrors(
				testDatasetAgeing0(
					phenListAgeing = phenListAgeing,
					method = method ,
					#### MM
					MM_fixed = MM_fixed,
					MM_random = MM_random,
					MM_lower = MM_lower,
					MM_weight = MM_weight,
					MM_direction = MM_direction,
					MM_checks   = MM_checks,
					#### FE
					FE_formula = FE_formula,
					RR_formula = RR_formula,
					RR_prop = RR_prop,
					FERR_rep = FERR_rep,
					debug = debug
				),
				debug = debug
			)
		},
		warning = function(war) {
			message0('The functions failed with a warning (see below): ')
			r$messages$warning <<- war
			warning(war)
			return(NULL)
		},
		error = function(err) {
			message0('The functions failed with an error (see below): ')
			r$messages$error <<- err
			warning(err)
			return(NULL)
		}
	)
	return(c(s,r))
}



testDatasetAgeing0 = function(phenListAgeing = NULL ,
															method                ,
															#### MM
															MM_fixed               ,
															MM_random              ,
															MM_lower               ,
															MM_weight              ,
															MM_direction = 'both'  ,
															MM_checks              ,
															##### FE or RR
															FE_formula,
															RR_formula,
															RR_prop,
															FERR_rep,
															##### Others
															debug = TRUE) {
	if (!is(phenListAgeing, 'PhenList') &&
			!is(phenListAgeing, 'PhenListAgeing'))
		stop('\n ~> function expects  "PhenList" or "PhenListAgeing" object \n')
	if (noVariation(data = phenListAgeing@datasetPL))
		stop('\n ~> There is no variation on Genotype.\n')
	
	
	if (method %in% 'MM') {
		output = M.opt(
			fixed = MM_fixed ,
			random = MM_random,
			object = phenListAgeing,
			lower = MM_lower,
			direction = MM_direction,
			weight = MM_weight,
			checks = MM_checks,
			trace = debug,
			method = 'MM'
		)
		# Important!
		if (!is.null(output))
			output$input$fixed = MM_fixed
		
	} else if (method %in% 'FE') {
		output = crunner(
			object = phenListAgeing,
			formula = MoveResponseToRightOfTheFormula(FE_formula),
			#expandDottedFormula(formula = FE_formula, data = phenListAgeing@datasetPL),
			rep = FERR_rep,
			method = 'FE',
			fullComparisions = TRUE
		)
		# Important!
		if (!is.null(output))
			output$input$formula = FE_formula
		
	} else if (method %in% 'RR') {
		output = RRrunner(
			object = phenListAgeing,
			formula = MoveResponseToRightOfTheFormula(RR_formula),
			#expandDottedFormula(formula = RR_formula, data = phenListAgeing@datasetPL),
			rep = FERR_rep,
			method = 'RR',
			RRprop = RR_prop
		)
		# Important!
		if (!is.null(output))
			output$input$formula = RR_formula
	} else{
		message0('No "method" is specified. ')
		output = NULL
	}
	return(output)
}
