testDatasetAgeing = function(phenListAgeing = NULL                 ,
														 method         = NULL                 ,
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
														 MM_checks    = c(TRUE, TRUE, TRUE),
														 MM_optimise  = c(TRUE, TRUE, TRUE),
														 ### FE or RR
														 FE_formula = category   ~ Genotype + Sex + LifeStage,
														 RR_formula = data_point ~ Genotype + Sex + LifeStage,
														 RR_prop    = 0.95,
														 FERR_rep   = 1500,
														 ##### Others
														 MMFERR_conf.level = 0.95 ,
														 seed       = NULL        ,
														 debug      = TRUE        ,
														 ...) {
	if(!is.null(seed))
		set.seed(seed)
	r = NULL
	s = tryCatch(
		expr = {
			suppressMessagesANDWarnings(
				testDatasetAgeing0(
					phenListAgeing = phenListAgeing,
					method = method ,
					#### MM
					MM_fixed = MM_fixed                   ,
					MM_random = MM_random                 ,
					MM_lower = MM_lower                   ,
					MM_weight = MM_weight                 ,
					MM_direction = MM_direction           ,
					MM_checks   = MM_checks               ,
					MM_optimise = MM_optimise             ,
					#### FE     
					FE_formula = FE_formula               ,
					RR_formula = RR_formula               ,
					RR_prop = RR_prop                     ,
					FERR_rep = FERR_rep                   ,
					debug = debug                         ,
					MMFERR_conf.level = MMFERR_conf.level ,
					...
				),
				sup.messages = !debug,
				sup.warnings = TRUE
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
			message(err)
			return(NULL)
		}
	)
	if (is.null(s))
		return(r)
	else
		return(s)
}

testDatasetAgeing0 = function(phenListAgeing = NULL     ,
															method                    ,
															#### MM   
															MM_fixed                  ,
															MM_random                 ,
															MM_lower                  ,
															MM_weight                 ,
															MM_direction = 'both'     ,
															MM_checks                 ,
															MM_optimise               ,
															##### FE or RR   
															FE_formula                ,
															RR_formula                ,
															RR_prop                   ,
															FERR_rep                  ,
															MMFERR_conf.level = 0.95  ,
															##### Others
															debug = TRUE              ,
															...) {
	message0('PhenStatAgeing loaded.')
	if(is.null(phenListAgeing))
		stop('\n ~> The input dataset cannot be NULL')
	if (!is(phenListAgeing, 'PhenList') &&
			!is(phenListAgeing, 'PhenListAgeing'))
		stop('\n ~> function expects "PhenList" or "PhenListAgeing" object \n')
	if (noVariation(data = phenListAgeing@datasetPL))
		stop('\n ~> There is no variation in Genotype.\n')
	if (MMFERR_conf.level >= 1   ||
			MMFERR_conf.level <= 0   ||
			length(MMFERR_conf.level) > 1)
		stop('\n ~> Confidence level must be a single value in (0,1) interval')
	if (length(method) > 1 || !method %in% c('MM', 'FE', 'RR'))
		stop('\n ~> `method` must be one of `MM`, `FE` or `RR`')
	
	if (method %in% 'MM') {
		message0('Linear Mixed Model (MM framework) in progress ...')
		output =    M.opt(
			fixed     = MM_fixed          ,
			random    = MM_random         ,
			object    = phenListAgeing    ,
			lower     = MM_lower          ,
			direction = MM_direction      ,
			weight    = MM_weight         ,
			checks    = MM_checks         ,
			optimise  = MM_optimise       ,
			trace     = FALSE             ,
			method    = 'MM'              ,
			ci_levels = MMFERR_conf.level ,
			...
		)
		# Important!
		if (!is.null(output$input))
			output$input$fixed = MM_fixed
		
	} else if (method %in% 'FE') {
		message0('Fisher Exact Test (FE framework) in progress ...')
		output = crunner(
			object = phenListAgeing                                                    ,
			formula = MoveResponseToRightOfTheFormula(FE_formula)                      ,
			#expandDottedFormula(formula = FE_formula, data = phenListAgeing@datasetPL),
			rep = FERR_rep                                                             ,
			method = 'FE'                                                              ,
			fullComparisions = TRUE                                                    ,
			ci_levels = MMFERR_conf.level                                              ,
			...
		)
		# Important!
		if (!is.null(output$input))
			output$input$formula = FE_formula
		
	} else if (method %in% 'RR') {
		message0('Reference Range Plus (RR framework) in progress ...')
		output = RRrunner(
			object = phenListAgeing                                                     ,
			formula = MoveResponseToRightOfTheFormula(RR_formula)                       ,
			#expandDottedFormula(formula = RR_formula, data = phenListAgeing@datasetPL) ,
			rep = FERR_rep                                                              ,
			method = 'RR'                                                               ,
			RRprop = RR_prop                                                            ,
			ci_levels = MMFERR_conf.level                                               ,
			...
		)
		# Important!
		if (!is.null(output$input))
			output$input$formula = RR_formula
	} else{
		message0('No "method" is specified. ')
		output = NULL
	}
	return(output)
}
