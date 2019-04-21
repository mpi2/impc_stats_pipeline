# Continuous OPT core
M.opt = function(object = NULL            ,
								 fixed  = NULL            ,
								 random = NULL            ,
								 lower  = ~ Genotype + 1  ,
								 direction   = 'both'     ,
								 trace       = TRUE       ,
								 weight                   ,
								 checks      = c(1, 1, 1) ,
								 method      = 'MM'       ,
								 optimise    = TRUE       ,
								 ...) {
	if (!method %in% c('MM')        ||
			is.null(all_vars0(fixed))   ||
			is.null(object)) {
		#####
		message0 ('Improper method (',
							method,
							') for the type of data, or the "formula" is left blank')
		return(NULL)
	}
	message0('MM framework in progress ...')
	# Future devs
	family      = poisson(link = "log")
	categorical = optimised = FALSE
	# Start from here
	sta.time    = Sys.time()
	lCont       = lmeControl (opt  = 'optim',
														maxIter = 1000,
														msMaxIter = 1000)
	gCont       = glsControl (
		opt  = 'optim',
		singular.ok = TRUE,
		maxIter = 1000,
		msMaxIter = 1000
	)
	glCont      = glm.control(epsilon = 10 ^ -36, maxit = 1000)
	G.Model     = FV.Model  = I.Model = SplitModels = F.Model = OutR = NULL
	VarHomo     = TRUE
	data        = object@datasetPL
	n           = nrow(data)
	fixed = ModelChecks(fixed = fixed,
											data = data,
											checks = checks)
	allVars     = all.vars(fixed)
	LifeStage   = 'LifeStage' %in% allVars
	Batch_exist = !categorical && !is.null(random) &&
		colExists(name = 'Batch', data = data)
	mdl         = ifelse(Batch_exist, 'lme', ifelse(categorical, 'glm', 'gls'))
	message0(mdl, ': Fitting the full model ... ')
	# Just for rounding the full model errors
	initialFixed = fixed
	fixedTerms   = formulaTerms(initialFixed)
	for (i in 1:length(fixedTerms)) {
		I.Model    = tryCatch(
			expr     = do.call(mdl,
												 listFun(
												 	list = list(
												 		model     = if (mdl == 'glm') {
												 			TRUE
												 		} else{
												 			fixed
												 		},
												 		fixed     = fixed  ,
												 		formula   = fixed  ,
												 		family    = family ,
												 		random    = random ,
												 		data      = data   ,
												 		na.action = na.omit,
												 		method    = ifelse(mdl == 'glm', 'glm.fit', 'ML'),
												 		weights   = if (mdl != 'glm') {
												 			weight
												 		} else{
												 			NULL
												 		},
												 		control   = if (Batch_exist) {
												 			lCont
												 		} else{
												 			if (categorical)
												 				glCont
												 			else
												 				gCont
												 		},
												 		...
												 	),
												 	FUN = ifelse(Batch_exist, 'lme', ifelse(categorical, 'glm', 'gls')),
												 	debug = TRUE
												 )),
			warning = function(war) {
				message0('* The full model failed with the warning (see below): ')
				message0('\t', war, breakLine = FALSE)
				return(NULL)
			},
			error = function(err) {
				message0('* The full model failed with the error (see below): ')
				message0('\t', err, breakLine = FALSE)
				return(NULL)
			}
		)
		###########
		if (is.null(I.Model)) {
			removedTerm = fixedTerms[length(fixedTerms) -	(i - 1)]
			ltemp = ModelInReference(model = lower, reference = fixed)
			if (!is.null(lower) &&
					!is.null(fixed) &&
					!is.null(ltemp) &&
					removedTerm %in% formulaTerms(ltemp)) {
				message0(
					'The following term did not removed: ',
					removedTerm,
					'\n\t the terms below will not remove from the model:\n\t ',
					printformula(formulaTerms(ltemp))
				)
				next
			}
			fixed       = update.formula(old = initialFixed, new = paste0('.~.-', removedTerm))
			message0(
				'Round ',
				i,
				' of fixing the error. The followeing term will be removed: ',
				removedTerm,
				'\n\tNew formula: ',
				printformula(fixed)
			)
		} else{
			message0('The full model successfully applied.')
			break
		}
	}
	###########
	if (is.null(I.Model))
		message0('Full model failed ...')
	###########
	message0('The specified "lower" model: \n\t',
					 ifelse(!is.null(lower), printformula(lower), 'Null lower'))
	lowerCorrected = ModelInReference(model = lower, reference = fixed)
	###########
	message0('The model optimisation is ',
					 ifelse(optimise, 'in progress ...', 'is not activated'))
	if (optimise && !is.null(I.Model) && !is.null(lowerCorrected)) {
		message0('\tOptimising the model ... ')
		message0('\tThe direction of  the optimisation (backward, forward, both): ', direction)
		F.Model = tryCatch(
			expr = stepAIC0(
				I.Model                     ,
				trace = trace               ,
				direction = direction       ,
				scope = list(lower = lowerCorrected) ,
				na.action = na.omit                  ,
				k = log(n)
			),
			warning = function(war) {
				message0('\t * The optimisation failed with the warning (see below): ')
				message0('\t', war, breakLine = FALSE)
				return(NULL)
			},
			error = function(err) {
				message0('\t * The optimisation failed with the error (see below): ')
				message0('\t', err, breakLine = FALSE)
				return(NULL)
			}
		)
		###########
		if (is.null(F.Model)) {
			message0('\t Optimisation did not apply')
			F.Model = I.Model
		} else{
			message0('\t Optimised model: ', printformula(formula(F.Model)))
			optimised = TRUE
		}
	} else{
		F.Model = I.Model
	}
	###########
	if (!is.null(F.Model)) {
		if (optimise && !is.null(weight) && !(mdl %in% 'glm')) {
			message0('Testing varHom ... ')
			FV.Model = tryCatch(
				expr  = update(F.Model, weights = NULL),
				warning = function(war) {
					message0('* Testing VarHom failed with the warning (see below): ')
					message0('\t', war, breakLine = FALSE)
					return(NULL)
				},
				error = function(err) {
					message0('* Testing VarHom failed with the error (see below): ')
					message0('\t', err, breakLine = FALSE)
					return(NULL)
				}
			)
			if (!is.null(F.Model) &&
					!is.null(FV.Model) &&
					AIC(FV.Model, k = log(n)) <= AIC(F.Model, k = log(n))) {
				# = is important
				F.Model = FV.Model
				VarHomo = FALSE
				message0('\tVarHom checked out ... ')
			}
		}
		# Batch test
		if (optimise && Batch_exist && !(mdl %in% 'glm')) {
			message0('Testing Batch ... ')
			G.Model =  tryCatch(
				expr = do.call(
					'gls',
					args        = list(
						model     = formula(F.Model)    ,
						data      = data                ,
						na.action = na.omit             ,
						method    = "ML"                ,
						weights   = F.Model$call$weights,
						control   = gCont
					)
				),
				
				warning = function(war) {
					message0('* Testing Batch failed with the warning (see below): ')
					message0('\t', war, breakLine = FALSE)
					return(NULL)
				},
				error = function(err) {
					message0('* Testing Batch failed with the error (see below): ')
					message0('\t', err, breakLine = FALSE)
					return(NULL)
				}
			)
			if (!is.null(G.Model) &&
					!is.null(F.Model) &&
					AIC(G.Model, k = log(n)) <= AIC(F.Model, k = log(n))) {
				F.Model = G.Model
				message0('\tBatch checked out ... ')
			}
		}
		SplitModels = SplitEffect(
			finalformula = formula(F.Model),
			fullModelFormula = fixed,
			F.Model = F.Model,
			data = data,
			depVariable = allVars[1]
		)
		message0('Estimating effect sizes ... ')
		EffectSizes = c(
			suppressMessages(
				AllEffSizes(
					object = F.Model,
					depVariable = allVars[1],
					effOfInd    = allVars[-1],
					data = data
				)
			),
			CombinedEffectSizes = lapply(SplitModels, function(x) {
				percentageChangeCont(
					model = x,
					data = getData(x),
					variable = NULL,
					depVar = allVars[1],
					individual = FALSE,
					mainEffsOnlyWhenIndivi = x$MainEffect
				)
			})
		)
		message0(
			'\tTotal effect sizes estimated: ',
			ifelse(
				!is.null(EffectSizes),
				length(EffectSizes),
				'Not possible because of errors!'
			)
		)
		message0('Quality tests in progress ... ')
		ResidualNormalityTest = QuyalityTests(F.Model, levels = allVars[-1])
	} else{
		message0('This process fully terminated. No success in recovering the initial model.')
		OutR$messages  = 'This process fully terminated. No success in recovering the initial model.'
		return(OutR)
	}
	message0('MM framework executed in ', round(difftime(Sys.time() , sta.time, units = 'sec'), 2), ' seconds')
	####
	OutR = list(
		output = list(
			Final.Model         = F.Model                  ,
			Initial.Model       = I.Model                  ,
			EffectSizes         = EffectSizes              ,
			ResidualNormalityTests = ResidualNormalityTest ,
			NoVarStr.Model      = FV.Model                 ,
			NoBatch.Model       = G.Model                  ,
			SplitModels         = SplitModels              ,
			Final.Model.Tag     = class(F.Model)[1]        ,
			Initial.Model.Tag   = mdl                      ,
			VarHomoIn           = VarHomo                  ,
			BatchIn             = Batch_exist  &&   !(mdl %in% 'glm')    ,
			SexIn               = termInTheModel(
				model = formula(F.Model)                                   ,
				term = 'Sex'                                               ,
				message = FALSE
			),
			LifeStageIn         = termInTheModel(
				model = formula(F.Model)                                   ,
				term = 'LifeStage'                                         ,
				message = FALSE
			)                                                            ,
			optimised = optimised
		)                                                              ,
		input = list(
			PhenListAgeing      = object                                 ,
			fixed               = initialFixed                           ,
			random              = random                                 ,
			data                = data                                   ,
			depVariable         = allVars[1]                             ,
			lower               = lower                                  ,
			direction           = direction                              ,
			LifeStage           = LifeStage                              ,
			weight              = weight                                 ,
			checks              = checks                                 ,
			optimise            = optimise                               ,
			others              = ...
		),
		extra = list(Cleanedformula  = fixed,
								 lowerCorrected  = lowerCorrected)
	)
	class(OutR) <- 'PhenStatAgeingMM'
	return(OutR)
}
