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
								 ...) {
	if (!method %in% c('MM')        ||
			is.null(all.vars0(fixed))   ||
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
	lCont       = lmeControl (opt  = 'optim')
	gCont       = glsControl (opt  = 'optim', singular.ok = TRUE)
	glCont      = glm.control(epsilon = 10 ^ -36)
	G.Model     = FV.Model  = I.Model = SplitModels = F.Model = NULL
	VarHomo     = TRUE
	data        = object@datasetPL
	fixed = ModelChecks(fixed = fixed,
											data = data,
											checks = checks)
	allVars     = all.vars(fixed)
	LifeStage   = 'LifeStage' %in% allVars
	Batch_exist = !categorical &&
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
												 		}
												 	),
												 	FUN = ifelse(Batch_exist, 'lme', ifelse(categorical, 'glm', 'gls')),
												 	debug = TRUE
												 )),
			warning = function(war) {
				message0('* The functions failed with a warning (see below): ')
				message0('\t', war)
				return(NULL)
			},
			error = function(err) {
				message0('* The functions failed with an error (see below): ')
				message0('\t', err)
				return(NULL)
			}
		)
		if (is.null(I.Model)) {
			removedTerm = fixedTerms[length(fixedTerms) -	(i - 1)]
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
			break
		}
	}
	message0('The specified "lower" model: \n\t', printformula(lower))
	lowerCorrected = ModelInReference(model = lower, reference = fixed)
	if (!is.null(lowerCorrected) && !is.null(I.Model)) {
		message0('Optimising the model ... ')
		F.Model = tryCatch(
			expr = stepAIC0(
				I.Model                     ,
				trace = trace               ,
				direction = direction       ,
				scope = list(lower = lowerCorrected) ,
				na.action = na.omit
			),
			warning = function(war) {
				message0('* The optimisation failed with a warning (see below): ')
				message0('\t', war)
				return(NULL)
			},
			error = function(err) {
				message0('* The optimisation failed with an error (see below): ')
				message0('\t', err)
				return(NULL)
			}
		)
		if (is.null(F.Model)) {
			message0('Optimisation did not apply')
			F.Model = I.Model
		} else{
			optimised = TRUE
		}
	} else{
		F.Model = I.Model
	}
	if (!is.null(F.Model)) {
		if (!is.null(weight) && !(mdl %in% 'glm')) {
			message0('Testing varHom ... ')
			FV.Model = update(F.Model, weights = NULL)
			if (AIC(FV.Model) <= AIC(F.Model)) {
				# = is important
				F.Model = FV.Model
				VarHomo = FALSE
				message0('VarHom checked out ... ')
			}
		}
		# Batch test
		if (Batch_exist && !(mdl %in% 'glm')) {
			message0('Testing Batch ... ')
			G.Model = do.call(
				'gls',
				args        = list(
					model     = formula(F.Model)    ,
					data      = data                ,
					na.action = na.omit             ,
					method    = "ML"                ,
					weights   = F.Model$call$weights,
					control   = gCont
				)
			)
			if (AIC(G.Model) <= AIC(F.Model)) {
				F.Model = G.Model
				message0('Batch checked out ... ')
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
		EffectSizes = suppressMessages(
			AllEffSizes(
				object = F.Model,
				depVariable = allVars[1],
				effOfInd    = allVars[-1],
				data = data
			)
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
		return(NULL)
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
			checks              = checks
		),
		extra = list(Cleanedformula  = fixed,
								 lowerCorrected  = lowerCorrected)
	)
	class(OutR) <- 'PhenStatAgeingMM'
	return(OutR)
}
