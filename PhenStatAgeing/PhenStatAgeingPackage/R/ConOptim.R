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
	categorical = FALSE
	# Start from here
	sta.time    = Sys.time()
	lCont       = lmeControl (opt  = 'optim')
	gCont       = glsControl (opt  = 'optim', singular.ok = TRUE)
	glCont      = glm.control(epsilon = 10 ^ -36)
	G.Model     = FV.Model  = I.Model = SplitModels = F.Model = NULL
	VarHomo     = TRUE
	data        = object@datasetPL
	if (any(checks > 0)) {
		if (checks[1])
			fixed  = checkModelTermsInData(formula = fixed,
																		 data = data,
																		 responseIsTheFirst = TRUE)
		if (checks[2])
			fixed = removeSingleLevelFactors(formula = fixed, data = data)
		if (checks[3])
			fixed = ComplementaryFeasibleTermsInContFormula(formula = fixed, data = data)
		message0('Checked model: ', printformula(fixed))
	}
	allVars     = all.vars(fixed)
	LifeStage   = 'LifeStage' %in% allVars
	Batch_exist = !categorical &&
		colExists(name = 'Batch', data = data)
	mdl         = ifelse(Batch_exist, 'lme', ifelse(categorical, 'glm', 'gls'))
	message0(mdl, ': Fitting the full model ... ')
	I.Model = do.call(mdl,
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
										))
	message0('The specified "lower" model: ', printformula(lower))
	lowerCorrected = ModelInReference(model = lower, reference = fixed)
	if (!is.null(lowerCorrected)) {
		message0('Optimising the model ... ')
		F.Model = MASS::stepAIC(
			I.Model                     ,
			trace = trace               ,
			direction = direction       ,
			scope = list(lower = lowerCorrected) ,
			na.action = na.omit
		)
	} else{
		F.Model = I.Model
	}
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
	EffectSizes = suppressMessages(AllEffSizes(
		object = F.Model,
		depVariable = allVars[1],
		effOfInd    = allVars[-1],
		data = data
	))
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
			)
		)                                                              ,
		input = list(
			PhenListAgeing      = object                                 ,
			fixed               = fixed                                  ,
			random              = random                                 ,
			data                = data                                   ,
			depVariable         = allVars[1]                             ,
			lower               = lower                                  ,
			direction           = direction                              ,
			LifeStage           = LifeStage                              ,
			weight              = weight                                 ,
			checks              = checks
		),
		extra = list(
			Cleanedformula  = fixed,
			lowerCorrected  = lowerCorrected)
	)
	class(OutR) <- 'PhenStatAgeingMM'
	return(OutR)
}


# #### Summary core for continues data
# summary.PhenListAgeingModel = function(object, ...) {
# 	if (is.null(object))
# 		stop('\n ~> NULL object\n')
# 	tmpobject  = object
# 	
# 	ai = formulaTerms(object$initial.model)
# 	af = formulaTerms(object$final.model)
# 	if (is.null(object$final.model$call$random)) {
# 		ar = NULL
# 	} else{
# 		ar = object$final.model$call$random
# 	}
# 	#notSigVars = ai[!(ai %in% af)]	#SigVars    = ai[(ai %in% af)]
# 	tmpobject$final.model$call$data = tmpobject$final.model$call$family = NULL
# 	F.Sum = summary(tmpobject$final.model) #main
# 	message0(
# 		'~> Final model: ',
# 		'\n\t Fixed: ' ,
# 		ifelse(is.null(af), '-', paste(af, collapse = ', ')),
# 		'\n\t Random: ' ,
# 		ifelse(is.null(ar), '-', format(ar)),
# 		'\n\t ----- ' ,
# 		'\n\t Removed term(s): ',
# 		ifelse(length(ai[!(ai %in% af)]) < 1, '-', paste(ai[!(ai %in% af)], collapse = ', ')),
# 		'\n\t ----- ' ,
# 		'\n\t Genotype effect size: ',
# 		ifelse(
# 			object$LifeStage,
# 			paste(
# 				c(
# 					'\n\t  Genotype for Early ~>',
# 					'\n\t  Genotype for Late ~>' ,
# 					'\n\t  LifeStage ~>'
# 				),
# 				object$effect,
# 				collapse = ', '
# 			),
# 			object$effect
# 		),
# 		'\n\t ----- ' ,
# 		'\n\t Variance homogeneity: ',
# 		ifelse(!is.null(object$VarHomo), object$VarHomo, 'Not specified'),
# 		'\n\t ----- ' ,
# 		'\n\n~> Model summary:\n'
# 	)
# 	print(F.Sum)
# 	outp = list(
# 		Initial.Terms = ai   ,
# 		Final.Terms = af     ,
# 		RandoEffect = ar     ,
# 		removed.terms = ifelse(length(ai[!(ai %in% af)]) < 1, '-', paste(ai[!(ai %in% af)], collapse = ', ')),
# 		object  = object     ,
# 		Summary = summary(object$final.model)
# 	)
# 	outp$JSON = toJSONI(outp$Summary)
# 	return(invisible(outp))
# }
# 
# # Plot for PhenListAgeingModel object
# plot.PhenListAgeingModel = function (x,
# 																		 ask = FALSE,
# 																		 mfrow = c(2, 2),
# 																		 ...) {
# 	if (is.null(x))
# 		stop('\n ~> NULL object\n')
# 	p = par()
# 	par(ask = ask, mfrow = mfrow)
# 	
# 	predR  = 	predict(x$final.model)
# 	residR = resid(x$final.model)
# 	plot(predR ,
# 			 residR,
# 			 xlab = 'Fitted values',
# 			 ylab = 'Residuals',
# 			 ...)
# 	abline(h = 0)
# 	hist(residR,
# 			 xlab = 'Residuals',
# 			 main = 'Histogram of residuals',
# 			 probability = TRUE,
# 			 ...)
# 	qqnorm(residR, main = 'Normal Q-Q plot of residuals', ...)
# 	qqline(residR, ...)
# 	plot(density(getData(x$final.model)[, x$depVariable]), main = 'Density of the response', ...)
# 	par(ask = p$ask, mfrow = p$mfrow)
# }

