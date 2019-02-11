# Continues OPT core
M.opt = function(fixed  = NULL,
								 random = NULL,
								 data   = NULL,
								 lower  = ~ Genotype + 1 ,
								 direction   = 'both',
								 family      = poisson(link = "log"),
								 trace       = TRUE ,
								 categorical = FALSE ,
								 weight              ,
								 PhenListAgeing      ,
								 ...) {
	sta.time    = Sys.time()
	lCont       = lmeControl (opt  = 'optim')
	gCont       = glsControl (opt  = 'optim', singular.ok = TRUE)
	glCont      = glm.control(epsilon = 10 ^ -22)
	G.Model     = FV.Model  = I.Model = SplitModels = F.Model = NULL
	VarHomo     = TRUE
	Inifixed    = fixed
	if (!is.null(fixed$correctd))
		fixed      = fixed$correctd
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
	lower = ModelInReference(model = lower, reference = fixed)
	if (!is.null(lower)) {
		message0('Optimising the model ... ')
		F.Model = MASS::stepAIC(
			I.Model                     ,
			trace = trace               ,
			direction = direction       ,
			scope = list(lower = lower) ,
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
	message0('Estimate effect sizes ... ')
	EffectSizes = suppressMessages(AllEffSizes(
		object = F.Model,
		depVariable = allVars[1],
		effOfInd    = allVars[-1],
		data = data
	))
	message0('Quality tests in progress ... ')
	ResidualNormalityTest = QuyalityTests(F.Model, levels = allVars[-1])
	message0('Finished (', round(Sys.time() - sta.time, 4), ' seconds )')
	return(list(
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
			fixed               = fixed                                  ,
			random              = random                                 ,
			data                = data                                   ,
			lower               = lower                                  ,
			direction           = direction                              ,
			family              = family                                 ,
			trace               = trace                                  ,
			categorical         = categorical                            ,
			LifeStage           = LifeStage                              ,
			weight              = weight                                 ,
			Fullfixed           = Inifixed                               ,
			PhenListAgeing      = PhenListAgeing
		),
		extra = list()
	))
}




#### Summary core for continues data
summary.PhenListAgeingModel = function(object, ...) {
	if (is.null(object))
		stop('\n ~> NULL object\n')
	tmpobject  = object
	ai = attr(terms(formula(object$initial.model)), "term.labels")
	af = attr(terms(formula(object$final.model)),   "term.labels")
	if (is.null(object$final.model$call$random)) {
		ar = NULL
	} else{
		ar = object$final.model$call$random
	}
	#notSigVars = ai[!(ai %in% af)]	#SigVars    = ai[(ai %in% af)]
	tmpobject$final.model$call$data = tmpobject$final.model$call$family = NULL
	F.Sum = summary(tmpobject$final.model) #main
	message0(
		'~> Final model: ',
		'\n\t Fixed: ' ,
		ifelse(is.null(af), '-', paste(af, collapse = ', ')),
		'\n\t Random: ' ,
		ifelse(is.null(ar), '-', format(ar)),
		'\n\t ----- ' ,
		'\n\t Removed term(s): ',
		ifelse(length(ai[!(ai %in% af)]) < 1, '-', paste(ai[!(ai %in% af)], collapse = ', ')),
		'\n\t ----- ' ,
		'\n\t Genotype effect size: ',
		ifelse(
			object$LifeStage,
			paste(
				c(
					'\n\t  Genotype for Early ~>',
					'\n\t  Genotype for Late ~>' ,
					'\n\t  LifeStage ~>'
				),
				object$effect,
				collapse = ', '
			),
			object$effect
		),
		'\n\t ----- ' ,
		'\n\t Variance homogeneity: ',
		ifelse(!is.null(object$VarHomo), object$VarHomo, 'Not specified'),
		'\n\t ----- ' ,
		'\n\n~> Model summary:\n'
	)
	print(F.Sum)
	outp = list(
		Initial.Terms = ai   ,
		Final.Terms = af     ,
		RandoEffect = ar     ,
		removed.terms = ifelse(length(ai[!(ai %in% af)]) < 1, '-', paste(ai[!(ai %in% af)], collapse = ', ')),
		object  = object     ,
		Summary = summary(object$final.model)
	)
	outp$JSON = toJSONI(outp$Summary)
	return(invisible(outp))
}

# Plot for PhenListAgeingModel object
plot.PhenListAgeingModel = function (x,
																		 ask = FALSE,
																		 mfrow = c(2, 2),
																		 ...) {
	if (is.null(x))
		stop('\n ~> NULL object\n')
	p = par()
	par(ask = ask, mfrow = mfrow)
	
	predR  = 	predict(x$final.model)
	residR = resid(x$final.model)
	plot(predR ,
			 residR,
			 xlab = 'Fitted values',
			 ylab = 'Residuals',
			 ...)
	abline(h = 0)
	hist(residR,
			 xlab = 'Residuals',
			 main = 'Histogram of residuals',
			 probability = TRUE,
			 ...)
	qqnorm(residR, main = 'Normal Q-Q plot of residuals', ...)
	qqline(residR, ...)
	plot(density(getData(x$final.model)[, x$depVariable]), main = 'Density of the response', ...)
	par(ask = p$ask, mfrow = p$mfrow)
}


AllEffSizes = function(object, depVariable, effOfInd, data) {
	lst       = flst = olst = NULL
	effOfInd  = effOfInd[effOfInd %in% names(data)]
	cats      = effOfInd[sapply(
		effOfInd,
		FUN = function(cl) {
			!is.numeric(data[, cl])
		}
	)]
	counter1 = counter2  = 1
	if (length(effOfInd)) {
		for (eff in effOfInd) {
			# 1. Main effect
			lstTmp = eff.size(
				object = object,
				depVariable = depVariable,
				effOfInd = eff,
				errorReturn = NULL
			)
			if (!is.null(lstTmp)) {
				lst[[counter1]] = lstTmp
				names(lst)[counter1] = eff
				counter1  = counter1  + 1
			}
			# 2. All subset interactions effect sizes
			if (length(cats) > 1) {
				for (j in 1:length(cats)) {
					cmbn = combn(x = cats, m = j)
					for (k in 1:ncol(cmbn)) {
						interact = interaction(data[, cmbn[, k] , drop = FALSE], sep = '.', drop = TRUE)
						for (lvl in levels(interact)) {
							olstTmp = eff.size(
								object = object,
								data = data[interact %in% lvl,],
								depVariable = depVariable,
								errorReturn = NULL,
								effOfInd = eff
							)
							if (!is.null(olstTmp)) {
								olst[[counter2]] = olstTmp
								names(olst)[counter2] = paste(eff, lvl, sep = '_')
								counter2 = counter2 + 1
							}
						}
					}
				}
			}
		}
		flst     = as.list(c(lst, olst))
		return(flst)
	} else{
		return(NULL)
	}
}

as.list0 = function(x, ...) {
	if (!is.null(x)) {
		return(as.list(x, ...))
	} else{
		return(NULL)
	}
}


QuyalityTests = function(object,
												 levels    = c('Genotype', 'Sex', 'LifeStage'),
												 list      = TRUE,
												 noDataLab = 'No data',
												 sep       = '_',
												 collapse  = '_') {
	if (is.null(object))
		return(NULL)
	r = resid(object)
	d = getData(object)
	levels = levels[levels %in% names(d)]
	levels = levels[sapply(levels, function(lvs) {
		!is.numeric(d[, lvs])
	})]
	counter = 1
	flst    = NULL
	if (length(levels) > 0) {
		for (i in 1:length(levels)) {
			cmb = combn(x = levels, i)
			for (j in 1:ncol(cmb)) {
				result = tapply(r, as.list(d[, cmb[, j], drop = FALSE]), function(x) {
					normality_test = if (!is.null(x) && !is.na(x) &&
															 length(x)           > 3  &&
															 length(unique(x)) > 3    &&
															 length(x)         < 5000 &&
															 var(x)            != 0) {
						shapiro.test(x)$p.value
					} else{
						'Not possible(Possible causes: <3 or >5000 unique data points'
					}
				})
				if (!is.null(result)) {
					flst[[counter]] = result
					names(flst)[counter] = paste(cmb[, j], collapse = collapse, sep = sep)
					counter = counter + 1
				}
			}
		}
		if (list) {
			flst = as.list(lapply(flst, function(f) {
				r = f
				if (length(dim(f)) > 1) {
					r = as.matrix(as.data.frame(f))
					r = unmatrix0(r)
				}
				r[is.na(r)] = noDataLab
				r = as.list(r)
				return(r)
			}))
		}
		return(flst)
	} else{
		return(NULL)
	}
}
