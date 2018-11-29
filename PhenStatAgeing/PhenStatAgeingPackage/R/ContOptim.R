# Continues OPT core
M.opt = function(fixed = NULL,
								 random = NULL,
								 data = NULL,
								 direction = 'both',
								 family    = poisson(link = "log"),
								 trace = FALSE,
								 categorical = FALSE ,
								 LifeStage,
								 weight,
								 #BoxCoxtransform = FALSE,
								 ...) {
	sta.time = Sys.time()
	VarHomo = NULL
	if (!categorical) {
		cat('\n => Linear Mixed Model : Model optimization in progress ... ')
		# if (LifeStage) {
		# 	weight  = varIdent(form =  ~ 1 |
		# 										 	LifeStage)
		# } else{
		# 	weight  = varIdent(form =  ~ 1 |
		# 										 	Genotype)
		# }
		I.Model = do.call(
			'lme',
			list(
				random    = random,
				fixed     = fixed,
				data      = data,
				na.action = na.omit,
				method    = 'ML',
				weights   = weight
			)
		)
		F.Model = MASS::stepAIC(
			I.Model,
			trace = trace,
			direction = direction,
			scope = list(lower = ~ Genotype + 1),
			na.action = na.omit
		)
		if (!is.null(weight)) {
			F2.Model = do.call(
				"lme",
				list(
					random    = random,
					fixed     = formula(F.Model),
					data      = data,
					na.action = na.omit,
					method    = 'ML',
					weights   = NULL
				)
			)
			if (AIC(F2.Model) <= AIC(F.Model)) {
				# = is important
				F.Model = F2.Model
				VarHomo = FALSE
			} else{
				VarHomo = TRUE
			}
		}
		
		G.Model = do.call(
			"gls",
			args = list(
				model = formula(F.Model),
				data = data,
				na.action = na.omit,
				method = "ML",
				weights   = if (!is.null(VarHomo) && VarHomo)
					weight
				else
					NULL
			)
		)
		if (AIC(G.Model) <= AIC(F.Model)) {
			F.Model = G.Model
		}
	} else{
		cat('\n => Generalized Linear Model : Model optimization in progress ... ')
		I.Model = do.call('glm',
											list(
												formula   = fixed ,
												family    = family,
												data      = data
											))
		F.Model = MASS::stepAIC(
			I.Model,
			trace = trace,
			direction = direction,
			scope = list(lower = ~ Genotype + 1)
		)
	}
	cat('(', round(Sys.time() - sta.time, 4), 'seconds )\n')
	return(list(
		Initial.Model = I.Model,
		Final.Model = F.Model,
		VarHomo = VarHomo
	))
}


#### Summary core for continues data
summary.PhenListAgeingModel = function(object, ...) {
	if (is.null(object))
		stop('\n => NULL object\n')
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
	cat(
		'=> Final model : ',
		'\n\t Fixed  : ' ,
		ifelse(is.null(af), '-', paste(af, collapse = ', ')),
		'\n\t Random : ' ,
		ifelse(is.null(ar), '-', format(ar)),
		'\n\t ----- ' ,
		'\n\t Removed term(s) : ',
		ifelse(length(ai[!(ai %in% af)]) < 1, '-', paste(ai[!(ai %in% af)], collapse = ', ')),
		'\n\t ----- ' ,
		'\n\t Genotype effect size : ',
		ifelse(
			object$LifeStage,
			paste(
				c(
					'\n\t  Genotype for Early =>',
					'\n\t  Genotype for Late =>' ,
					'\n\t  LifeStage =>'
				),
				object$effect,
				collapse = ', '
			),
			object$effect
		),
		'\n\t ----- ' ,
		'\n\t Variance homogeneity:',
		ifelse(!is.null(object$VarHomo), object$VarHomo, 'Not specified'),
		'\n\t ----- ' ,
		'\n\n=> Model summary :\n'
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
plot.PhenListAgeingModel = function (x, ask = FALSE, ...) {
	if (is.null(x))
		stop('\n => NULL object\n')
	p = par()
	par(ask = ask)
	plot(
		predict(x$final.model),
		resid(x$final.model),
		xlab = 'Fitted values',
		ylab = 'Residuals',
		...
	)
	abline(h = 0)
	hist(
		resid(x$final.model),
		xlab = 'Residuals',
		main = 'Histogram of residuals',
		probability = TRUE,
		...
	)
	qqnorm(resid(x$final.model), main = 'Normal Q-Q plot of residuals', ...)
	plot(density(x$final.model$data[, x$depVariable]), main = 'Density of the response', ...)
	par(ask = p$ask)
}
