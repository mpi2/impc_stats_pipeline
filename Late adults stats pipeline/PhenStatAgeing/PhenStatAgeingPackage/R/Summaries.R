summary.PhenStatAgeingRR = function(object, ...) {
	if (!is.null(object$messages) || is.null(object)) {
		message0('Due to error(s), no plot available')
		message0(object$messages)
		stop()
	}
	summaryCore(object, ...)
}

summary.PhenStatAgeingFE = function(object, ...) {
	if (!is.null(object$messages) || is.null(object)) {
		message0('Due to error(s), no plot available')
		message0(object$messages)
		stop()
	}
	summaryCore(object, ...)
}

summary.PhenStatAgeingMM = function(object, ...) {
	if (!is.null(object$messages) || is.null(object)) {
		message0('Due to error(s), no plot available')
		message0(object$messages)
		stop()
	}
	summaryCore(object, ...)
}


summaryCore = function(x, ...) {
	vo = vectorOutputAgeing(object = x,
													JSON = FALSE,
													Null = FALSE)
	out = list(
		'Method'                        = vo$Method,
		'Formula'                       = vo$`Additional information`$Formula$input,
		#'Dependent variable'           = vo$`Dependent variable`,
		'Tested Gene'                   = vo$`Gp2 genotype`,
		'Reference Gene'                = vo$`Gp1 genotype`,
		'Genotype contribution overal'  = vo$`Genotype contribution`$Overal,
		'Genotype contribution Females' = vo$`Genotype contribution`$`Sex FvKO p-val`,
		'Genotype contribution Males'   = vo$`Genotype contribution`$`Sex MvKO p-val`,
		'Sexual dimorphism detected?'    = vo$`Genotype contribution`$`Sexual dimorphism detected`,
		'Sex pvalue'                    = vo$`Sex p-val`,
		'Body weight p-value'            = vo$`Weight p-val`,
		'LifeStage p-value'              = vo$`LifeStage p-val`
	)
	outT = as.matrix(out)
	outT = as.matrix(cbind(rownames(outT), outT))
	rownames(outT) = NULL
	print(kable(
		outT,
		format = 'rst',
		col.names = c('Statistics', 'Value'),
		...
	))
	return(invisible(outT))
}

