summary.PhenStatAgeingRR = function(object, format = 'rst', ...) {
	if (!is.null(object$messages) || is.null(object)) {
		message0('Due to error(s), no plot available')
		message0(object$messages)
		stop()
	}
	summaryCore(object, procedure = 'RR', format = format, ...)
}

summary.PhenStatAgeingFE = function(object, format = 'rst', ...) {
	if (!is.null(object$messages) || is.null(object)) {
		message0('Due to error(s), no plot available')
		message0(object$messages)
		stop()
	}
	summaryCore(object, procedure = 'FE', format = format, ...)
}

summary.PhenStatAgeingMM = function(object, format = 'rst', ...) {
	if (!is.null(object$messages) || is.null(object)) {
		message0('Due to error(s), no plot available')
		message0(object$messages)
		stop()
	}
	summaryCore(object, procedure = 'MM', format = format, ...)
}


summaryCore = function(x,
											 procedure = 'MM',
											 format = 'rst',
											 ...) {
	requireNamespace("knitr")
	vo = vectorOutputAgeing(object = x,
													JSON = FALSE,
													Null = FALSE)
	out = list(
		'Method'                         = vo$Method,
		'Model'                          = vo$`Additional information`$Formula$input,
		'----------------------------'   = '----------------------------',
		'Tested Gene'                    = vo$`Gp2 genotype`,
		'Reference Gene'                 = vo$`Gp1 genotype`,
		'----------------------------'   = '----------------------------',
		'Sex in the optimised model?'    = vo$`Genotype contribution`$`Sexual dimorphism detected`,
		'----------------------------'   = ifelse(
			procedure == 'RR',
			'- Separate p-values for (Low vs NormalHigh) and (LowNormal vs High) -',
			'----------------------------'
		),
		'Genotype contribution overal'   = vo$`Genotype p-val`,
		'Genotype contribution Females'  = vo$`Sex FvKO p-val`,
		'Genotype contribution Males'    = vo$`Sex MvKO p-val`,
		'----------------------------'   = '----------------------------',
		'LifeStage contribution'         = vo$`LifeStage p-val`,
		'Genotype contribution Early'    = vo$`LifeStage EvKO p-val`,
		'Genotype contribution Late'     = vo$`LifeStage LvKO p-val`,
		'----------------------------'   = '----------------------------',
		'Sex contribution'               = vo$`Sex p-val`,
		'Body weight contribution'       = vo$`Weight p-val`
	)
	outT = prepareSummaryOutput(out)
	print(kable(
		outT,
		format = format,
		col.names = c('Statistic', 'Value'),
		...
	))
	return(invisible(outT))
}

prepareSummaryOutput = function(out, nullMessage = 'Not applicable') {
	outT = as.matrix(out)
	outT = as.matrix(cbind(rownames(outT), outT))
	rownames(outT) = NULL
	return(outT)
}
