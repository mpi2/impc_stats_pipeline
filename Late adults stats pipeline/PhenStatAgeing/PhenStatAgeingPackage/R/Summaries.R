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
		'Sexual dimorphism detected?'    = vo$`Genotype contribution`$`Sexual dimorphism detected`,
		'----------------------------'   = ifelse(
			procedure == 'RR',
			'~ Separate p-values for (Low vs NormalHigh) and (LowNormal vs High) ~',
			'----------------------------'
		),
		'Genotype contribution overal'   = vo$`Genotype contribution`$Overal,
		'Genotype contribution Females'  = vo$`Genotype contribution`$`Sex FvKO p-val`,
		'Genotype contribution Males'    = vo$`Genotype contribution`$`Sex MvKO p-val`,
		'----------------------------'   = '----------------------------',
		'LifeStage contribution'              = vo$`LifeStage p-val`,
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
		col.names = c('Statistics', 'Value'),
		...
	))
	return(invisible(outT))
}


prepareSummaryOutput = function(out, nullMessage = 'Not applicable') {
	# out = lapply(out, function(x) {
	# 	if (length(x) < 1 ||
	# 			is.null(x) ||
	# 			is.character(x)) {
	# 		if (is.character(x) && nchar(x) > 3) {
	# 			x = gsub(
	# 				x = x,
	# 				pattern = 'NULL',
	# 				replacement =  nullMessage,
	# 				fixed = TRUE
	# 			)
	# 		}
	# 	} else if (is.null(x)) {
	# 		x = nullMessage
	# 	}else{
	# 		x = x
	# 	}
	# 	return(x)
	# })
	outT = as.matrix(out)
	outT = as.matrix(cbind(rownames(outT), outT))
	rownames(outT) = NULL
	return(outT)
}
