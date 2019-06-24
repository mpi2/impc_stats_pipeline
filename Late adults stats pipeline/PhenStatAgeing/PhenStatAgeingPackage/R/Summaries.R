summary.NULL = function(object, ...) {
	message0('No summary available for a NULL object')
}
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
	vo = vectorOutputAgeing(object           = x    ,
													JSON             = FALSE,
													ReportNullSchema = FALSE,
													RemoveNullKeys   = FALSE )
	pasteComma2 = function(...) {
		requireNamespace("rlist")
		inp = as.list(...)
		inp = list.clean(inp)
		if (is.null(inp)           ||
				length (inp) < 1)
			r = ' - '
		else
			r = pasteCommaJustForSummary(
				...,
				replaceNull   = TRUE  ,
				truncate      = FALSE ,
				replaceNullby = ' - '
			)
		return(r)
	}
	out = list(
		'Method'                         = vo$Method,
		'Final model'                    = if (procedure %in% 'MM') {
			if (!is.null(x$output$Final.Model))
				formula(x$output$Final.Model)
			else
				NULL
		} else if (procedure %in% 'FE') {
			RightFormula2LeftFormula(x$extra$Cleanedformula)
		}	else if (procedure %in% 'RR') {
			RightFormula2LeftFormula(x$extra$Cleanedformula)
		} else{
			NULL
		},
		'............................'   = '............................',
		'Tested Gene'                    = vo$`Gp2 genotype`,
		'Reference Gene'                 = vo$`Gp1 genotype`,
		'............................'   = '............................',
		'Sexual dimorphism detected?'    = vo$`Genotype contribution`$`Sexual dimorphism detected`,
		'............................'   = ifelse(
			procedure == 'RR',
			'* Separate p-values for (Low vs NormalHigh), (LowNormal vs High) and (details) ',
			'............................'
		),
		'Genotype contribution overal'   = pasteComma2(vo$`Genotype p-val`),
		'Genotype contribution Females'  = pasteComma2(vo$`Sex FvKO p-val`),
		'Genotype contribution Males'    = pasteComma2(vo$`Sex MvKO p-val`),
		'............................'   = '............................',
		'LifeStage contribution'         = pasteComma2(vo$`LifeStage p-val`),
		'Genotype contribution Early'    = pasteComma2(vo$`LifeStage EvKO p-val`),
		'Genotype contribution Late'     = pasteComma2(vo$`LifeStage LvKO p-val`),
		'............................'   = '............................',
		'Sex contribution'               = pasteComma2(vo$`Sex p-val`),
		'Body weight contribution'       = pasteComma2(vo$`Weight p-val`)
	)
	outT = prepareSummaryOutput(out)
	print(kable(
		outT,
		format = format,
		col.names = c('Statistic', 'Value'),
		...
	))
	outTRemoved = outT[!apply(outT, 1, function(x) {
		any(grepl(
			pattern = '..........',
			x = x,
			fixed = TRUE
		))
	}), ]
	return(invisible(outTRemoved))
}

prepareSummaryOutput = function(out, nullMessage = 'Not applicable') {
	outT = as.matrix(out)
	outT = as.matrix(cbind(rownames(outT), outT))
	rownames(outT) = NULL
	return(outT)
}

pasteCommaJustForSummary = function(...,
											replaceNull   = TRUE   ,
											truncate      = TRUE   ,
											width         = 100    ,
											trailingSpace = TRUE   ,
											replaceNullby = 'NULL') {
	sep = ifelse(trailingSpace, ', ', ',')
	if (replaceNull)
		r = paste(
			replaceNull(as.list(...), replaceBy = replaceNullby),
			sep      = sep,
			collapse = sep
		)
	else
		r = paste(...,
							sep      = sep,
							collapse = sep)
	
	if (truncate)
		r = truncate_text(r, width)
	
	return(r)
}
