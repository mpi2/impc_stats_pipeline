PhenListAgeing = function(object,
													DOE = 'Batch',
													DOB = 'date_of_birth',
													d.threshold = 16 * 7,
													debug       = TRUE) {
	#### Function to add age and LifeStage to the data
	#### Negative age will be removed
	object@datasetPL           = droplevels(object@datasetPL)
	age.in.day = as.Date(object@datasetPL[, DOE]) - as.Date(object@datasetPL[, DOB])
	
	object@datasetPL$Age       = as.numeric(age.in.day)
	object@datasetPL$LifeStage = ifelse    (age.in.day > d.threshold, 'Late', 'Early')
	object@datasetPL$LifeStage = as.factor (object@datasetPL$LifeStage)
	object@datasetPL           = object@datasetPL[object@datasetPL$Age > 0, ]
	
	LL = levels(object@datasetPL$LifeStage)
	LS = levels(object@datasetPL$Sex)
	LG = levels(object@datasetPL$Genotype)
	
	message0 ('Age range: ', paste0(range(age.in.day), collapse = '-'))
	
	if (length(LL) != 2 && debug) {
		message0('Ageing pipeline requires two levels in the LifeStage. Levels: ',
						 pasteComma(LL,replaceNull = FALSE),
						 '\nNormal pipeline would apply to this data')
		#return(NULL)
	}
	if (length(LS) != 2 && debug) {
		message0('There should be two levels in Sex. Levels: ', pasteComma(LS,replaceNull = FALSE))
		#return(NULL)
	}
	if ('Weight' %in% names(object@datasetPL) &&
			sum(is.na(object@datasetPL[, 'Weight'])) > 0 && debug)
		message0('There are ', sum(is.na(object@datasetPL[, 'Weight'])), ' NAs in body weights.')
	
	if (length(LG) < 2 && debug) {
		message0('Genotype must have two levels. Levels: ', pasteComma(LG,replaceNull = FALSE))
		return(NULL)
	}
	object        = unclass(object)
	class(object) = 'PhenListAgeing'
	return(object)
}


summary.PhenListAgeing = function(object,
																	vars = NULL,
																	...) {
	data = object@datasetPL
	if ('LifeStage' %in% names(object@datasetPL))
		cnames = c('Batch',
							 'Genotype',
							 'Sex',
							 'Age',
							 'LifeStage')
	else
		cnames = c('Batch',
							 'Genotype',
							 'Sex')
	ndata = data[, if (is.null(vars)) {
		cnames
	} else{
		vars
	}]
	r = describe(ndata, ...)
	
	return(r)
}


plot.PhenListAgeing = function(x,
															 vars = NULL,
															 ...) {
	data = x@datasetPL
	if ('LifeStage' %in% names(data))
		cnames = c('Batch',
							 'Genotype',
							 'Sex',
							 'Age',
							 'LifeStage')
	else
		cnames = c('Batch',
							 'Genotype',
							 'Sex')
	ndata = data[, if (is.null(vars)) {
		cnames
	} else{
		vars
	}]
	r = describe(ndata, ...)
	
	plot = plot(r, ...)
	return(plot)
}
