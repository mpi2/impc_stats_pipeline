PhenListAgeing = function(PhenListobject,
													DOE = NULL,
													DOB = NULL,
													d.threshold = 16 * 7,
													debug       = TRUE) {
	#### Function to add age and LifeStage to the data
	#### Negative age will be removed
	PhenListobject@datasetPL           = droplevels(PhenListobject@datasetPL)
	
	Ageing = !is.null(DOE) &&
		       !is.null(DOB) &&
		       all(c(DOE, DOB) %in% names(PhenListobject@datasetPL))
	
	if (Ageing) {
		age.in.day = as.Date(PhenListobject@datasetPL[, DOE]) - as.Date(PhenListobject@datasetPL[, DOB])
		PhenListobject@datasetPL$Age       = as.numeric(age.in.day)
		PhenListobject@datasetPL$LifeStage = ifelse    (age.in.day > d.threshold, 'Late', 'Early')
		PhenListobject@datasetPL$LifeStage = as.factor (PhenListobject@datasetPL$LifeStage)
		PhenListobject@datasetPL           = PhenListobject@datasetPL[PhenListobject@datasetPL$Age > 0, ]
		message0 ('Age range: ', paste0(range(age.in.day), collapse = '-'))
	}
	
	LL = levels(PhenListobject@datasetPL$LifeStage)
	LS = levels(PhenListobject@datasetPL$Sex)
	LG = levels(PhenListobject@datasetPL$Genotype)
	
	if (length(LL) != 2 && debug && Ageing) {
		message0(
			'Ageing pipeline requires two levels in the LifeStage. Levels: ',
			pasteComma(LL, replaceNull = FALSE),
			'\nNormal pipeline will apply to this data'
		)
		#return(NULL)
	}
	if (length(LS) != 2 && debug) {
		message0('There should be two levels in Sex. Levels: ',
						 pasteComma(LS, replaceNull = FALSE))
		#return(NULL)
	}
	if ('Weight' %in% names(PhenListobject@datasetPL) &&
			sum(is.na(PhenListobject@datasetPL[, 'Weight'])) > 0 && debug)
		message0('There are ', sum(is.na(PhenListobject@datasetPL[, 'Weight'])), ' NAs in body weights.')
	
	if (length(LG) < 2 && debug) {
		message0('Genotype must have two levels. Levels: ',
						 pasteComma(LG, replaceNull = FALSE))
		return(NULL)
	}
	PhenListobject        = unclass(PhenListobject)
	class(PhenListobject) = 'PhenListAgeing'
	return(PhenListobject)
}


summary.PhenListAgeing = function(PhenListobject,
																	vars = NULL,
																	...) {
	data = PhenListobject@datasetPL
	if ('LifeStage' %in% names(PhenListobject@datasetPL))
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
