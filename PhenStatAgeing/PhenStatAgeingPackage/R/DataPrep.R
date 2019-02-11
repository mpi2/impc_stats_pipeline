PhenListAgeing = function(object,
													DOE = 'Batch',
													DOB = 'date_of_birth',
													d.threshold = 16 * 7) {
	#### Function to add age and LifeStage to the data
	#### Negative age will be removed
	object@datasetPL           = droplevels(object@datasetPL)
	age.in.day = as.Date(object@datasetPL[, DOE]) - as.Date(object@datasetPL[, DOB])
	
	object@datasetPL$Age       = as.numeric(age.in.day)
	object@datasetPL$LifeStage = ifelse    (age.in.day > d.threshold, 'Late', 'Early')
	object@datasetPL$LifeStage = as.factor (object@datasetPL$LifeStage)
	object@datasetPL           = object@datasetPL[object@datasetPL$Age > 0,]
	
	message ('~> Age range: ',paste0(range(age.in.day),collapse = '-'))
	
	if (length(levels(object@datasetPL$LifeStage)) < 2)
		message('~> There is only one level in the LifeStage')
	
	if (length(levels(object@datasetPL$Sex)) < 2)
		message('~> There is only one level in Sex')
	
	if('Weight' %in% names(object@datasetPL) && sum(is.na(object@datasetPL[,'Weight']))>0)
		message('~> There are ',sum(is.na(object@datasetPL[,'Weight'])),' NAs in body weights.')
		
	if (length(levels(object@datasetPL$Genotype)) < 2)
		stop('~> Genotype must have two levels')
	
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
