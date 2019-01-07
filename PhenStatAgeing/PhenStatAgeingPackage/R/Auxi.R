# Startup message
.onAttach <- function(lib, pkg) {
	packageStartupMessage(
		paste0(
			'\n >===============================================================================<',
			'\n PhenStatAgeing is developed by International Mouse Phenotyping Consortium (IMPC) ',
			'\n More details on https://www.mousephenotype.org/                                            ',
			'\n Contact us on hamedhm@ebi.ac.uk                                                  ',
			'\n >===============================================================================<'
		),
		domain = NULL,
		appendLF = TRUE
	)
}



# Test type (decimal [Continues]/not decimal [Categorical])
IsCategorical = function(num = NA, max.consider = Inf) {
	if (length(num) < 1 || is.null(num) || length(na.omit(num))<1) # do not change the order
		stop('\n => Null vector or all NAs \n')
	
	num = na.omit(num)
	if (!is.numeric(num)) {
		cat = TRUE
		method = 'FE'
	} else if (all(num %% 1 == 0) && length(unique(num)) < max.consider) {
		cat = TRUE
		method = 'GLM'
	} else{
		cat = FALSE
		method = 'MM'
	}
	return(
		list(cat = cat, method = method)
	)
}

# Typical fixed effect
TypicalModel = function(depVariable,
												withWeight = TRUE,
												Sex        = TRUE,
												LifeStage  = TRUE,
												mandatory  = 'Genotype'
												) {
	fixed =  reformulate(termlabels  = mandatory,	response = depVariable)
	if (Sex)
		fixed = update(fixed,  ~ . * Sex)
	if (LifeStage)
		fixed = update(fixed,  ~ . * LifeStage)
	if (withWeight)
		fixed = update(fixed,  ~ . + Weight)
	return(fixed)
}

# check for missings
CheckMissing = function(data, formula) {
	org.data = data
	new.data = data[complete.cases(data[, all.vars(formula)]), ]
	missings = ifelse(all(dim(org.data) == dim(new.data)),	0, dim(org.data)[1] -
											dim(new.data)[1])
	if (missings)
		message('\n => The data contain ',
						missings,
						' missing(s). \n => Missing data removed. \n ')
	return(invisible(
		list(
			org.data = org.data,
			new.data = new.data,
			missings = missings
		)
	))
	
}

# Cohens effect size
eff.size = function(object,
										data = NULL,
										depVariable,
										effOfInd = 'Genotype') {
	if (all(is.null(data)))
		data = getData(object)
	if (length(levels(data[, effOfInd])) != 2)
		stop('\n => There must be two levels in ', paste(effOfInd,collapse = ' and '), ' \n')
	f    = reformulate(termlabels = effOfInd, depVariable)
	agr  = aggregate(f, data = data, FUN = mean)
	if (any(dim(agr) != 2))
		stop(
			'\n => No variation in ',
			paste(effOfInd, collapse = ' or '),
			'; or more than one variable in effOfInd',
			'\n'
		)
	MDiff = agr[1, 2] - agr[2, 2]
	NModel = update(object,  as.formula(paste('~', effOfInd)), data = data)
	r      = resid(NModel)
	# This can be the second approach
	# NullModel = update(object,  ~1, data = data)
	# ef.size   = anova(NullModel,NModel)$L.Ratio[2]
	return(abs(MDiff) / sd(r))
}

# Check whether there is any variation in the data for a certain variables
noVariation = function(data, f = '~Genotype') {
	xtb = xtabs(formula =  f,
							drop.unused.levels = TRUE,
							data = data)
	if (any(dim(xtb) < 2)) {
		return(TRUE)
	}	else{
		return(FALSE)
	}
}

# Flatting a table
ConvDf2Flat = function(dframe,
											 ch1 = '*',
											 ch2 = ':',
											 chend = ';') {
	out = apply(as.data.frame(dframe), 1, function(x) {
		if (length(x) > 2) {
			paste(paste(paste(x[1:(length(x) - 1)], collapse = ch1), trimws(x[length(x)]), sep  = ch2), collapse = chend)
		} else if (length(x) == 2) {
			paste(paste(x[1], ch2, x[2]), collapse = chend)
		} else{
			paste(paste(x), collapse =  chend)
		}
	})
	return(out)
}

#### List of all procedures
lop = function() {
	return(
		c(
			'TYPICAL',
			'ABR',
			'ACS',
			'ALZ',
			'BLK',
			'BWT',
			'CAL',
			'CBC',
			'CHL',
			'CSD',
			'DXA',
			'ECG',
			'ECH',
			'ELZ',
			'EVL',
			'EVM',
			'EVO',
			'EVP',
			'EYE',
			'FER',
			'GEL',
			'GEM',
			'GEO',
			'GEP',
			'GPL',
			'GPM',
			'GPO',
			'GRS',
			'HEM',
			'HIS',
			'HWT',
			'IMM',
			'INS',
			'IPG',
			'OFD',
			'PAT',
			'VIA',
			'XRY'
		)
	)
}
