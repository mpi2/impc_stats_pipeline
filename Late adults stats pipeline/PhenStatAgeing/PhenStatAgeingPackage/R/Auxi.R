# Startup message
.onAttach <- function(lib, pkg) {
	packageStartupMessage(
		paste0(
			'\n >===============================================================================<',
			'\n PhenStatAgeing is developed by International Mouse Phenotyping Consortium (IMPC) ',
			'\n More details           : https://www.mousephenotype.org/                         ',
			'\n Source code and issues : https://git.io/fj9W3                                    ',
			'\n Contact us             : hamedhm@ebi.ac.uk                                       ',
			'\n >===============================================================================<'
		),
		domain = NULL,
		appendLF = TRUE
	)
}

Matrix2List = function(x, ...) {
	if (is.null(x))
		return(NULL)
	
	if (length(x) == 1 || class(x) == 'numeric') {
		return(as.list(x))
	}
	
	if (!is(x, 'matrix'))
		x = as.matrix(x)
	r = as.list(unmatrix0(x, ...))
	return(r)
}

replaceNull = function(x, replaceBy = 'NULL') {
	x = sapply(x, function(xx) {
		if (is.null(xx)) {
			replaceBy
		} else{
			xx
		}
	})
	return(unlist(x))
}
pasteComma = function(...,
											replaceNull   = TRUE   ,
											truncate      = TRUE   ,
											width         = 100    ,
											trailingSpace = TRUE   ,
											replaceNullby = 'NULL',
											sep = ',') {
	sep = ifelse(trailingSpace && sep %in% ',', paste0(sep, ' '), sep)
	if (replaceNull)
		r = paste(
			replaceNull(list(...), replaceBy = replaceNullby),
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

truncate_text = function(x, width) {
	ifelse(nchar(x) > width, paste0(strtrim(x, width), ' ...'), x)
}

pasteUnderscore = function(...) {
	paste(..., sep = '_', collapse = '_')
}
pastedot = function(...) {
	paste(..., sep = '.', collapse = '.')
}
checkModelTermsInData = function(formula,
																 data,
																 responseIsTheFirst  = TRUE,
																 pattern = '[.~+-()]') {
	formula = as.formula(formula)
	vars    = all_vars0(formula, functions = FALSE)
	vars    = vars[!grepl(pattern = pattern,
												x = vars,
												fixed = TRUE)]
	if (responseIsTheFirst) {
		if (!(vars[1] %in% names(data))) {
			message0('Response does not exist in the data!\n\tFormula: ',
							 printformula(formula))
			stop(
				'Response has not been specified properly. Please check that the response exists in the data'
			)
		}
	}
	In      = vars %in% names(data)
	if (any(!In)) {
		message0(
			'Some terms in the model are not included in the data. See: \n\t  ',
			pasteComma(vars[!In], replaceNull = FALSE, truncate = FALSE),
			'\n\t Initial  model: ',
			printformula(formula)
		)
		ft  = vars [!In]
		formula = update.formula(formula,
														 reformulate0(
														 	termlabels = c('.', ft, paste0(ft, ':.')),
														 	response   = NULL,
														 	intercept  = TRUE,
														 	sep        = '-'
														 ))
		message0('\t Polished model: ', printformula(formula))
	}
	return(formula)
}


reformulate0 = function (termlabels,
												 response = NULL,
												 intercept = TRUE,
												 sep = '+')
{
	if (!is.character(termlabels) || !length(termlabels))
		stop("'termlabels' must be a character vector of length at least one")
	has.resp <- !is.null(response)
	termtext <- paste(if (has.resp)
		"response",
		"~",
		paste(termlabels, collapse = sep),
		collapse = "")
	if (!intercept)
		termtext <- paste(termtext, "- 1")
	rval <- eval(parse(text = termtext, keep.source = FALSE)[[1L]])
	if (has.resp)
		rval[[2L]] <- if (is.character(response))
			as.symbol(response)
	else
		response
	environment(rval) <- parent.frame()
	rval
}

suppressMessagesANDWarnings = function(exp,
																			 sup.messages = TRUE,
																			 sup.warnings = FALSE) {
	if (sup.messages && sup.warnings) {
		suppressMessages(suppressWarnings(exp))
	} else if (sup.messages && !sup.warnings) {
		suppressMessages(exp)
	} else if (!sup.messages && sup.warnings) {
		suppressWarnings(exp)
	} else{
		exp
	}
}
# Typical fixed effect
TypicalModel = function(depVariable,
												withWeight = TRUE,
												Sex        = TRUE,
												LifeStage  = TRUE,
												mandatory  = 'Genotype',
												data                   ,
												others     = NULL      ,
												debug      = TRUE) {
	colNames = colnames(data)
	if (!mandatory %in% colNames)
		stop('Genotype does not found in the dataset!')
	
	fixed =  reformulate(termlabels  = mandatory,	response = depVariable)
	if (Sex && ('Sex' %in% colNames))
		fixed = update(fixed,  ~ . * Sex)
	if (LifeStage && ('LifeStage' %in% colNames))
		fixed = update(fixed,  ~ . * LifeStage)
	if (withWeight && ('Weight' %in% colNames))
		fixed = update(fixed,  ~ . + Weight)
	if (!is.null(others))
		fixed = update(fixed,  reformulate(response = NULL,	termlabels = c('.', others)))
	####
	if (debug)
		message0('Initial  model: ', printformula(fixed))
	return(fixed)
}



FeasibleTermsInContFormula = function(formula, data) {
	if (is.null(formula) || is.null(data)) {
		message0('Null data or the formula. Check the data and/or formula')
		stop()
	}
	Allvars = all_vars0(formula)[all_vars0(formula) %in% names(data)]
	isCat = !sapply(data[, Allvars, drop = FALSE], is.numeric)
	vars  = Allvars[isCat]
	lvars = length(vars)#min(length(vars), sapply(strsplit(formulaTerms(formula = formula), split = ':'), length), na.rm = TRUE)
	names = r = NULL
	if (getResponseFromFormula(formula = formula) %in% vars) {
		message0('\tResponse is included in the checks ....')
	}
	if (lvars > 0) {
		for (i in 1:lvars) {
			message0('\t',
							 i,
							 ' of ',
							 lvars,
							 '. Checking for the feasibility of terms and interactions ...')
			cmb = combn(vars, i)
			for (j in 1:ncol(cmb)) {
				message0('\t\t Checking ', pasteComma(cmb[, j]))
				xtb = xtabs(
					formula = paste0('~', paste0(cmb[, j], collapse = '+')),
					data = data,
					drop.unused.levels = FALSE
				)
				r = c(r, if (all(dim(xtb) >= 2)) {
					min(xtb, na.rm = TRUE)
				} else{
					0
				})
				names = c(names, paste0(cmb[, j], collapse = ':'))
			}
		}
		return(data.frame(
			names = names,
			min.freq = r,
			stringsAsFactors = FALSE
		))
	} else{
		return(NULL)
	}
}

variablesInData = function(df, names, debug = TRUE) {
	if (is.null(df) || is.null(names) || sum(names %in%  names(df)) < 1)
		return(NULL)
	newNames = names[names %in%  names(df)]
	if (debug)
		message0('Variables that being found in data: ',
						 pasteComma(newNames, truncate = FALSE))
	return(newNames)
}

ComplementaryFeasibleTermsInContFormula = function(formula, data) {
	message0(
		'Checking for the feasibility of terms and interactions ...\n\t Formula: ',
		printformula(formula)
	)
	fbm = FeasibleTermsInContFormula(formula = formula, data = data)
	if (!is.null(fbm)          &&
			(min(fbm$min.freq, na.rm = TRUE) < 1 ||
			 length(formulaTerms(formula)) != nrow(fbm))) {
		formula = update.formula(
			old = formula,
			new =
				reformulate0(
					termlabels = c('.', fbm$names[fbm$min.freq <= 0]),
					response = '.',
					intercept = TRUE,
					sep = '-'
				)
		)
		if (min(fbm$min.freq, na.rm = TRUE) < 1)
			message0(
				'The following term(s) removed because there is either "no data" or "no data for the interactions":\n\t ** Note. Not all terms necessarily in the initial model \n\t ',
				pasteComma(fbm[fbm$min.freq <= 0, c('names')], replaceNull = FALSE, truncate = FALSE)
			)
	}
	return(formula)
}

sign0 = function(x) {
	if (is.null(x))
		return(NULL)
	if (sign(x) > 0)
		return('positive')
	else if (sign(x) == 0)
		return('neutral')
	else if (sign(x) < 0)
		return('negative')
	else
		return(NULL)
}

dist0 = function(x, func = lower.tri) {
	if (is.null(x))
		return(NULL)
	out = outer (x, x, `-`)
	r   = out[func(out)]
	return(r)
}

CheckMissing = function(data, formula) {
	if (is.null(formula) || is.null(data)) {
		message0('Null data or the formula. Check the data and/or formula')
		stop()
	}
	org.data = data
	new.data = data[complete.cases(data[, all_vars0(formula)]),]
	missings = ifelse(all(dim(org.data) == dim(new.data)),	0, dim(org.data)[1] -
											dim(new.data)[1])
	if (missings)
		message0(
			'The data (variable(s) = ',
			pasteComma(all_vars0(formula), truncate = FALSE),
			') contain ',
			missings,
			' missing(s) ...\n\tMissing data removed'
		)
	return(invisible(
		list(
			org.data = org.data,
			new.data = droplevels0(new.data),
			missings = missings
		)
	))
}

droplevels0 = function(x, ...) {
	if (is.null(x)                        ||
			class(x) %in% c('matrix', 'integer', 'double', 'numeric'))
		return(x)
	return(droplevels(x, ...))
}

range0 = function(x, ...) {
	ran = range(x, na.rm = TRUE)
	if (length(ran) == 2) {
		return(diff(ran))
	} else{
		return(NULL)
	}
}

order0 = function(x, levels = FALSE) {
	if (is.null(x))
		return(NULL)
	
	if (levels) {
		r = x[order(levels(x))]
	} else{
		r = x[order(x)]
	}
	return(r)
}



percentageChangeCont = function(model                ,
																data                 ,
																variable             ,
																depVar               ,
																individual = TRUE    ,
																mainEffsOnlyWhenIndivi = 'Sex',
																FUN = range0,
																sep = ' ') {
	if (!(
		!is.null(data)                                           &&
		(!is.null(variable) || !is.null(mainEffsOnlyWhenIndivi)) &&
		!is.null(model)                                          &&
		!is.null(FUN(data[, depVar]))                            &&
		FUN(data[, depVar]) != 0
	))
		return(NULL)
	####
	model       =
		tryCatch(
			expr = update(model                          ,
										data = NormaliseDataFrame(data)),
			error = function(e) {
				message0(
					'\t\tError(s) in the (combined) effect size estimation for',
					pasteComma(variable),
					'. See: '
				)
				message0('\t\t', e, breakLine = FALSE)
				return(NULL)
			} ,
			warning = function(w) {
				message0(
					'\t\tWarning(s) in the (combined) effect size estimation for',
					pasteComma(variable),
					'. See: '
				)
				message0('\t\t', w, breakLine = FALSE)
				return(NULL)
			}
		)
	if (is.null(model))
		return(NULL)
	coefs = unlist(
		suppressMessagesANDWarnings(
			summary(model, verbose = FALSE)$tTable[, 1],
			sup.messages = TRUE,
			sup.warnings = TRUE
		)
	)
	######
	ran   = FUN(data[, depVar])
	if (is.null(coefs) || is.null(ran) || is.nan(ran))
		return(NULL)
	######
	if (individual) {
		out   = sapply(variable, function(x) {
			message0('\tCalculating the percentage change for: ',
							 pasteComma(x))
			if (is.numeric(data[, x])) {
				r = coefs[-1]
				names(r) = NULL
			} else{
				r  = coefs
				lr = length(r)
				ol = order0(levels(data[, x]))
				##########
				if (lr != nlevels(data[, x])) {
					message0(
						'\tCare reguired for the percentage change calculation for: ',
						x,
						'. Levels = ',
						names(r),
						'Coefficients: ',
						r
					)
					names(r) = paste(ol[1:lr], '_CareRequired')
				} else{
					names(r) = ol
				}
			}
			return(r / ran * 100)
		}, USE.NAMES = TRUE)
		return(Matrix2List(out, sep = sep))
	} else{
		out = extractCoefOfInterest(coefs  = coefs,
																main   = mainEffsOnlyWhenIndivi,
																data   = data)
		return(
			list(
				'Value'               = as.list(out)                             ,
				'Variable'            = mainEffsOnlyWhenIndivi                   ,
				'Model'               = printformula(formula(model))             ,
				'Type'                = 'Standardized interaction coefficients ' ,
				'Percentage change' = Matrix2List(out / ran * 100, sep = sep)
			)
		)
	}
}

optimM = function(optimise) {
	paste0(c('Fixed term = ', 'Weight term = ', 'Random term = '),
				 optimise,
				 collapse = ', ')
}

applyFormulaToData = function(formula = NULL, data, add = FALSE) {
	if (is.null(formula) || is.null(all_vars0(formula)))
		return(data)
	if (is.null(data))
		return(NULL)
	nms =
		trimws(scan(
			text = paste(unlist(as.list(
				attr(terms(as.formula(formula)), "variables")
			))[-1], sep = '+', collapse = ' + '),
			what = "",
			sep = "+",
			quiet = TRUE
		))
	m <- sapply(nms, function(x)
		eval(parse(text = x), data))
	if (!is.null(m) && add)
		m = cbind(data, m)
	return(list(data = m, names = nms))
}

optimiseMessage = function(optimise) {
	message0('The model optimisation is ',
					 ifelse(
					 	all(optimise),
					 	'in progress ...',
					 	paste0('in the following order:\n\t', optimM(optimise))
					 ))
}

extractCoefOfInterest = function(coefs, main = 'Sex', data) {
	lvls = lapply(
		main,
		FUN = function(x) {
			if (!is.numeric(data[, x])) {
				r = paste(x, levels(data[, x]), sep = '')
				r = c(paste0(r, ':'), paste0(':', r))
			} else{
				r = x
			}
			return(r)
		}
	)
	
	Cnames = names(coefs)
	r = coefs[grepl(pattern = pasteBracket(
		unlist(lvls),
		left = '',
		right = '',
		col = '|'
	),
	x = Cnames)]
	return(r)
}

pasteBracket = function(...,
												replaceNull = TRUE,
												col = '|',
												right = ']:',
												left = '[') {
	if (replaceNull)
		paste(
			left,
			replaceNull(list(...), replaceBy = 'NULL'),
			right,
			sep = '',
			collapse = col
		)
	else
		paste(left, ..., right, sep = ' ', collapse = col)
}

eff.size = function(object,
										data        = NULL        ,
										depVariable = 'data_point',
										effOfInd    = 'Genotype'  ,
										errorReturn = NULL        ,
										debug       = FALSE) {
	if (all(is.null(data)))
		data = NormaliseDataFrame(getData(object))
	f    = reformulate(termlabels = effOfInd, depVariable)
	# Remove NAs
	data = data[complete.cases(data[, all_vars0(f)]), ]
	if (!(depVariable %in% names(data)          &&
				length(na.omit(data[, depVariable])) > 1))
		return(NULL)
	
	agr  = aggregate(f, data = data, FUN = function(x){mean(x,na.rm = TRUE)})
	if (debug) {
		cat('\n\n')
		print(agr)
		print(dim(agr))
	}
	if (any(dim(agr) < 2)) {
		message0(
			' \t\t(standardized) Effect size estimation: No variation or less than two levels in ',
			pasteComma(effOfInd, collapse = ',')
		)
		return(errorReturn)
	}
	
	NModel       =
		tryCatch(
			expr = update(
				object,
				reformulate(
					termlabels = effOfInd,
					response = '.',
					intercept = TRUE
				),
				data = data
			),
			error = function(e) {
				message0('\t\tError(s) in the effect size estimation for',
								 pasteComma(effOfInd),
								 '. See: ')
				message0('\t\t', e, breakLine = FALSE)
				return(NULL)
			} ,
			warning = function(w) {
				message0('\t\tWarning(s) in the effect size estimation for',
								 pasteComma(effOfInd),
								 '. See: ')
				message0('\t\t', w, breakLine = FALSE)
				return(NULL)
			}
		)
	if (!is.null(NModel)) {
		CoefEffSizes = sapply(data[, effOfInd, drop = FALSE], FUN = is.numeric)
		PerChange = percentageChangeCont(
			model    = NModel    ,
			data     = data      ,
			variable = effOfInd  ,
			depVar   = depVariable
		)
		if (sum(CoefEffSizes)) {
			# For continues covariates it is the coefficient
			efSi         = list(
				'Value'               = as.list(coef(NModel))[[effOfInd]][1],
				'Variable'            = effOfInd                            ,
				'Model'               = printformula(formula(NModel))       ,
				'Type'                = 'Standardized coefficient'          ,
				'Percentage change' = PerChange
			)
		} else{
			# For categorical covariates it is the mean difference
			MDiff        = max(dist(agr[, depVariable, drop = FALSE]), na.rm = TRUE)
			r            = resid(NModel)
			sd           = sd0(r, na.rm = TRUE)
			efSi         = list(
				'Value'               = ifelse(!is.na(sd) &&
																			 	sd > 0, abs(MDiff) / sd, NA),
				'Variable'            = effOfInd                            ,
				'Model'               = printformula(formula(NModel))       ,
				'Type'                = 'Mean difference'                   ,
				'Percentage change' = PerChange
			)
		}
	} else{
		efSi = errorReturn
	}
	return(efSi)
}

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

diff0 = function(x) {
	r = if (length(x) < 2) {
		x
	}	else{
		diff(x)
	}
	return(r)
}

summary0 = function(x, ...) {
	if (is.null(x) || length(x) < 1) {
		message0('Null column found in the data!')
		return(x)
	}
	if (is.numeric(x)) {
		r  = summary(diff0(x), ...)
	} else{
		xx = as.factor(x)
		r  = c(sort(levels(xx)), summary(diff0(as.integer(xx)), ...))
	}
	return(r)
}

RemoveDuplicatedColumnsFromDf = function(x, formula = NULL) {
	x         = as.data.frame(x)
	if (!is.null(formula)                      ||
			sum(all_vars0(formula) %in% names(x)) > 1)
		vars    = all_vars0(formula)
	else
		vars    = names(x)
	colVars   = names(x)  %in% vars
	if (sum(colVars)) {
		message0('Checking duplications in the data model:\n\t ',
						 pasteComma(names(x)[colVars],
						 					 truncate = FALSE))
		subX    = x[, colVars, drop = FALSE]
		#numCols = sapply(subX, is.numeric)
		#ConCols = subX[,  numCols, drop = FALSE]
		#CatCols = subX[, !numCols, drop = FALSE]
		dcols   = duplicated(lapply(subX, summary0))
		if (any(dcols)) {
			message0(
				'\tDuplicated columns found (and removed) in the input data. Removed variables:\n\t ',
				pasteComma(names(subX)[dcols], truncate = FALSE)
			)
		} else{
			message0('\tNo duplicate found')
		}
		uniqCols  = subX   [, !dcols  , drop = FALSE]
		r         = cbind(x[, !colVars, drop = FALSE], uniqCols)
		return(r)
	} else{
		message0('Formula terms do not exist in the input data')
		return(x)
	}
}

colExists = function(name, data) {
	if ((name %in% names(data)) &&
			length(complete.cases(data[, name])) > 0) {
		return(TRUE)
	} else{
		return(FALSE)
	}
}

listFun = function(list, FUN, debug = FALSE) {
	if (debug)
		message0('\tApplied model: ', FUN)
	fArgs = names(list) %in% formalArgs(args(FUN))
	l     = list[names(list)[fArgs]]
	return(l)
}
RandomEffectCheck = function(formula, data) {
	if (is.null(formula) || is.null(data))
		return(NULL)
	message0('Checking the random effect term ...\n\tFormula: ',
					 printformula(formula))
	difTerms = setdiff(all_vars0(formula), names(data))
	if (length(difTerms)) {
		message0(
			'\tSome terms in the random effect do not exist in the data. See:\n\t ',
			difTerms,
			'\n\tRandom effect is set to NULL'
		)
		return(NULL)
	} else{
		return(formula)
	}
}
ModelChecks = function(fixed,
											 data,
											 checks = c(0, 0, 0),
											 responseIsTheFirst = TRUE) {
	if (length(checks) != 3) {
		message0('"checks" must be a vector of 3 values. Example c(1,1,1) or c(1,1,0) or c(0,0,0)')
		return(fixed)
	}
	if (any(checks > 0)) {
		if (checks[1])
			fixed  = checkModelTermsInData(formula = fixed,
																		 data = data,
																		 responseIsTheFirst = TRUE)
		if (checks[2])
			fixed = removeSingleLevelFactors(formula = fixed, data = data)
		if (checks[3])
			fixed = ComplementaryFeasibleTermsInContFormula(formula = fixed, data = data)
		message0('\tChecked model: ', printformula(fixed))
	}
	return(fixed)
}

termInTheModel = function(model, term, message = FALSE) {
	if (message)
		message0(
			'Search for Term:\n Term: ',
			paste(term, collapse = ',') ,
			'\n\t Model: ' ,
			paste(formula(model), collapse = ', ')
		)
	return(all(term %in% all_vars0(formula(model))))
}

SplitEffect = 	function(finalformula     ,
												fullModelFormula  ,
												F.Model           ,
												data              ,
												depVariable       ,
												mandatoryVar = 'Genotype',
												ci_levels    = .95) {
	Allargs  = all_vars0(fullModelFormula)[!all_vars0(fullModelFormula) %in% c(depVariable, mandatoryVar)]
	if (is.null(Allargs)) {
		message0 ('Nothing to split on ...')
		return(NULL)
	}
	isCat    = !sapply(data[, Allargs, drop = FALSE], is.numeric)
	args     = Allargs[isCat]
	argsCon  = if (length(Allargs[!isCat]) > 0) {
		Allargs[!isCat]
	} else{
		#	NULL
		1
	}
	largs    = length(args)
	dname    = names(data)
	l        = NULL
	names    = c()
	counter = 1
	if (largs > 0) {
		for (i in 1:largs) {
			argComb = combn(args, i, simplify = TRUE)
			for (j in 1:ncol(argComb)) {
				arg = as.vector(argComb[, j])
				if (arg %in% dname &&
						!is.numeric(data[, arg]) &&
						# model = finalformula
						termInTheModel(model = fullModelFormula, term = arg)) {
					message0(counter,
									 '. Split on ',
									 paste0(arg, collapse = ' & '),
									 ' ...')
					newModel = reformulate(
						response = depVariable,
						termlabels = c(
							arg              ,
							paste(
								mandatoryVar   ,
								arg            ,
								collapse = ':' ,
								sep = ':'
							),
							argsCon
						),
						intercept = TRUE
					)
					message0('Checking the split model:\n\t',
									 printformula(newModel))
					
					l0    = tryCatch(
						update(F.Model,
									 newModel),
						error = function(e) {
							message0(e, breakLine = FALSE)
							return(NULL)
						} ,
						warning = function(w) {
							message0(w, breakLine = FALSE)
							return(NULL)
						}
					)
					message0(
						'\tTested model: ',
						printformula(newModel),
						ifelse(!is.null(l0), ' [Successful]', ' [Failed]'),
						breakLine = FALSE
					)
					if (!is.null(l0)) {
						l0  = intervalsCon (object = l0, lvls = ci_levels)
						l0$MainEffect   = arg
						l0$SplitFormula = printformula(newModel)
						l[[counter]]    = l0
						names[counter]  = paste(
							paste0(mandatoryVar, collapse = '_'),
							paste0(arg         , collapse = '.'),
							collapse = '_',
							sep      = '_'
						)
						counter = counter + 1
					}
				}
			}
		}
		if (length(names) > 0) {
			message0('SplitEffects. Output names: ',
							 paste0(names,  collapse = ', '))
			names(l) = paste0(names)
		}
	}
	return(l)
}

dim0 <- function(...) {
	args <- list(...)
	sapply(args , function(x) {
		if (is.null(dim(x)))
			return(length(x))
		dim(x)
	})
}

printformula = function(formula, message = TRUE) {
	if (!is.null(formula)) {
		r = paste01(format(formula, trim = TRUE, width = 0), collapse = '')
	} else{
		if (message)
			message0('Ops! the formula is blank')
		r = NULL
	}
	return(r)
}
# Categorical effect size
# https://github.com/mpi2/stats_working_group/raw/master/PhenStatUserGuide/PhenStatUsersGuide.pdf p122
cat.eff.size = function(xtb,
												varName = NULL,
												formula = NULL) {
	if (any(dim(xtb) < 1)) {
		r = NULL
	} else{
		r = max(apply(prop.table(xtb, margin = 2), 1, function(x) {
			max(dist(x, method = 'maximum', diag = TRUE), na.rm = TRUE)
		}), na.rm = TRUE)
	}
	out = list(
		value               = r                                               ,
		variable            = ifelse(!is.null(varName)                        ,
																 varName                                  ,
																 'Variable does not exist')               ,
		model               = printformula(RightFormula2LeftFormula(formula, removeResIfExists = TRUE)),
		type                = 'Proportion change'                             ,
		'percentage change' = NULL
	)
	return(out)
}

printVarFreqfromTable = function(tbl) {
	if (is.null(tbl) || any(dim(tbl) < 1))
		return(NULL)
	r = pasteComma(paste0(names(tbl), '[', tbl, ']'))
	return(r)
}

GenderIncludedInAnalysis = function(x, sexCol = 'Sex') {
	r = if (!sexCol %in% colnames(x)) {
		paste0(sexCol, ' does not included in the input data')
	}	else if (nlevels(x[, sexCol]) > 1) {
		paste0('Both sexes included; ', printVarFreqfromTable(table(x$Sex)))
	} else{
		paste0('Only one sex included in the analysis; ',
					 printVarFreqfromTable(table(x$Sex)))
	}
	return(r)
}

# Test engine
ctest = function(x,
								 formula   = NULL  ,
								 asset     = NULL  ,
								 rep       = 1500  ,
								 ci_levels = 0.95  ,
								 RRextraResults   = NULL             ,
								 overalTableName  = 'Complete table' ,
								 ...) {
	xtb = xtabs(
		formula = formula        ,
		data = x                 ,
		drop.unused.levels = TRUE,
		na.action = 'na.omit'
	)
	
	checkRes = checkTableForFisherTest(xtb = xtb, asset = asset)
	if (!checkRes$passed) {
		r        = setNames(list(overal = list(
			p.value = NULL, effect = NULL
		)), overalTableName)
	} else{
		if (any(dim(xtb) < 2)) {
			r        = setNames(list(overal = list(
				p.value = 1, effect = 1
			)), overalTableName)
		} else{
			r = fisher.test1(
				x           = xtb            ,
				formula     = formula        ,
				ci_levels   = ci_levels      ,
				simulate.p.value = rep > 0   ,
				conf.int    = TRUE           ,
				conf.level  = ci_levels      ,
				B           = rep            ,
				overalTableName = overalTableName ,
				...
			)
		}
	}
	return(
		list(
			result      = r               ,
			note        = checkRes$message,
			table       = xtb             ,
			input       = x               ,
			formula     = formula         ,
			RRextra     = RRextraResults
		)
	)
}

# check table for fisher.test
checkTableForFisherTest = function(xtb, asset = NULL) {
	message = NULL
	if (!is.null(asset)) {
		if (asset <= dim(xtb)[3]) {
			xtb = xtb[, , asset]
		} else{
			message = c(message, 'There are some empty levels in data')
			return(list(passed     = FALSE       ,
									note       = message))
		}
	}
	if (length(dim0(xtb)) < 2                      ||
			sum(margin.table(xtb, margin = 2) > 0) < 2 ||
			sum(margin.table(xtb, margin = 1) > 0) < 2) {
		message = c(message,
								'No variation in al least two levels of one or more categories')
		return(list(passed   = FALSE  ,
								note     = message))
	}
	return(list(passed   = TRUE  ,
							note     = message))
}

# Bulletproof fisher.test
fisher.test0 = function(x, formula, ci_levels, ...) {
	r = tryCatch(
		fisher.test(x,
								...),
		error = function(e) {
			message0(e, breakLine = FALSE)
			return(NULL)
		} ,
		warning = function(w) {
			message0(w, breakLine = FALSE)
			return(NULL)
		}
	)
	r           = intervalsCat(r, ci_levels)
	r$effect    = cat.eff.size(x,
														 varName = pasteComma(all_vars0(formula)[-1], truncate = FALSE),
														 formula = formula)
	r$formula   = formula
	r$table     = x
	return(r)
}

# Fisher test with broken table
fisher.test1 = function(x,
												formula,
												ci_levels,
												overalTableName = 'Complete table',
												...) {
	if (is.null(x))
		return(NULL)
	
	nrx = nrow(x)
	outList = NULL
	if (nrx > 2) {
		message0(
			'\t\t testing sub tables in progress ...\n\t\t\t Total tests: ',
			ncombn(n = nrx, x = 2:nrx),
			'; ',
			pasteComma(names(attributes(x)$dimnames))
		)
		for (i in 2:nrx) {
			cbn = combn(nrx, i)
			subtbl = lapply(1:ncol(cbn), function(j) {
				r = suppressMessages(fisher.test0(
					x = x[cbn[, j], ],
					formula = formula,
					ci_levels = ci_levels,
					...
				))
				if (nrow(cbn[, j, drop = FALSE]) == nrow(x)) {
					# this name is used in more places!
					r$data.name = overalTableName
				} else{
					r$data.name = paste(rownames(x[cbn[, j], ]), sep = '.', collapse = '.')
				}
				return(r)
			})
			
			outList = c(outList, setNames(object = subtbl, lapply(subtbl, function(x)
				x$data.name)))
		}
	} else{
		outList = setNames(list(overal = fisher.test0(
			x = x,
			formula = formula,
			ci_levels = ci_levels,
			...
		)), overalTableName)
	}
	return(outList)
}

ncombn = function(n,x,total = TRUE){
	r = lfactorial(n)- (lfactorial(x)+lfactorial(n-x))
	if(total)
		return(sum(exp(r)))
	else
		return(exp(r))
}

# Change with super extra care
# This is a complicated function for modeling all variation of the variables
AllTables = function(dframe        = NULL,
										 # dataframe
										 vars          = NULL,
										 # list of categorical variables
										 cl            = 0   ,
										 # lock the number of columns: works only if shrinke = TRUE
										 response.name = NULL,
										 # response name in the beginning of each variable name [only]
										 cols          = NULL,
										 # filter on the columns of interest : works only if shrinke = TRUE
										 shrink        = FALSE) {
	# remove no variation levels (columns)
	if (is.null(dframe)) {
		message0('Null data frame ')
		return(NULL)
	}
	
	cat  = vars
	lcat = length(cat)
	# Make all tables
	l2     = list()
	message0('\tSplitting in progress ...')
	for (i in unique(pmax(1, 1:(lcat - 1)))) {
		cb = combn(cat, i)
		for (j in 1:ncol(cb)) {
			message0('\tspliting on ', pasteComma(cb[, j], replaceNull = FALSE))
			out = split(dframe,
									interaction(dframe[cb[, j]]), drop = TRUE)
			l2   = c(l2, out)
		}
	}
	# Remove fixed value columns
	message0('\tShrinking in progress ...')
	if (shrink)
		l2 = lapply(l2, function(x) {
			# remove fixed value columns
			NonZeroFreq = c(apply(x, 2, function(xx) {
				length(unique(na.omit(xx))) > 1
			}))
			NonZeroFreq [names(NonZeroFreq)  %in% c('Freq', response.name)] = TRUE
			r = x[, NonZeroFreq , drop = FALSE]
			return(r)
		})
	
	message0('\tKeeping variables of interest ...')
	if (!is.null(cols))
		l2 = l2[as.logical(lapply(
			# keep certain columns
			l2,
			FUN = function(x) {
				all(cols %in% colnames(x))
			}
		))]
	
	# You do not want to get split on all combinations!
	if (all(cl > 0))
		# Fixed number of columns in each element of the list
		l2 = l2[which(lapply(
			l2,
			FUN = function(x) {
				ncol(x) %in% cl
			}
		) == TRUE)]
	
	# Split the big tables into 2x2 tables
	message0('\tDichotomising the final tables ...')
	oblCols = c(response.name, 'Freq')
	l3 = lapply(
		names(l2),
		FUN = function(z) {
			x  = l2[[z]]
			r  = list()
			nl = names(x)[!names(x) %in% oblCols]
			if (length(nl) > 1) {
				cmbn = combn(nl, length(nl) - 1)
				for (i in 1:ncol(cmbn)) {
					r[[i]] = lapply(i, function(y) {
						x[, -which(names(x) %in% cmbn[, y]), drop = FALSE]
					})
					names(r[[i]]) = paste(z, paste0(cmbn[, i], collapse = '....'), sep = '....')
				}
				return(unlist(r, recursive = FALSE))
			} else{
				return(x)
			}
		}
	)
	names(l3)     = names(l2)
	# Which sublists are sublevels?
	message0('\tFinalising the tables ....')
	k = unlist(lapply(l3, function(x) {
		all(sapply(x, is.list))
	}), recursive = FALSE)
	f = function(l) {
		names(l) = NULL
		r = unlist(l, recursive = FALSE, use.names = TRUE)
		return(r)
	}
	if (any(k)) {
		l3           = c(l3[!k], f(l3[k]))
	}
	
	return(l3)
}


FormulaHasIntercept = function(formula) {
	if (is.null(formula))
		return(NULL)
	attr(terms(as.formula(formula))    , which = 'intercept') > 0
}

FormulaHasResponse = function(formula) {
	if (is.null(formula))
		return(NULL)
	attr(terms(as.formula(formula))    , which = 'response') > 0
}
getResponseFromFormula = function(formula) {
	if (is.null(formula))
		return(NULL)
	if (attr(terms(as.formula(formula))    , which = 'response'))
		all_vars0(formula)[1]
	else
		NULL
}


formulaTerms = function(formula,
												response  = FALSE,
												intercept = FALSE) {
	if (!is.null(formula)) {
		r = c(
			if (response && FormulaHasResponse(formula))
				getResponseFromFormula(formula)
			else
				NULL,
			attr(terms(as.formula(formula))        , which = 'term.labels'),
			if (intercept &&
					FormulaHasIntercept(formula))
				1
			else
				NULL
		)
	}
	else{
		r = NULL
	}
	return(r)
}

expand.formula = function(formula) {
	reformulate(termlabels = labels(terms(formula)),
							response =
								if (attr(terms(formula), "response") > 0)
									formula[[2]]
							else
								NULL)
}

UnlistCall = function(x) {
	as(unlist(x), 'character')
}

ListOperation <- function(x, FUN = NULL)
{
	cnames <- names(x)
	if (is.null(cnames))
		return (x)
	x1 <- lapply(cnames, function(y)
		ListOperation(x[[y]]))
	x1 <- FUN(x1)
	
	return(x1)
}

# You may want to use list.clean from rlist package
prunelist = function(x) {
	message0('Pruning list in progress ...')
	r = lapply(x, function(y) {
		if (is.null(y) || length(y) < 1) {
			NULL
		} else{
			y
		}
	})
	return(r)
}


ModelInReference = function(model,
														reference,
														responseIncluded = FALSE,
														veryLower = ~ Genotype + 1) {
	mo  = formulaTerms(formula = model, intercept = TRUE)
	re  = formulaTerms(formula = reference, intercept = TRUE)
	r   = re[re %in% mo]
	if (length(r) > 0) {
		out = reformulate(
			termlabels = r,
			response   = if (responseIncluded && FormulaHasResponse(model))
				all_vars0(model)[1]
			else
				NULL	,
			intercept  = TRUE
		)
		if (length(mo[!(mo %in% re)]) > 0) {
			message0('Some terms in the "lower" model are ignored. See:\n\t',
							 pasteComma(mo[!(mo %in% re)]))
			message0('The polished "lower": ', printformula(out))
		}
	} else{
		message0('An invalid "lower". It is set to: ', printformula(veryLower))
		out = veryLower
	}
	return(out)
}

TermInFormulaReturn = function(formula,
															 term,
															 return,
															 active,
															 not = NA,
															 debug = TRUE) {
	terms  = formulaTerms(formula = formula)
	if (debug)
		message0(printformula(terms))
	if (active && (term %in% terms)) {
		return(return)
	} else{
		return(not)
	}
}

paste01 = function(...) {
	r = paste0(...)
	while (any(grepl(pattern = '  ', x = r))) {
		r = gsub(
			pattern = '  ',
			x = r,
			replacement =  ' ',
			fixed = TRUE
		)
	}
	return(r)
}

greplM = function(x = NULL, pattern = NULL, ...) {
	if (is.null(x))
		return(NULL)
	if (is.null(pattern))
		return(x)
	
	r = rep(TRUE, length(x))
	for (p in pattern)
		r = r & grepl(pattern = p, x = x, ...)
	return(r)
}

modelContrasts = function(formula, data, ...) {
	if (is.null(formula) || is.null(data))
		return(NULL)
	r = colnames(model.matrix(formula , data, ...))
	return(r)
}

multiBatch = function(data) {
	if ('Batch' %in% colnames(data)) {
		batchColumn <- na.omit(data[, "Batch"])
		if (length(levels(batchColumn)) > 1)
			TRUE
		else
			FALSE
	}
	else
		FALSE
}

unmatrix0 = function (x, byrow = FALSE, sep = ' ')
{
	if (is.null(x))
		return(NULL)
	
	rnames <- rownames(x)
	cnames <- colnames(x)
	if (is.null(rnames))
		rnames <- paste("r", 1:nrow(x), sep = "")
	if (is.null(cnames))
		cnames <- paste("c", 1:ncol(x), sep = "")
	nmat <- outer(rnames, cnames, paste, sep = sep)
	if (byrow) {
		vlist <- c(t(x))
		names(vlist) <- c(t(nmat))
	}
	else {
		vlist <- c(x)
		names(vlist) <- c(nmat)
	}
	return(vlist)
}

MoveResponseToRightOfTheFormula = function(formula) {
	if (is.null(formula)) {
		message0('Ops! the formula is blank')
		return(NULL)
	}
	newFormula = update.formula(old = formula,
															new = reformulate(
																response   = NULL,
																termlabels = c(all_vars0(formula)[1], '.')
															))
	out = formula(delete.response(terms(newFormula)))
	message0('The input formula: ', printformula(formula))
	message0('The reformatted formula for the algorithm: ',
					 printformula(out))
	return(out)
}

renameVariableNameInList = function(list,
																		name,
																		replace = NULL,
																		prefix,
																		not = FALSE) {
	if (is.null(list))
		return(NULL)
	
	if (is.null(replace))
		if (!not) {
			names(list)[names(list) %in% name] = paste(prefix, names(list)[names(list) %in% name], sep = '_')
		} else{
			names(list)[!names(list) %in% name] = paste(prefix, names(list)[!names(list) %in% name], sep = '_')
		}
	else
		if (!not) {
			names(list)[names(list) %in% name] = replace
		} else{
			names(list)[!names(list) %in% name] = replace
		}
	
	return(list)
}

decimalplaces <- function(x) {
	if ((x %% 1) != 0) {
		nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed = TRUE)[[1]][[2]])
	} else {
		return(0)
	}
}

dataSignature = function(formula, data, digits = 10) {
	a.vars = all_vars0(formula)
	if (!is.null(formula) &&
			!is.null(data)    &&
			!is.null(a.vars)  &&
			nrow(data) > 1    &&
			sum(a.vars %in% names(data))) {
		v.vars = a.vars[a.vars %in% names(data)]
		res0 = lapply(v.vars, function(name) {
			d = data[, name]
			if (is.numeric(d))
				paste0(
					name,
					':[',
					pasteComma(
						paste0('n='   , length(d))                            ,
						paste0('mean=', round(mean(d, na.rm = TRUE), digits)) ,
						paste0('sd='  , round(sd0 (d, na.rm = TRUE), digits)) ,
						truncate      = FALSE                                 ,
						trailingSpace = FALSE
					),
					']'
				)
			else
				paste0(name,
							 ':[',
							 pasteComma(
							 	paste0(levels(d), '=', table(d)),
							 	truncate           = FALSE      ,
							 	trailingSpace      = FALSE
							 ),
							 ']')
		})
		
		overal  = mean(rowMeans(apply(data[, v.vars, drop = FALSE], 2, function(x) {
			r = if (!is.numeric(x)) {
				as.integer(as.factor(x))
			} else{
				x
			}
			return(r)
		}), na.rm = TRUE))
		
		res0$overal    = paste0('Overall:['  , overal , ']')
		res0$precision = paste0('Precision:[', digits , ']')
		res = pasteComma(
			sort(unlist(res0), decreasing = FALSE) ,
			truncate      = FALSE                  ,
			trailingSpace = FALSE
		)
	} else{
		res = paste0(
			'Something is wrong with the data/model. Make sure that the data is not null and the formula matches the data. [This is a random number ',
			randRegSeed(
				n       = 1   ,
				decimal = TRUE,
				round   = 10
			),
			']'
		)
	}
	return(res)
}

randRegSeed = function(n = 1,
											 max = .Machine$double.xmax,
											 decimal = TRUE,
											 round = 10) {
	r = runif(n, 1, as.numeric(Sys.time()) * 10000) %% max
	if (decimal)
		r = r - (r %% 1)
	r = round(r, digits = round)
	return(r)
}

as.numeric01 = function(x) {
	if (!is.null(x))
		return(suppressWarnings(as.numeric(x)))
	else
		return(NULL)
}

PhenListAgeingLevels = function(object, sep = '') {
	l = NULL
	SexLab    = 'Sex'
	FemaleLab = ifelse(
		is.na(object$input$PhenListAgeing@dataset.values.female),
		'Female',
		object$input$PhenListAgeing@dataset.values.female
	)
	MaleLab = ifelse(
		is.na(object$input$PhenListAgeing@dataset.values.male),
		'Male',
		object$input$PhenListAgeing@dataset.values.male
	)
	SexLevels = as.list(paste(SexLab, sort(c(
		FemaleLab, MaleLab
	), decreasing = FALSE),
	sep =	sep))
	names(SexLevels) = unlist(SexLevels)
	##########
	GenotypeLab = 'Genotype'
	ControlLab = ifelse(
		is.na(object$input$PhenListAgeing@refGenotype),
		'control',
		object$input$PhenListAgeing@refGenotype
	)
	MutantLab = ifelse(
		is.na(object$input$PhenListAgeing@testGenotype),
		'experimental',
		object$input$PhenListAgeing@testGenotype
	)
	GenotypeLevels = as.list(paste(GenotypeLab, sort(c(
		ControlLab, MutantLab
	), decreasing = TRUE),
	sep =	sep))
	names(GenotypeLevels) = unlist(GenotypeLevels)
	
	##########
	BatchLab = 'Batch'
	##########
	WeightLab = 'Weight'
	##########
	l = list(
		response  = object$input$depVariable	,
		LifeStage = list(
			LifeStage    = 'LifeStage'          ,
			Early        = 'Early'              ,
			Late         = 'Late'               ,
			Levels = list(LifeStageEarly = 'LifeStageEarly', LifeStageLate = 'LifeStageLate')
		),
		Sex = list(
			Sex    = SexLab   ,
			Female = FemaleLab,
			Male   = MaleLab  ,
			Levels = as.list(SexLevels)
		),
		Genotype   = list(
			Genotype = GenotypeLab,
			Control  = ControlLab ,
			Mutant   = MutantLab  ,
			Levels   = as.list(GenotypeLevels)
		),
		Batch  = BatchLab       ,
		Weight = WeightLab
	)
	return(l)
}

CombineLevels = function(..., debug = TRUE) {
	r = c(paste(..., sep = ':'))
	if (debug)
		print(r)
	return(r)
}

RRNewObjectAndFormula = function(object              ,
																 RRprop              ,
																 formula             ,
																 labels       = NULL ,
																 depVarPrefix = NULL ,
																 refLevel     = NULL ,
																 ### see right in the cut function
																 right        = TRUE) {
	allTerms         = all_vars0(formula)
	newobject        = object
	RRcutObject  = RRCut(
		object         = object        ,
		prob           = RRprop        ,
		depVariable    = allTerms[1]   ,
		labels         = labels        ,
		depVarPrefix   = depVarPrefix  ,
		right          = right         ,
		refLevel       = refLevel      ,
		lower          = allTerms[2]
	)
	if (is.null(RRcutObject))
		return(NULL)
	
	newobject        = RRcutObject$discObject
	newFormula       = replaceElementInFormula(
		formula = formula,
		pattern = allTerms[1] ,
		replace = RRcutObject$newdepVariable
	)
	return(
		list(
			newobject      = newobject                  ,
			newFormula     = newFormula                 ,
			newDepVariable = RRcutObject$newdepVariable ,
			object         = object                     ,
			formula        = formula                    ,
			depVariable    = allTerms[1]                ,
			###
			RRprop         = RRprop                     ,
			labels         = labels                     ,
			depVarPrefix   = depVarPrefix               ,
			refLevel       = RRcutObject$refLevel       ,
			empiricalQuantiles = RRcutObject$percentages
		)
	)
}

jitter0 = function(x,
									 factor = 1,
									 amount = NULL,
									 upper = 1,
									 lower = .5,
									 maxtry = 1500) {
	for (i in 1:maxtry) {
		if (i >= maxtry)
			message0('\tNo solusion found for the specified RR_prop.')
		xx = jitter(x = x,
								amount = amount,
								factor = factor)
		if (min(xx, na.rm = TRUE) > lower && max(xx, na.rm = TRUE) < upper)
			break
	}
	return(xx)
}

addJitterToTheEntireData = function(x, min = -1, max = 0) {
	if (is.null(x) || length(na.omit(x)) < 1)
		return(x)
	lx     = length(x)
	jitter = sort(runif(n = lx, min = min, max = max), decreasing = FALSE)
	message0(
		'A small jitter (max value = ',
		min(jitter, na.rm = TRUE) -
			max(jitter, na.rm = TRUE), 
		') is added to the data')
	o      = order(x)
	x      = x + jitter[o]
	return(x)
}

ExpandTails = function(x,
											 amount,
											 upper = TRUE,
											 lower = TRUE) {
	xx = na.omit(x)
	if (is.null (x)   ||
			is.null (xx)  ||
			length  (xx) < 1)
		return(x)
	if (upper)
		x[x == max(xx, na.rm = TRUE)] = max(xx, na.rm = TRUE) + amount
	if (lower)
		x[x == min(xx, na.rm = TRUE)] = min(xx, na.rm = TRUE) - amount
	return(x)
}

ReplaceTails = function(x         ,
												lower = 0 ,
												upper = 0 ,
												add   = FALSE) {
	xx = na.omit(x)
	if (is.null (x)   ||
			is.null (xx)  ||
			length  (xx) < 1)
		return(x)
	if (!add) {
		if (upper)
			x[which.max(x)] = upper[1]
		if (lower)
			x[which.min(x)] = lower[1]
	} else{
		x = c(lower, x, upper)
	}
	return(sort(unname(x)))
}

CheckTheValidityOfTheRRLower = function(lower, data, depvar, minLevels = 2) {
	if (is.null(data) || is.null(lower)) {
		stop('~> Null data or the `Reference variable`. Please check the input data or the `Reference variable`')
	}
	if (is.null(depvar) || !depvar %in% names(data)) {
		stop('~> dependent variable does not exist in the data')
	}
	if (!is.numeric(data[, depvar])) {
		stop('~> dependent variable must be numeric')
	}
	if (is.null(lower) || !lower %in% names(data)) {
		stop('~> `Reference variable` does not exist in the data')
	}
	if (is.numeric(data[, lower])) {
		stop('~> `Reference variable` must be a factor')
	}
	if (nlevels(as.factor(data[, lower])) < minLevels) {
		stop('~> `Reference variable` must have at least two levels')
	}
	return(TRUE)
}

RRGetTheLeveLWithMaxFrequency = function(data, lower) {
	tbl = table(data[, lower])
	maxTbl = tbl %in% max(tbl, na.rm = TRUE)[1]
	r      = names(tbl[maxTbl])
	if (length(r) > 1) {
		message0(
			'\tMore than one variable with the highest frequency detected. See:\n\tLevel(frequency): ',
			pasteComma(paste0(r, '(', tbl[maxTbl], ')'), truncate = FALSE),
			'\n\t\tThe first one (`',
			r[1],
			'`) would be used.'
		)
	} else{
		message0('\t\tNominated level(frequency) = ', 	pasteComma(paste0(r, '(', tbl[maxTbl], ')')))
	}
	return(r[1])
}

ExtractDatasetPLIfPossible = function(x) {
	if (class(x) %in% c('PhenList', 'PhenListAgeing')) {
		x     = x@datasetPL
	}
	return(x)
}

RRCut = function(object                     ,
								 prob         = .95         ,
								 depVariable  = 'data_point',
								 labels       = NULL,
								 depVarPrefix = NULL,
								 right        = TRUE,
								 refLevel     = NULL,
								 lower        = 'Genotype'
) {
	data     = object
	if (class(data) %in% c('PhenList', 'PhenListAgeing')) {
		refLevel = data@refGenotype
		data     = data@datasetPL
	}
	if (!CheckTheValidityOfTheRRLower(lower = lower, data = data,depvar = depVariable)) {
		return(NULL)
	}
	if (prob == .5) {
		message0('\t`prop` must be different from 0.5')
		return(NULL)
	}
	if(is.null(refLevel)){
		message0('\tReference level left blank, then the dominate level will be set as the reference level')
		refLevel = RRGetTheLeveLWithMaxFrequency(data = data, lower = lower)
	}else{
		message0('Reference level is set to `', refLevel,'`')
	}
	# Preparation ...
	JitterPrecision = 4 + 1 * decimalplaces(min(data[, depVariable], na.rm = TRUE))
	data$data_point_discretised  = NA
	controls = subset(data,
										data[,lower] %in% refLevel)
	mutants  = subset(data,
										!(data[,lower] %in% refLevel))
	prb      = unique(c(0,  prob, 1))
	qntl     = ReplaceTails(
		x      = quantile(x = controls[, depVariable],
											probs = prb                ,
											na.rm = TRUE               ),
		upper = max(data[, depVariable], na.rm = TRUE)[1] + 1,
		lower = min(data[, depVariable], na.rm = TRUE)[1] - 1,
		add   = FALSE
	)
	if (sum(duplicated(qntl))) {
		message0(
			'\t* duplicates in quantiles detected, then small (dependes on the data precision) jitters will be added to quantiles.\n\t\tQauntiles: ',
			pasteComma(round(qntl,5), replaceNull = FALSE)
		)
		message0('\tJitter precision (decimal) = ', JitterPrecision)
		qntl[duplicated(qntl)] = jitter0(
			x = qntl[duplicated(qntl)],
			amount = 10 ^	-JitterPrecision,
			upper  = 1,
			lower  = 0.5
		)
		qntl = sort(qntl)
	}
	message0(
		'\tInitial quantiles for cutting the data ',
		'[probs = ('                               ,
		pasteComma(round(prb, 3))                  ,
		'), n.reference = '                        ,
		length(controls[, depVariable]),
		']: '                                    ,
		pasteComma(round(qntl, 3), replaceNull = FALSE)
	)
	if (length(unique(qntl)) < 2) {
		message0('\tThe algorithm cannot specify non-unique quantiles.')
		return(NULL)
	}
	controls$data_point_discretised  = cut(
		x = controls[, depVariable],
		breaks = qntl,
		labels = if (is.null(labels))
			FALSE
		else
			labels,
		include.lowest = TRUE ,
		right          = right 
	)
	###
	tbc  = prop.table(table(controls$data_point_discretised))
	tbpc = setNames(as.list(tbc),names(tbc))
	message0('\tDetected percentile in the data: ',pasteComma(paste(names(tbc),'=',round(tbc,6))))
	###
	mutants$data_point_discretised   = cut(
		x = mutants [, depVariable],
		breaks = qntl,
		labels = if (is.null(labels))
			FALSE
		else
			labels,
		include.lowest = TRUE,
		right          = right
	)
	newObj                          = rbind(mutants, controls)
	newdepVariable                  = pasteUnderscore(c(depVarPrefix, depVariable, 'discretised'))
	names(newObj)[names(newObj) %in% 'data_point_discretised'] = newdepVariable
	newObj[, newdepVariable]        = as.factor(newObj[, newdepVariable])
	#message0('\tA new column is added to the data object: ', newdepVariable)
	return(
		list(
			object      = object             ,
			discObject  = newObj             ,
			newdepVariable = newdepVariable  ,
			depVariable = depVariable        ,
			percentages = tbpc               ,
			refLevel    = refLevel           ,
			lower       = lower    
		)
	)
}

RRDiscretizedEngine = function(data,
															 formula      = data_point ~ Genotype + Sex + zygosity ,
															 depVar       = 'data_point'                           ,
															 lower        = 'Genotype'                             ,
															 refLevel     = NULL                                   ,
															 labels       = c('Low', 'NormalHigh')                 ,
															 depVarPrefix = 'Low'                                  ,
															 right        = TRUE                                   ,
															 prob         = .95) {
	requireNamespace("rlist")
	l1 = l2 = l3 = l4 = NULL
	if (!CheckTheValidityOfTheRRLower(lower  = lower,
																		depvar = depVar,
																		data   = data)) {
		return(NULL)
	}
	vars  =  names(data) %in% all.vars(formula)
	df    =        data[, vars]
	cat   = !sapply(df[, all.vars(formula)], is.numeric)
	extra = all.vars(formula)[cat &
															!all.vars(formula) %in% c(depVar, lower)]
	lextra = length(extra)
	message0('Preparing the reference ranges ...')
	message0('Preparing the data for the variable: ', lower)
	
	l1 = RRNewObjectAndFormula(
		object  = data                ,
		RRprop  = prob                ,
		formula = formula             ,
		labels  = labels              ,
		depVarPrefix =  depVarPrefix  ,
		refLevel     = refLevel       ,
		right        = right
	)
	
	if (lextra > 0) {
		for (i in 1:lextra) {
			cbn = combn(extra, i)
			for (j in 1:ncol(cbn)) {
				message0('\tspliting on ', pasteComma(cbn[, j], replaceNull = FALSE))
				out = split(df,
										interaction(df[cbn[, j]]), drop = TRUE)
				l2   = c(l2, out)
			}
		}
	}
	if (!is.null(l2)) {
		l3 = lapply(names(l2), function(name) {
			message0('Preparing the data for the combined effect: ', name)
			x = l2[[name]]
			out = RRNewObjectAndFormula(
				object       = x            ,
				RRprop       = prob         ,
				formula      = formula      ,
				labels       = labels       ,
				depVarPrefix = depVarPrefix ,
				refLevel     = refLevel     ,
				right        = right
			)
		})
		names(l3) = names(l2)
	}
	l4        = c(list(l1), l3)
	#### must improve in future
	if (!is.null(l4)) {
		names(l4) =	ifelse(
			nzchar(names(l4))                         ,
			paste(depVar, lower, names(l4), sep = '.'),
			paste(depVar, lower,            sep = '.')
		)
	}
	l4        = list.clean(l4)
	return(l4)
}

RRextra = function(object,
									 prob = .95,
									 depVariable = 'data_point') {
	if (prob <= .5) {
		message0('"prob" must be greater than 0.5')
		return(NULL)
	}
	# Preparation ...
	controls = subset(object@datasetPL,
										object@datasetPL$Genotype %in% object@refGenotype)
	mutants  = subset(object@datasetPL,
										object@datasetPL$Genotype %in% object@testGenotype)
	prb      = unique(c(0, 1 - prob, prob, 1))
	qntl     = quantile(x = controls[, depVariable], probs = prb)
	cutsC    = cut(x = controls[, depVariable], breaks = qntl)
	
	message0('Creating control cuts ...')
	XclassValue = tapply(controls[, depVariable], cutsC , function(x) {
		(x)
	})
	message0('Creating mutant cuts ...')
	MclassValue = lapply(
		XclassValue,
		FUN = function(xx) {
			sum(mutants$data_point %in% xx)
		}
	)
	CclassValue = lapply(XclassValue, length)
	message0('Creating output tables ...')
	# Overal Table
	tbl = rbind(unlist(CclassValue), unlist(MclassValue))
	dimnames(tbl) = list(c('Control', 'Mutant'), c('Low', 'Normal', 'High'))
	# Table Low
	tblLow = cbind(tbl[, 1], rowSums(tbl[, 2:3]))
	dimnames(tblLow) = list(c('Control', 'Mutant'), c('Low', 'Normal/High'))
	
	tblHigh = cbind(tbl[, 1], rowSums(tbl[, 1:2]))
	dimnames(tblHigh) = list(c('Control', 'Mutant'), c('Low/Normal', 'High'))
	
	return(list(
		overall = tbl,
		tblLow = tblLow,
		tblHigh = tblHigh
	))
}

TermInModelAndnLevels = function(model,
																 term = 'LifeStage',
																 data,
																 threshold = 1) {
	if (is.null(model) ||
			is.null(data) || is.null(term) || !(term %in% colnames(data)))
		return(FALSE)
	# if (!is.null(model$correctd))
	# 	model      = model$correctd
	r = termInTheModel(model = model ,
										 term = term,
										 message = FALSE) &&
		colLevelsSimple(data, all_vars0(model)[1]) >		threshold
	return(r)
}
#######################
# From PhenStat
columnChecks0 = function(dataset,
												 columnName,
												 dataPointsThreshold = 4) {
	presence <- TRUE
	numeric <- FALSE
	levelsCheck <- 0
	variabilityThreshold <- 10
	# Test: dependent variable presence
	if (!(columnName %in% colnames(dataset))) {
		presence <- FALSE
	}
	else {
		columnOfInterest <- na.omit(dataset[, c(columnName)])
		
		if (all(sapply(columnOfInterest, is.numeric))) {
			numeric <- TRUE
		}
		
		dataPointsSummary <- columnLevels(dataset, columnName)
		
		NoCombinations <- dataPointsSummary[3]
		variabilityThreshold <- NoCombinations
		for (i in 1:NoCombinations) {
			if (dataPointsSummary[3 + i] >= dataPointsThreshold)
				levelsCheck <- levelsCheck + 1
		}
	}
	
	values <-
		c(presence, numeric, (levelsCheck >= variabilityThreshold))
	
	return (values)
}

colLevelsSimple = function(dataset, colName) {
	return(length(unique(dataset[, colName])))
}

columnLevels = function(dataset, columnName) {
	columnOfInterest <- na.omit(dataset[, c(columnName)])
	
	
	values <- c(length(columnOfInterest))
	
	#Test for the data points quantity for Genotype/sex combinations
	Genotype_levels <- levels(factor(dataset$Genotype))
	Sex_levels <- levels(factor(dataset$Sex))
	values <- append(values, length(levels(factor(columnOfInterest))))
	
	values <-
		append(values, length(Genotype_levels) * length(Sex_levels))
	
	for (i in 1:length(Genotype_levels)) {
		GenotypeSubset <-
			subset(dataset, dataset$Genotype == Genotype_levels[i])
		for (j in 1:length(Sex_levels)) {
			GenotypeSexSubset <- subset(GenotypeSubset,
																	GenotypeSubset$Sex == Sex_levels[j])
			
			columnOfInterestSubset <-
				na.omit(GenotypeSexSubset[, c(columnName)])
			
			values <- append(values, length(columnOfInterestSubset))
			
		}
	}
	return (values)
}

# From  capitalize {Hmisc}
capitalise = function (string)
{
	capped <- grep("^[A-Z]", string, invert = TRUE)
	substr(string[capped], 1, 1) <- toupper(substr(string[capped],
																								 1, 1))
	return(string)
}

sort0 = function(x, ...) {
	if (length(x) > 0)
		x = sort(x = x, ...)
	return(x)
}

message0 = function(...,
										breakLine  = TRUE,
										capitalise = TRUE,
										appendLF   = TRUE,
										active     = TRUE) {
	if (active) {
		x = paste0(..., collapse = '')
		if (breakLine)
			nmessage = unlist(strsplit(x = x, split = '\n'))
		else
			nmessage = x
		if (capitalise)
			nmessage = capitalise(nmessage)
		message(paste(Sys.time(), nmessage, sep = '. ', collapse = '\n'),
						appendLF = appendLF)
	}
}

warning0 = function(...,
										breakLine = TRUE,
										capitalise = TRUE) {
	x = paste0(..., collapse = '')
	if (breakLine)
		nmessage = unlist(strsplit(x = x, split = '\n'))
	else
		nmessage = x
	if (capitalise)
		nmessage = capitalise(nmessage)
	warning(paste(Sys.time(), nmessage, sep = '. ', collapse = '\n'))
}


extractFisherSubTableResults = function(x, what = 'p.value') {
	r = lapply0(x, function(y)
		y[what])
	return(r)
}

lapply0 = function(X, FUN, ...) {
	if (is.null(X))
		return(NULL)
	r = lapply(X = X, FUN = FUN, ...)
	return(r)
}

all_vars0 = function(x, ...) {
	if (is.null(x))
		return(NULL)
	fif = all.vars(formula(x), ...)
	if (length(fif) > 0) {
		return(fif)
	} else{
		return(NULL)
	}
}

intervalsCon = function(object, lvls, ...) {
	if (is.null(object))
		return(object)
	
	message0('\tComputing the confidence intervals at the level of ',
					 pasteComma(lvls),
					 ' ...')
	ci = lapply(lvls, function(x) {
		citerms = if (is(object, 'lme')) {
			c('all', 'fixed', 'var-cov')
		} else if (is(object, 'gls')) {
			c('all', 'coef', 'var-cov')
		} else{
			c('all')
		}
		for (citerm in citerms) {
			intv = tryCatch(
				expr = intervals(
					object = object,
					level  = x,
					which  = citerm,
					...
				),
				error = function(e) {
					message0('\t\t ~> Error in estimating the confidence intervals for `',
									 citerm,
									 '` term(s)')
					#message0('\t\t ~> ', e, breakLine = FALSE)
					return(NULL)
				} ,
				warning = function(w) {
					message0('\t\t ~> Error in estimating the confidence intervals for `',
									 citerm,
									 '` term(s)')
					#message0('\t\t ~> ', w, breakLine = FALSE)
					return(NULL)
				}
			)
			if (!is.null(intv)) {
				message0('\t CI for `', citerm, '` term(s) successfully estimated')
				break
			} else if (!citerm %in% tail(citerms, 1)) {
				message0('\tAdjustment applied. Retrying ....')
			} else{
				message0('CI estimation failed.')
			}
		}
		if (is.null(intv)) {
			return(NULL)
		} else{
			return(list(intervals = intv, level = x))
		}
	})
	if (!is.null(ci))
		names(ci) = paste('CI_', lvls, sep = '')
	object$intervals = ci
	return(object)
}



intervalsCat = function(object, lvls = .95, ...) {
	if (is.null(object))
		return(object)
	
	# message0('\t\tReformatting the confidence interval in the level of ',
	# 				 pasteComma(lvls),
	# 				 ' ...')
	c0 = object$conf.int
	c2 = NULL
	if (!is.null(c0)) {
		ci = list(
			intervals = list(
				lower = c0[1]                     ,
				est.  = as.vector(object$estimate),
				upper = c0[2]
			)                                   ,
			level   = attr(c0, "conf.level")
		)
	} else{
		ci = NULL
	}
	c2$intervals = ci
	if (!is.null(ci))
		names(c2)  = paste('CI_', lvls, sep = '')
	object$interval = c2
	return(object)
}

expandDottedFormula = function(formula, data) {
	if (is.null(formula) || is.null(data)) {
		message0('Null formula or data')
		return(formula)
	}
	
	newformula  =  formula(terms.formula(x = formula, data = data))
	message0('Extended formula (if (*.:) characters included):\n\t ',
					 printformula(newformula))
	return(newformula)
}

removeSingleLevelFactors = function(formula, data) {
	cat    = all_vars0(formula)[!sapply(data[, all_vars0(formula)], is.numeric)]
	if (length(cat)) {
		FactsThatMustBeRemoved = cat[lapply(data[, cat, drop = FALSE], function(x) {
			length(unique(na.omit(x)))
		}) <= 1]
		if (length(FactsThatMustBeRemoved)) {
			message0(
				'The following terms from the model are removed because they only contain one level: ',
				pasteComma(
					FactsThatMustBeRemoved,
					replaceNull = FALSE,
					truncate = FALSE
				)
			)
			formula = update.formula(formula,
															 reformulate0(
															 	termlabels = c('.', FactsThatMustBeRemoved),
															 	response   = NULL,
															 	intercept  = TRUE,
															 	sep        = '-'
															 ))
		}
	}
	return(formula)
}


modelSummaryPvalueExtract = function(x,
																		 variable   = 'Genotype'                          ,
																		 anova      = TRUE                                ,
																		 what       = c('Pr(>|z|)', 'Pr(>Chi)', 'p-value'),
																		 debug      = TRUE,
																		 ci_display = FALSE) {
	if (is.null(x))
		return(NULL)
	if (anova) {
		if (any(class(x) %in% 'glm')) {
			mOrg = anova(x, test = 'LRT')
		} else{
			mOrg = anova(x)
		}
	} else{
		mOrg  = coef(summary(x))
	}
	mOrg  = as.data.frame(mOrg)
	if (debug)
		print(mOrg)
	mSum  =  mOrg[, colnames(mOrg) %in% what, drop = FALSE]
	if (!is.null(variable) && any(variable %in% rownames(mOrg))) {
		mSumFiltered  = mSum[rownames(mOrg) %in% variable, , drop = FALSE]
	} else{
		mSumFiltered = NULL # NA?
	}
	if (ci_display && !is.null(mSumFiltered)) {
		fixedEffInters = CheckWhetherNamesExistInListExtractCI(
			x        = x$intervals[[1]]$intervals,
			lnames   = c('coef', 'fixed')        ,
			variable = variable
		)
		return(list(
			'Value'      = as.vector(unlist(mSumFiltered))     ,
			'Confidence' =  Matrix2List(x = fixedEffInters)    ,
			'Level'      = ifelse(is.null(fixedEffInters), NULL, x$intervals[[1]]$level)
		))
	} else{
		return(as.vector(unlist(mSumFiltered)))
	}
}

CheckWhetherNamesExistInListExtractCI = function(x, lnames, variable, minusCol = 2) {
	if (!is.null(x) && any(names(x) %in% lnames)) {
		Toplst = x[names(x) %in% lnames][[1]]
		Toplst = Toplst[rownames(Toplst) %in% variable , -minusCol,
										drop = FALSE]
		return(Toplst)
	} else{
		return(NULL)
	}
}

CatEstimateAndCI = function(object) {
	if (is.null(object))
		return(NULL)
	v = object$interval[[1]]
	if (!is.null(v)) {
		return(list(
			'Value' = v$intervals$est.    ,
			'Confidence' = list(
				'Lower' = v$intervals$lower ,
				'Upper' = v$intervals$upper
			)                             ,
			'Level'   = v$level
		))
	} else{
		return('Only available for 2x2 tables')
	}
}

NormaliseDataFrame = function(data           ,
															colnames = NULL) {
	message0('Normalising the data.frame in progress ...')
	
	if (!is.null(colnames)) {
		colnames = colnames[colnames %in% names(data)]
		
		if (length(colnames) < 1) {
			message0('No variable name found in the data. Please check "colnames" ...')
			return(data)
		}
	} else{
		message0(
			'No variable selected for normalisation. All numerical variables will be normalised.'
		)
		colnames = names(data)
	}
	######
	data  [, colnames] = as.data.frame(lapply(
		data[, colnames, drop = FALSE],
		FUN = function(x) {
			if (is.numeric(x) && length(unique(x)) > 1) {
				sdx = sd0(x, na.rm = TRUE)
				r   =    (x - mean(x, na.rm = TRUE)) / ifelse(!is.na(sdx) &&
																												sdx > 0, sdx, 1)
			} else{
				r = x
			}
			return(r)
		}
	))
	return(data)
}

getlmeObjectConfigs = function(obj) {
	l    = list(
		'Formula' = printformula(if (class(obj) %in% 'lme') {
			obj$call$fixed
		} else{
			obj$call$model
		}, message = FALSE),
		'Random effect'         = printformula(obj$call$random, message = FALSE)
	)
	return(l)
}

extractLmeTerms = function(object) {
	if (is.null(object) || !is.null(object$messages))
		return(NULL)
	r = list(
		'Initial model' = getlmeObjectConfigs(object$output$Initial.Model),
		'Final model'   = getlmeObjectConfigs(object$output$Final.Model)
	)
	return(r)
}

extractFERRTerms = function(object) {
	if (is.null(object) || !is.null(object$messages))
		return(NULL)
	lnitial = final = list(
		formula               = printformula(object$input$formula, message = FALSE),
		random_effect         = NULL
	)
	final$formula = printformula(RightFormula2LeftFormula(object$extra$Cleanedformula),
															 message = FALSE)
	out = list('Initial model' = lnitial,
						 'Final model'   = final)
	
	if (!is.null(object$input$RRprop))
		out$'RR quantile' = unlist(object$input$RRprop)
	
	return(out)
}

RightFormula2LeftFormula = function(formula, removeResIfExists = FALSE) {
	r = if (!is.null(formula) &&
					length(all_vars0(formula)) > 0) {
		if (removeResIfExists && FormulaHasResponse(formula)) {
			formula = update(formula, NULL ~ .)
		}
		
		update(
			as.formula(formula),
			reformulate0(
				termlabels = c('.', all_vars0(formula)[1]),
				response = all_vars0(formula)[1],
				sep = '-'
			)
		)
	} else{
		NULL
	}
	return(r)
}

sd0 = function(x, ...) {
	if (!is.numeric(x))
		return(NA)
	r = if (length(na.omit(x)) > 1)
		sd(x, ...)
	else
		0
	return(r)
}

SummaryStats = function(x,
												formula                              ,
												#label = 'raw_data_summary_statistics',
												lower = FALSE                        ,
												drop = TRUE                          ,
												sep = '_'                            ,
												removeSpecialChars = FALSE           ,
												replace = '_') {
	r       = NULL
	formula = checkModelTermsInData(formula = formula,
																	data = x,
																	responseIsTheFirst = TRUE)
	depVar  = all_vars0(formula)[1]
	if (is.null(depVar)) {
		message0('Null response! check the formula and the data')
		return(NULL)
	}
	# do not move me
	if (any(dim(x) == 0))
		return('empty dataset')
	
	cat    = all_vars0(formula)[!sapply(x[, all_vars0(formula)], is.numeric)]
	if (length(cat) > 0)
		lvls   = interaction(x[, cat], sep = sep, drop = drop)
	else
		lvls = rep(depVar, nrow(x))
	
	isNumeric  = is.numeric(x[, depVar])
	summaryT   = as.list(tapply(x[, depVar], INDEX = lvls, function(xx) {
		if (isNumeric) {
			c  = ifelse(length(na.omit(xx)) > 0, length(na.omit(xx)), 0)
			m  = ifelse(length(na.omit(xx)) > 0, mean(xx, na.rm = TRUE), NA)
			sd = ifelse(length(na.omit(xx)) > 0, sd0 (xx, na.rm = TRUE), NA)
			r = list(
				'Count' = c                       ,
				'Mean' = m                        ,
				'SD' = sd                         ,
				'Normality test' = head(shapiro.test0(xx), 3)
			)
		} else{
			c = ifelse(length(na.omit(xx)) > 0, length(na.omit(xx)), 0)
			r = list('Count' = c,
							 'Mean' = NULL,
							 'SD' = NULL)
		}
		return(r)
	}, default = -999.991233210123))
	##
	fTmp = function(isNum) {
		if (isNum) {
			r = list('Count' = 0,
							 'Mean' = NA,
							 'SD' = NA)
		} else{
			r = list('Count' = 0)
		}
		return(r)
	}
	summaryT[summaryT %in% c(-999.991233210123)] = fTmp(isNum = isNumeric)
	
	if (lower)
		nnames = tolower(names(summaryT))
	else
		nnames =  names(summaryT)
	
	if (removeSpecialChars) {
		nnames = RemoveSpecialChars(nnames, replaceBy = replace)
	}
	
	#r = list(lbl = summaryT)
	#names(r) = label
	r = summaryT
	return(r)
}

replaceElementInFormula = function(formula,
																	 pattern  ,
																	 replace) {
	as.formula(gsub(
		pattern = pattern,
		replacement = replace,
		x = formula,
		fixed = TRUE
	))
}

RemoveSpecialChars = function(x,
															what = '[^0-9A-Za-z]',
															replaceBy = '_',
															message = FALSE) {
	what = gsub(
		pattern = ']',
		x = what,
		replacement = paste0(replaceBy, ']')
	)
	if (message)
		message0('pattern: ',
						 what ,
						 '; replaced by ',
						 replaceBy)
	r = gsub(what, replaceBy , x , ignore.case = TRUE)
	###
	if (any(nchar(r) < 1)) {
		RN = RandomRegardSeed(1, stringOutput = TRUE, round = 8)
		r[nchar(r) < 1] = paste('no_name',
														RN,
														sep = '_',
														collapse = '_')
	}
	return(r)
}

RandomRegardSeed = function(n = 1,
														max = 10 ^ 9,
														decimal = TRUE,
														round = 5,
														stringOutput = FALSE,
														what = '[^0-9A-Za-z]',
														replaceBy = '_') {
	r = runif(n, 1, as.numeric(Sys.time()) * 10000) %% max
	if (decimal)
		r = r %% 1
	r = round(r, digits = round)
	if (stringOutput) {
		r = RemoveSpecialChars(as.character(r), what = what, replaceBy = replaceBy)
	}
	return(r)
}


####
rndProce = function(procedure = NULL) {
	if (length(procedure) < 1   ||
			is.null(procedure)      ||
			#!(procedure %in% lop()) ||
			length(procedure) > 1) {
		message0 (
			'Error in inputing the procedure symbol. Current input: ',
			ifelse(is.null(procedure), 'NULL', procedure)              ,
			'. The default random effect, (1 |  Batch), is selected.'
		)
		return(reformulate  (' 1 |  Batch', response = NULL, intercept = TRUE))
	}
	##############################################################
	if (procedure == 'TYPICAL') {
		random = reformulate  (' 1 |  Batch', response = NULL, intercept = TRUE)
	} else if (procedure == 'ABR') {
		random = reformulate  (' 1 |  Batch', response = NULL, intercept = TRUE)
	} else if (procedure == 'ACS') {
		random = reformulate  (' 1 |  Batch', response = NULL, intercept = TRUE)
	} else if (procedure == 'ALZ') {
		random = reformulate  (' 1 |  Batch', response = NULL, intercept = TRUE)
	} else if (procedure == 'BLK') {
		random = reformulate  (' 1 |  Batch', response = NULL, intercept = TRUE)
	} else if (procedure == 'BWT') {
		random = reformulate  (' 1 |  Batch', response = NULL, intercept = TRUE)
	} else if (procedure == 'CAL') {
		random = reformulate  (' 1 |  Batch', response = NULL, intercept = TRUE)
	} else if (procedure == 'CBC') {
		random = reformulate  (' 1 |  Batch', response = NULL, intercept = TRUE)
	} else if (procedure == 'CHL') {
		random = reformulate  (' 1 |  Batch', response = NULL, intercept = TRUE)
	} else if (procedure == 'CSD') {
		random = reformulate  (' 1 |  Batch', response = NULL, intercept = TRUE)
	} else if (procedure == 'DXA') {
		random = reformulate  (' 1 |  Batch', response = NULL, intercept = TRUE)
	} else if (procedure == 'ECG') {
		random = reformulate  (' 1 |  Batch', response = NULL, intercept = TRUE)
	} else if (procedure == 'ECH') {
		random = reformulate  (' 1 |  Batch', response = NULL, intercept = TRUE)
	} else if (procedure == 'ELZ') {
		random = reformulate  (' 1 |  Batch', response = NULL, intercept = TRUE)
	} else if (procedure == 'EVL') {
		random = reformulate  (' 1 |  Batch', response = NULL, intercept = TRUE)
	} else if (procedure == 'EVM') {
		random = reformulate  (' 1 |  Batch', response = NULL, intercept = TRUE)
	} else if (procedure == 'EVO') {
		random = reformulate  (' 1 |  Batch', response = NULL, intercept = TRUE)
	} else if (procedure == 'EVP') {
		random = reformulate  (' 1 |  Batch', response = NULL, intercept = TRUE)
	} else if (procedure == 'EYE') {
		random = reformulate  (' 1 |  Batch', response = NULL, intercept = TRUE)
	} else if (procedure == 'FER') {
		random = reformulate  (' 1 |  Batch', response = NULL, intercept = TRUE)
	} else if (procedure == 'GEL') {
		random = reformulate  (' 1 |  Batch', response = NULL, intercept = TRUE)
	} else if (procedure == 'GEM') {
		random = reformulate  (' 1 |  Batch', response = NULL, intercept = TRUE)
	} else if (procedure == 'GEO') {
		random = reformulate  (' 1 |  Batch', response = NULL, intercept = TRUE)
	} else if (procedure == 'GEP') {
		random = reformulate  (' 1 |  Batch', response = NULL, intercept = TRUE)
	} else if (procedure == 'GPL') {
		random = reformulate  (' 1 |  Batch', response = NULL, intercept = TRUE)
	} else if (procedure == 'GPM') {
		random = reformulate  (' 1 |  Batch', response = NULL, intercept = TRUE)
	} else if (procedure == 'GPO') {
		random = reformulate  (' 1 |  Batch', response = NULL, intercept = TRUE)
	} else if (procedure == 'GRS') {
		random = reformulate  (' 1 |  Batch', response = NULL, intercept = TRUE)
	} else if (procedure == 'HEM') {
		random = reformulate  (' 1 |  Batch', response = NULL, intercept = TRUE)
	} else if (procedure == 'HIS') {
		random = reformulate  (' 1 |  Batch', response = NULL, intercept = TRUE)
	} else if (procedure == 'HWT') {
		random = reformulate  (' 1 |  Batch', response = NULL, intercept = TRUE)
	} else if (procedure == 'IMM') {
		random = reformulate  (' 1 |  Batch', response = NULL, intercept = TRUE)
	} else if (procedure == 'INS') {
		random = reformulate  (' 1 |  Batch', response = NULL, intercept = TRUE)
	} else if (procedure == 'IPG') {
		random = reformulate  (' 1 |  Batch', response = NULL, intercept = TRUE)
	} else if (procedure == 'OFD') {
		random = reformulate  (' 1 |  Batch', response = NULL, intercept = TRUE)
	} else if (procedure == 'PAT') {
		random = reformulate  (' 1 |  Batch', response = NULL, intercept = TRUE)
	} else if (procedure == 'VIA') {
		random = reformulate  (' 1 |  Batch', response = NULL, intercept = TRUE)
	} else if (procedure == 'XRY') {
		random = reformulate  (' 1 |  Batch', response = NULL, intercept = TRUE)
	} else{
		random = reformulate  (' 1 |  Batch'      , response = NULL, intercept = TRUE)
	}
	return(random)
}

shapiro.test0 = function(x, ...) {
	if (!is.null(x)                     &&
			is.numeric(x)                   &&
			length(x)           > 3         &&
			length(x)         < 5000        &&
			length(unique(na.omit(x))) > 3  &&
			!is.na(sd0(x, na.rm = TRUE))    &&
			sd0(x, na.rm = TRUE)             > 0) {
		r = list(
			'P-value'  = shapiro.test(x)$p.value   ,
			'N'        = length(x)                 ,
			'Unique n' = length(unique(na.omit(x))),
			'SD'       = sd0(x, na.rm = TRUE)      ,
			'Test'     = 'Shapiro'
		)
	} else{
		r = list(
			'P-value'  = NULL                                                               ,
			'N'        = length(x)                                                          ,
			'Unique n' = length(unique(na.omit(x)))                                         ,
			'SD'       = sd0(x, na.rm = TRUE)                                               ,
			'Test'     = 'Not possible (Possible cause:n < 3 or n > 5000 unique data points'
		)
	}
	return(r)
}


QuyalityTests = function(object,
												 levels    = c('Genotype', 'Sex', 'LifeStage'),
												 list      = TRUE,
												 noDataLab = 'No data',
												 sep       = '_',
												 collapse  = '_') {
	if (is.null(object) ||
			length (levels) < 1)
	{
		message0('\tSkipped . No model or the levels in the data for the quality test')
		return(NULL)
	}
	
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
					shapiro.test0(x)
				})
				if (!is.null(result)) {
					flst[[counter]] = result
					names(flst)[counter] = paste(cmb[, j], collapse = collapse, sep = sep)
					counter = counter + 1
				}
			}
		}
		flst$Overall =  shapiro.test0(x = r)
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

as.list0 = function(x, ...) {
	if (!is.null(x)) {
		return(as.list(x, ...))
	} else{
		return(NULL)
	}
}

AllEffSizes = function(object, depVariable, effOfInd, data) {
	if (length(effOfInd) < 1 ||
			effOfInd ==   1) {
		message0('\tSkipped . No variable found for the effect size ...')
		return(NULL)
	}
	lst       = flst = olst = NULL
	data      = NormaliseDataFrame(data)
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
			message0('\tLevel:', pasteComma(unlist(eff)))
			# 1. Main effect
			lstTmp = eff.size(
				object      = object,
				depVariable = depVariable,
				effOfInd    = eff,
				errorReturn = NULL,
				data        = data
			)
			if (!is.null(lstTmp)) {
				lst[[counter1]]      = lstTmp
				names(lst)[counter1] = eff
				counter1  = counter1  + 1
			}
			# 2. All subset interactions effect sizes
			if (length(cats) > 1) {
				for (j in 1:length(cats)) {
					cmbn = combn(x = cats, m = j)
					for (k in 1:ncol(cmbn)) {
						interact = interaction(data[, cmbn[, k] , drop = FALSE], sep = ' ', drop = TRUE)
						for (lvl in levels(interact)) {
							message0('\tLevel:', pasteComma(unlist(lvl)))
							olstTmp = eff.size(
								object      = object,
								data        = droplevels0(data[interact %in% lvl, , drop = FALSE]),
								depVariable = depVariable,
								errorReturn = NULL,
								effOfInd    = eff
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

extractAICc = function(fit, scale = FALSE, k = NULL, ...) {
	requireNamespace("AICcmodavg")
	# k=0 is just loglikelihood or -2*logLik(fit)
	if (k == 0)
		extractAIC(fit = fit,
							 scale = scale,
							 k = 0,
							 ...)
	else
		edf  = AICc(mod = fit, return.K = TRUE, ...)
	
	AIC  = AICc(mod = fit, return.K = FALSE, ...)
	l    = list(edf = edf, AIC = AIC)
	return(unlist(l))
}

stepAIC0 = function (object,
										 scope,
										 scale = 0,
										 direction = c("both", "backward",
										 							"forward"),
										 trace = 1,
										 keep = NULL,
										 steps = 1000,
										 use.start = FALSE,
										 k = 2,
										 ...)
{
	mydeviance <- function(x, ...) {
		dev <- deviance(x)
		if (!is.null(dev))
			dev
		else
			extractAICc(x, k = 0)[2L]
	}
	cut.string <- function(string) {
		if (length(string) > 1L)
			string[-1L] <- paste("\n", string[-1L], sep = "")
		string
	}
	re.arrange <- function(keep) {
		namr <- names(k1 <- keep[[1L]])
		namc <- names(keep)
		nc <- length(keep)
		nr <- length(k1)
		array(unlist(keep, recursive = FALSE), c(nr, nc), list(namr,
																													 namc))
	}
	step.results <- function(models, fit, object, usingCp = FALSE) {
		change <- sapply(models, "[[", "change")
		rd <- sapply(models, "[[", "deviance")
		dd <- c(NA, abs(diff(rd)))
		rdf <- sapply(models, "[[", "df.resid")
		ddf <- c(NA, abs(diff(rdf)))
		AIC <- sapply(models, "[[", "AIC")
		heading <- c(
			"Stepwise Model Path \nAnalysis of Deviance Table",
			"\nInitial Model:",
			deparse(formula(object)),
			"\nFinal Model:",
			deparse(formula(fit)),
			"\n"
		)
		aod <- if (usingCp)
			data.frame(
				Step = change,
				Df = ddf,
				Deviance = dd,
				`Resid. Df` = rdf,
				`Resid. Dev` = rd,
				Cp = AIC,
				check.names = FALSE
			)
		else
			data.frame(
				Step = change,
				Df = ddf,
				Deviance = dd,
				`Resid. Df` = rdf,
				`Resid. Dev` = rd,
				AIC = AIC,
				check.names = FALSE
			)
		attr(aod, "heading") <- heading
		class(aod) <- c("Anova", "data.frame")
		fit$anova <- aod
		fit
	}
	Terms <- terms(object)
	object$formula <- Terms
	if (inherits(object, "lme"))
		object$call$fixed <- Terms
	else if (inherits(object, "gls"))
		object$call$model <- Terms
	else
		object$call$formula <- Terms
	if (use.start)
		warning("'use.start' cannot be used with R's version of 'glm'")
	md <- missing(direction)
	direction <- match.arg(direction)
	backward <- direction == "both" | direction == "backward"
	forward <- direction == "both" | direction == "forward"
	if (missing(scope)) {
		fdrop <- numeric()
		fadd <- attr(Terms, "factors")
		if (md)
			forward <- FALSE
	}
	else {
		if (is.list(scope)) {
			fdrop <- if (!is.null(fdrop <- scope$lower))
				attr(terms(update.formula(object, fdrop)), "factors")
			else
				numeric()
			fadd <- if (!is.null(fadd <- scope$upper))
				attr(terms(update.formula(object, fadd)), "factors")
		}
		else {
			fadd <- if (!is.null(fadd <- scope))
				attr(terms(update.formula(object, scope)), "factors")
			fdrop <- numeric()
		}
	}
	models <- vector("list", steps)
	if (!is.null(keep))
		keep.list <- vector("list", steps)
	n <- nobs(object, use.fallback = TRUE)
	fit <- object
	bAIC <- extractAICc(fit, scale, k = k, ...)
	edf <- bAIC[1L]
	bAIC <- bAIC[2L]
	if (is.na(bAIC))
		stop("AIC is not defined for this model, so 'stepAIC' cannot proceed")
	if (bAIC == -Inf)
		stop("AIC is -infinity for this model, so 'stepAIC' cannot proceed")
	nm <- 1
	Terms <- terms(fit)
	if (trace) {
		cat("Start:  AICc=",
				format(round(bAIC, 2)),
				"\n",
				cut.string(deparse(formula(fit))),
				"\n\n",
				sep = "")
		utils::flush.console()
	}
	models[[nm]] <- list(
		deviance = mydeviance(fit),
		df.resid = n -
			edf,
		change = "",
		AIC = bAIC
	)
	if (!is.null(keep))
		keep.list[[nm]] <- keep(fit, bAIC)
	usingCp <- FALSE
	while (steps > 0) {
		steps <- steps - 1
		AIC <- bAIC
		ffac <- attr(Terms, "factors")
		if (!is.null(sp <-
								 attr(Terms, "specials")) && !is.null(st <- sp$strata))
			ffac <- ffac[-st, ]
		scope <- factor.scope(ffac, list(add = fadd, drop = fdrop))
		aod <- NULL
		change <- NULL
		if (backward && length(scope$drop)) {
			aod <- dropterm0(
				fit,
				scope$drop,
				scale = scale,
				trace = max(0,
										trace - 1),
				k = k,
				...
			)
			rn <- row.names(aod)
			row.names(aod) <- c(rn[1L], paste("-", rn[-1L], sep = " "))
			if (any(aod$Df == 0, na.rm = TRUE)) {
				zdf <- aod$Df == 0 & !is.na(aod$Df)
				nc <- match(c("Cp", "AIC"), names(aod))
				nc <- nc[!is.na(nc)][1L]
				ch <- abs(aod[zdf, nc] - aod[1, nc]) > 0.01
				if (any(is.finite(ch) & ch)) {
					warning("0 df terms are changing AIC")
					zdf <- zdf[!ch]
				}
				if (length(zdf) > 0L)
					change <- rev(rownames(aod)[zdf])[1L]
			}
		}
		if (is.null(change)) {
			if (forward && length(scope$add)) {
				aodf <- addterm(
					fit,
					scope$add,
					scale = scale,
					trace = max(0, trace - 1),
					k = k,
					...
				)
				rn <- row.names(aodf)
				row.names(aodf) <- c(rn[1L], paste("+", rn[-1L],
																					 sep = " "))
				aod <- if (is.null(aod))
					aodf
				else
					rbind(aod, aodf[-1, , drop = FALSE])
			}
			attr(aod, "heading") <- NULL
			if (is.null(aod) || ncol(aod) == 0)
				break
			nzdf <- if (!is.null(aod$Df))
				aod$Df != 0 | is.na(aod$Df)
			aod <- aod[nzdf, ]
			if (is.null(aod) || ncol(aod) == 0)
				break
			nc <- match(c("Cp", "AIC"), names(aod))
			nc <- nc[!is.na(nc)][1L]
			o <- order(aod[, nc])
			if (trace) {
				print(aod[o, ])
				utils::flush.console()
			}
			if (o[1L] == 1)
				break
			change <- rownames(aod)[o[1L]]
		}
		usingCp <- match("Cp", names(aod), 0) > 0
		fit <- update(fit, paste("~ .", change), evaluate = FALSE)
		fit <- eval.parent(fit)
		nnew <- nobs(fit, use.fallback = TRUE)
		if (all(is.finite(c(n, nnew))) && nnew != n)
			message0("\t Warning! number of rows in use has changed: remove missing values ...")
		Terms <- terms(fit)
		bAIC <- extractAICc(fit, scale, k = k, ...)
		edf <- bAIC[1L]
		bAIC <- bAIC[2L]
		if (trace) {
			cat("\nStep:  AIC=",
					format(round(bAIC, 2)),
					"\n",
					cut.string(deparse(formula(fit))),
					"\n\n",
					sep = "")
			utils::flush.console()
		}
		if (bAIC >= AIC + 1e-07)
			break
		nm <- nm + 1
		models[[nm]] <- list(
			deviance = mydeviance(fit),
			df.resid = n -
				edf,
			change = change,
			AIC = bAIC
		)
		if (!is.null(keep))
			keep.list[[nm]] <- keep(fit, bAIC)
	}
	if (!is.null(keep))
		fit$keep <- re.arrange(keep.list[seq(nm)])
	step.results(models = models[seq(nm)], fit, object, usingCp)
}

dropterm0 =
	function (object,
						scope,
						scale = 0,
						test = c("none", "Chisq"),
						k = 2,
						sorted = FALSE,
						trace = FALSE,
						...)
	{
		tl <- attr(terms(object), "term.labels")
		if (missing(scope))
			scope <- drop.scope(object)
		else {
			if (!is.character(scope))
				scope <- attr(terms(update.formula(object, scope)),
											"term.labels")
			if (!all(match(scope, tl, 0L)))
				stop("scope is not a subset of term labels")
		}
		ns <- length(scope)
		ans <-
			matrix(nrow = ns + 1L,
						 ncol = 2L,
						 dimnames = list(c("<none>",
						 									scope), c("df", "AIC")))
		ans[1, ] <- extractAICc(object, scale, k = k, ...)
		n0 <- nobs(object, use.fallback = TRUE)
		env <- environment(formula(object))
		for (i in seq_len(ns)) {
			tt <- scope[i]
			if (trace) {
				message(gettextf("trying - %s", tt), domain = NA)
				utils::flush.console()
			}
			nfit <- update(object, as.formula(paste("~ . -", tt)),
										 evaluate = FALSE)
			nfit <- eval(nfit, envir = env)
			ans[i + 1, ] <- extractAICc(nfit, scale, k = k, ...)
			nnew <- nobs(nfit, use.fallback = TRUE)
			if (all(is.finite(c(n0, nnew))) && nnew != n0)
				message0("\t Warning! number of rows in use has changed: remove missing values ...")
		}
		dfs <- ans[1L, 1L] - ans[, 1L]
		dfs[1L] <- NA
		aod <- data.frame(Df = dfs, AIC = ans[, 2])
		o <- if (sorted)
			order(aod$AIC)
		else
			seq_along(aod$AIC)
		test <- match.arg(test)
		if (test == "Chisq") {
			dev <- ans[, 2L] - k * ans[, 1L]
			dev <- dev - dev[1L]
			dev[1L] <- NA
			nas <- !is.na(dev)
			P <- dev
			P[nas] <- safe_pchisq0(dev[nas], dfs[nas], lower.tail = FALSE)
			aod[, c("LRT", "Pr(Chi)")] <- list(dev, P)
		}
		aod <- aod[o, ]
		head <-
			c("Single term deletions", "\nModel:", deparse(formula(object)))
		if (scale > 0)
			head <- c(head, paste("\nscale: ", format(scale), "\n"))
		class(aod) <- c("anova", "data.frame")
		attr(aod, "heading") <- head
		aod
	}

safe_pchisq0 <- function(q, df, ...)
{
	df[df <= 0] <- NA
	pchisq(q = q, df = df, ...)
}

lowHighList = function(x, y, ...) {
	if (!is.null(x) || !is.null(y))
		r = list('Low' = x, 'High' = y , ...)
	else
		r = list('Low' = x, 'High' = y)
	return(r)
}

factorise = function(dataset, variables) {
	vars = variables[variables %in% names(dataset)]
	if (length(vars) > 0) {
		for (v in vars) {
			dataset[, v] = as.factor(dataset[, v])
		}
	}
	return(dataset)
}

PhenListAgeingRelabling = function(dataset, col, l1, rel1, l2, rel2) {
	if (col %in% names(dataset)) {
		dataset[, col] = as.factor(dataset[, col])
		levels(dataset[, col]) = capitalise(levels(dataset[, col]))
		if (!is.null(l1))
			levels(dataset[, col])[levels(dataset[, col]) %in% l1] =
			rel1
		if (!is.null(l2))
			levels(dataset[, col])[levels(dataset[, col]) %in% l2] =
			rel2
		lbls = c(rel1, rel2)
		if (!is.null(lbls) && any(!levels(dataset[, col]) %in% lbls)) {
			message0(
				'There are some unused levels in `',
				col,
				'` that will be removed. \n\t Levels: ',
				pasteComma(levels(dataset[, col]))
			)
			RowIndFromCol = dataset[, col] %in% lbls
			if (sum(RowIndFromCol) < 1) {
				message0('\t  Preprocessing variable leads to an empty column \n\t   The variable renamed to `',
								 paste0(col, '_labels`'))
				names(dataset)[names(dataset) %in% col] = paste0(col, '_labels')
			} else{
				dataset = droplevels(subset(dataset, RowIndFromCol))
			}
		}
	}
	return(dataset)
}

checkSummary = function(dataset, var, ...) {
	lvls = NULL
	if (is.null(dataset) ||
			is.null(var)     ||
			length (na.omit(dataset[, var])) < 1)
		return(lvls)
	
	if (is.factor(dataset[, var]) || is.character(dataset[, var]))
		lvls = paste0(
			'\t Levels (Total levels = ',
			nlevels(as.factor(dataset[, var])),
			'): \n\t  ',
			pasteComma(levels(as.factor(dataset[, var])), width = 250, sep = '\n\t  ')
		)
	else
		lvls = paste0(
			'\t Summary:\n\t  mean = ',
			mean(dataset[, var], na.rm = TRUE),
			'\n\t  sd   = ',
			sd0(dataset[, var], na.rm = TRUE)
		)
	message0(lvls)
	return(invisible(lvls))
}

checkPhenlistColumns = function(dataset, vars) {
	allExist = NULL
	if (length(vars)) {
		for (v in vars) {
			r = v %in% names(dataset)
			message0('checking whether variable `',
							 v,
							 '` exists in the data ... \n\tResult = ',
							 r)
			if (r)
				checkSummary(dataset = dataset, var = v)
			allExist = c(allExist , r)
		}
	} else{
		allExist = FALSE
	}
	return(allExist)
}


fIsBecome = function(dataset,
										 is = 'Assay.Date',
										 rename = 'Batch',
										 ifexist = NULL) {
	if (is.null(is) || is.null(rename))
		return(dataset)
	
	
	for (iss in is) {
		if (iss %in% colnames(dataset) &&
				ifelse(is.null(ifexist), TRUE, !iss %in% ifexist)) {
			####
			if (!is.null(ifexist)               &&
					ifexist  %in% colnames(dataset) &&
					!ifexist %in% is) {
				colnames(dataset)[colnames(dataset) %in% rename] = paste('Original', ifexist, sep = '.')
			}
			###
			names(dataset)[names(dataset) %in% iss] =
				rename
			message0('Variable `',
							 iss,
							 '` renamed to `',
							 rename, '`')
			break
		}
	}
	return(dataset)
}

is.df.empty = function(x) {
	if (is.null(x) || any(dim(x) < 1))
		return(TRUE)
	else
		return(FALSE)
}

CleanEmptyRecords = function(x, vars) {
	if(is.df.empty(x))
		return(NULL)
	vars1 = vars[vars %in% names(x)]
	if (length(vars1)) {
		x = x[all(!is.na(x[, vars1])) && all(x[, vars1] != ""), ]
	}
	return(x)
}


fastBarnardextest  = function(x,
															tail = 2 ,
															prob = seq(10 ^ -10, 1 - 10 ^ -10, length.out = 101),
															plot = FALSE) {
	if (is.null(x)  ||
			!(is.table(x) ||
				is.matrix(x)) ||
			any(dim(x) != 2))
		stop('The input must be a 2 by 2 table/matrix')
	
	fprob = function(i, j, c1, c2) {
		n  = c1 + c2
		pa = i / c1
		pb = j / c2
		px = (i + j) / n
		if (px == 0 || pa == pb) {
			return(0)
		} else
			return((pa - pb) / sqrt(px * (1 - px) * ((1 / c1) + (1 / c2))))
	}
	c  = colSums(x)
	r  = rowSums(x)
	c1 = c[1]
	c2 = c[2]
	n = sum(c)
	pao = x[1, 1] / c1
	pbo = x[1, 2] / c2
	pxo = r[1] / n
	TXO = abs(pao - pbo) / sqrt(pxo * (1 - pxo) * (1 / c1 + 1 / c2))
	cbn = matrix(c(rep(0:c1, each = c2 + 1), rep(0:c2, c1 + 1)), nrow = 2, byrow = TRUE)
	
	n1     = lfactorial(c1)
	n2     = lfactorial(c2)
	lprob  = log(prob)
	clprob = log(1 - prob)
	Fact   = lfactorial(0:max(c1, c2, na.rm = TRUE)[1])
	###################################################
	T = sapply(1:ncol(cbn), function(col) {
		i   = cbn[1, col]
		j   = cbn[2, col]
		s = n1 + n2 +
			(i + j) * lprob +
			(n - (i + j)) * clprob - sum(Fact[c(i + 1,
																					j + 1,
																					c1 - i + 1,
																					c2 - j + 1)])
		t = fprob(i = i,
							j = j,
							c1 = c1,
							c2 = c2)
		return(c(t = t, s = exp(s)))
	})
	r = t(cbind(P = apply(T[-1,], 1, function(x) {
		sum(x[T[1,] >= TXO])
	}), prob = prob))
	###################################################
	Nuisance.parameter = unlist(r[2, ][which.max(r[1, ])])
	p.value            = unlist(r[1, ][which.max(r[1, ])])
	if (plot) {
		plot(
			r[2, ],
			r[1, ],
			type = "l",
			main = "Barnard's exact P-value",
			xlab = "Nuisance parameter",
			ylab = "P.value"
		)
		abline(v = Nuisance.parameter, col = 2, lwd = 2)
	}
	return(
		list(
			p.value            = unname(min(tail * p.value, 1)),
			Nuisance.parameter = unname(Nuisance.parameter),
			Wald.Statistic     = unname(TXO)               ,
			tail               = tail                      ,
			seq                = r
		)
	)
}

