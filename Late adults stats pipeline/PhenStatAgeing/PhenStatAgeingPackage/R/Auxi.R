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

Matrix2List = function(x, ...) {
	if (is.null(x))
		return(NULL)
	
	if (length(x) == 1) {
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
											replaceNull = TRUE,
											truncate = TRUE   ,
											width = 100) {
	if (replaceNull)
		r = paste(replaceNull(list(...), replaceBy = 'NULL'),
							sep = ', ',
							collapse = ', ')
	else
		r = paste(..., sep = ', ', collapse = ', ')
	if (truncate)
		r = truncate_text(r,width)
}

truncate_text = function(x,width){
	ifelse(nchar(x) > width, paste0(strtrim(x, width), '...'), x)
}

pasteUnderscore = function(...) {
	paste(..., sep = '_', collapse = '_')
}

checkModelTermsInData = function(formula,
																 data,
																 responseIsTheFirst  = TRUE,
																 pattern = '[.~+-()]') {
	formula = as.formula(formula)
	vars    = all.vars(formula, functions = FALSE)
	vars    = vars[!grepl(pattern = pattern,
												x = vars,
												fixed = TRUE)]
	if (responseIsTheFirst) {
		if (!(vars[1] %in% names(data))) {
			message0('Response does not exist in the data!')
			stop(
				'Response has not been specified properly. Please check that the response exists in the data'
			)
		}
	}
	In      = vars %in% names(data)
	if (any(!In)) {
		message0(
			'Some terms in the model are not included in the data. See: \n\t ',
			pasteComma(vars[!In], replaceNull = FALSE),
			'\n\t Initial  model: ',printformula(formula)
		)
		ft  = vars [!In]
		formula = update.formula(formula,
														 reformulate0(
														 	termlabels = c('.', ft, paste0(ft, ':.')),
														 	response   = NULL,
														 	intercept  = TRUE,
														 	sep        = '-'
														 ))
		message0('\t Polished model: ',printformula(formula))
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
	if ( sup.messages && sup.warnings) {
		suppressMessages(suppressWarnings(exp))
	} else if ( sup.messages && !sup.warnings) {
		suppressMessages(exp)
	} else if ( !sup.messages && sup.warnings) {
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
	inF = fixed
	fixed = ComplementaryFeasibleTermsInContFormula(formula = fixed, data = data)
	
	
	if (debug)
		message0(
			'Initial  model: ',
			printformula(inF),
			'\n',
			'Polished model: ',
			printformula(fixed),
			' [any removal: ',
			!identical(fixed, inF),
			']'
		)
	#return(list(correctd = fixed, initial = inF))
	return(fixed)
}



FeasibleTermsInContFormula = function(formula, data) {
	Allvars = all.vars(formula)[all.vars(formula) %in% names(data)]
	isCat = !sapply(data[, Allvars, drop = FALSE], is.numeric)
	vars  = Allvars[isCat]
	lvars = length(vars)#min(length(vars), sapply(strsplit(formulaTerms(formula = formula), split = ':'), length), na.rm = TRUE)
	names = r = NULL
	if(getResponseFromFormula(formula = formula) %in% vars){
		message0('\tResponse is included in the checks ....')
	}
	if (lvars > 0) {
		for (i in 1:lvars) {
			message0('\t',i,
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
		message0('Variables that being found in data: ', pasteComma(newNames))
	return(newNames)
}

ComplementaryFeasibleTermsInContFormula = function(formula, data) {
	message0(
		'Checking for the feasibility of terms and interactions. Formula:\n\t ',
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
				'The following term(s) removed because there is either "no data" or "no data in the interactions":\n\t **not all terms necessarily in the initial model \n\t ',
				pasteComma(fbm[fbm$min.freq <= 0, c('names')], replaceNull = FALSE)
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
	org.data = data
	new.data = data[complete.cases(data[, all.vars(formula)]), ]
	missings = ifelse(all(dim(org.data) == dim(new.data)),	0, dim(org.data)[1] -
											dim(new.data)[1])
	if (missings)
		message0(
			'The data (variable(s) = ',
			pasteComma(all.vars(formula)),
			') contain ',
			missings,
			' missing(s)...\n\tMissing data removed.'
		)
	return(invisible(
		list(
			org.data = org.data,
			new.data = new.data,
			missings = missings
		)
	))
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
																sep = '.') {
	if (!(
		!is.null(data)                                           &&
		(!is.null(variable) || !is.null(mainEffsOnlyWhenIndivi)) &&
		!is.null(model)                                          &&
		!is.null(FUN(data[, depVar]))                            &&
		FUN(data[, depVar]) != 0
	))
		return(NULL)
	####
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
				r = coefs
				names(r) = order0(levels(data[, x]))
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
				value            = out                                      ,
				percentageChange = Matrix2List(out / ran * 100, sep = sep)  ,
				variable         = mainEffsOnlyWhenIndivi                   ,
				type             = 'interaction coefficients '
			)
		)
	}
}

extractCoefOfInterest = function(coefs, main = 'Sex',data) {
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
	r = coefs[grepl(pattern = pasteBracket( unlist(lvls),left = '',right = '',col = '|'), x = Cnames)]
	return(r)
}

pasteBracket = function(..., replaceNull = TRUE,col='|',right=']:',left='[') {
	if (replaceNull)
		paste(left,replaceNull(list(...), replaceBy = 'NULL'),right,
					sep = '',
					collapse = col)
	else
		paste(left,...,right, sep = ' ', collapse = col)
}

eff.size = function(object,
										data        = NULL        ,
										depVariable = 'data_point',
										effOfInd    = 'Genotype'  ,
										errorReturn = NULL        ,
										debug       = FALSE) {
	if (all(is.null(data)))
		data = getData(object)
	f    = reformulate(termlabels = effOfInd, depVariable)
	agr  = aggregate(f, data = data, FUN = mean)
	if (debug) {
		cat('\n\n')
		print(agr)
		print(dim(agr))
	}
	if (any(dim(agr) < 2)) {
		message0(
			' \t\tEffect size estimation: No variation or less than two levels in ',
			pasteComma(effOfInd, collapse = ',')
		)
		return(errorReturn)
	}
	
	NModel       =
		tryCatch(
			expr = update(object,  reformulate(termlabels = effOfInd,response = NULL,intercept = TRUE), data = data),
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
			model = NModel,
			data = data,
			variable = effOfInd,
			depVar = depVariable
		)
		if (sum(CoefEffSizes)) {
			# For continues covariates it is the coefficient
			efSi         = list(value            = as.list(coef(NModel))[[effOfInd]][1],
													percentageChange = PerChange,
													variable         = effOfInd,
													type             = 'coefficient')
		} else{
			# For categorical covariates it is the mean difference
			MDiff        = max(dist(agr[, depVariable, drop = FALSE]), na.rm = TRUE)
			r            = resid(NModel)
			efSi         = list(value            = ifelse(sd(r) > 0, abs(MDiff) / sd(r), NA),
													percentageChange = PerChange,
													variable         = effOfInd,
													type             = 'Mean difference')
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
		message0('Used model: ', FUN)
	fArgs = names(list) %in% formalArgs(args(FUN))
	l     = list[names(list)[fArgs]]
	return(l)
}

ModelChecks = function(fixed, data, checks = c(0, 0, 0)) {
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
	return(all(term %in% all.vars(formula(model))))
}

SplitEffect = 	function(finalformula,
												fullModelFormula,
												F.Model,
												data,
												depVariable,
												mandatoryVar = 'Genotype') {
	Allargs  = all.vars(fullModelFormula)[!all.vars(fullModelFormula) %in% c(depVariable, mandatoryVar)]
	isCat    = !sapply(data[, Allargs, drop = FALSE], is.numeric)
	args     = Allargs[isCat]
	argsCon  = if (length(Allargs[!isCat]) > 0) {
		Allargs[!isCat]
	} else{
		#	NULL
		1
	}
	largs    = length(args)
	l        = NULL
	names    = c()
	counter = 1
	if (largs > 0) {
		for (i in 1:largs) {
			argComb = combn(args, i, simplify = TRUE)
			for (j in 1:ncol(argComb)) {
				arg = as.vector(argComb[, j])
				if (arg %in% names(data) &&
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
							arg,
							paste(
								mandatoryVar,
								arg,
								collapse = ':',
								sep = ':'
							),
							argsCon
						),
						intercept = TRUE
					)
					message0('Check the split model:\n\t', printformula(newModel))
					# we can safely remove the line below!
					# newModel = ModelChecks(fixed = newModel,
					# 											 data = data,
					# 											 checks = c(1, 1, 0))
					
					l0    = tryCatch(
						update(F.Model,
									 newModel),
						error = function(e) {
							message0(e,breakLine=FALSE)
							return(NULL)
						} ,
						warning = function(w) {
							message0(w,breakLine=FALSE)
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
						l0$MainEffect = arg
						l[[counter]]   = l0
						names[counter] = paste(
							paste0(mandatoryVar, collapse = '_'),
							paste0(arg, collapse = '.'),
							collapse = '_',
							sep = '_'
						)
						counter = counter + 1
					}
				}
			}
		}
		if (length(names) > 0) {
			message0('SplitEffects. Output names: ',
							 paste0(names, collapse = ', '))
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

printformula = function(formula) {
	if (!is.null(formula)) {
		r = paste01(format(formula, trim = TRUE, width = 0), collapse = '')
	} else{
		r = NULL
	}
	return(r)
}
# Categorical effect size
# https://github.com/mpi2/stats_working_group/raw/master/PhenStatUserGuide/PhenStatUsersGuide.pdf p122
cat.eff.size = function(xtb,varName = NULL) {
	if (any(dim(xtb) < 1)) {
		r = NULL
	} else{
		r = max(apply(prop.table(xtb, margin = 2), 1, function(x) {
			max(dist(x, method = 'maximum', diag = TRUE), na.rm = TRUE)
		}), na.rm = TRUE)
	}
	
	out = list(
		value = r                                         ,
		percentageChange = NULL                           ,
		variable         = ifelse(!is.null(varName)       ,
															varName                 , 
															'Variable does not exist'),
		type = 'Proportion change'
	)
	return(out)
}
# Test engine
ctest = function(x,
								 formula,
								 asset = NULL,
								 rep = 1500,
								 ...) {
	message = c()
	xtb = xtabs(
		formula = formula,
		data = x,
		drop.unused.levels = TRUE,
		na.action = 'na.omit'
	)
	if (!is.null(asset)) {
		if (asset <= dim(xtb)[3]) {
			xtb = xtb[, , asset]
		} else{
			message = c(message, 'There are some empty levels in data')
			return(
				list(
					result     = NULL,
					effectSize = NULL,
					note       = message,
					table      = xtb,
					input      = x  ,
					formula    = formula
				)
			)
		}
	}
	if (length(dim0(xtb)) < 2                      ||
			sum(margin.table(xtb, margin = 2) > 0) < 2 ||
			sum(margin.table(xtb, margin = 1) > 0) < 2) {
		message = c(message,
								'No variation in al least two levels of one or more categories')
		return(
			list(
				result = NULL,
				effectSize = NULL,
				note    = message,
				table   = xtb,
				input   = x,
				formula = formula
			)
		)
	}
	
	if (any(dim(xtb) < 2)) {
		r = list(p.value = 1)
		effect = 1
	} else{
		r = tryCatch(
			fisher.test(
				xtb,
				simulate.p.value = rep > 0,
				conf.int = FALSE,
				B = rep,
				...
			),
			error = function(e) {
				message0(e,breakLine=FALSE)
				return(NULL)
			} ,
			warning = function(w) {
				message0(w,breakLine=FALSE)
				return(NULL)
			}
		)
		effect     = cat.eff.size(xtb, varName = pasteComma(all.vars0(formula)[-1]))
		r$formula   = formula
		r$table     = xtb
	}
	return(
		list(
			result      = r      ,
			effectSize  = effect ,
			note        = message,
			table       = xtb    ,
			input       = 'Only have values when the function fails',
			formula     = formula
		)
	)
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
		all.vars(formula)[1]
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
								NULL
	)
}

UnlistCall = function(x){
	as(unlist(x),'character')
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
				all.vars(model)[1]
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
	if(is.null(x))
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
	newFormula = update.formula(old = formula, new = paste0('~', all.vars(formula)[1], '+.'))
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

RRNewObjectAndFormula = function(object,
																 RRprop,
																 formula,
																 labels = NULL,
																 depVarPrefix = NULL) {
	allTerms         = all.vars0(formula)
	newobject        = object
	RRcutObject      = RRCut(
		object = object          ,
		prob   = RRprop          ,
		depVariable = allTerms[1],
		labels = labels,
		depVarPrefix = depVarPrefix
	)
	newobject@datasetPL = RRcutObject$discObject
	newFormula          = replaceElementInFormula(
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
			depVariable    = allTerms[1]
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
		if(i == maxtry)
			message0('No solusion found for the RR_prop. You may want to revise the RR_prop?')
		xx = jitter(x = x,
								amount = amount,
								factor = factor)
		if (min(xx,na.rm = TRUE) > lower && max(xx,na.rm = TRUE) < upper)
			break
	}
	return(xx)
}


RRCut = function(object                     ,
								 prob         = .95         ,
								 depVariable  = 'data_point',
								 labels       = NULL,
								 depVarPrefix = NULL) {
	if (prob == .5) {
		message0('"prob" must be different from 0.5')
		return(NULL)
	}
	# Preparation ...
	object@datasetPL$data_point_discretised  = NA
	controls = subset(object@datasetPL, object@datasetPL$Genotype %in% object@refGenotype)
	mutants  = subset(object@datasetPL, object@datasetPL$Genotype %in% object@testGenotype)
	prb      = unique(c(0,  prob, 1))
	qntl     = quantile(x = controls[, depVariable],
											probs = prb,
											na.rm = TRUE)
	if (sum(duplicated(qntl))) {
		message0(
			'duplicates in quantiles detected, then small (dependes ont he data precision) jitter will be added to quantiles. Qauntiles: ',
			pasteComma(qntl, replaceNull = FALSE)
		)
		JitterPrecision = 1 + decimalplaces(min(controls[, depVariable], na.rm = TRUE))
		message0('Jitter precision (decimal) = ', JitterPrecision)
		qntl[duplicated(qntl)] = jitter0(
			x = qntl[duplicated(qntl)],
			amount = 10 ^	-JitterPrecision,
			upper = 1,
			lower = .5
		)
		qntl = sort(qntl)
	}
	message0('quantile(s) for cutting the data '      ,
					 '[probs = '                                 ,
					 pasteComma(round(prb,3))                 ,
					 ']: '                                    ,
					 pasteComma(qntl, replaceNull = FALSE))
	controls$data_point_discretised  = cut(
		x = controls[, depVariable],
		breaks = qntl,
		labels = if (is.null(labels))
			FALSE
		else
			labels,
		include.lowest = TRUE
	)
	mutants$data_point_discretised   = cut(
		x = mutants [, depVariable],
		breaks = qntl,
		labels = if (is.null(labels))
			FALSE
		else
			labels,
		include.lowest = TRUE
	)
	newObj                          = rbind(mutants, controls)
	newdepVariable                  = pasteUnderscore(depVarPrefix, depVariable, 'discretised')
	names(newObj)[names(newObj) %in% 'data_point_discretised'] = newdepVariable
	newObj[, newdepVariable] = as.factor(newObj[, newdepVariable])
	message0('a new column is added to the data object: ', newdepVariable)
	return(
		list(
			object = object,
			discObject = newObj,
			newdepVariable = newdepVariable ,
			depVariable = depVariable
		)
	)
}


RRextra = function(object,
									 prob = .95,
									 depVariable = 'data_point') {
	if (prob <= .5) {
		message0('"prob" must be greater than 0.5')
		return(NULL)
	}
	# Preparation ...
	controls = subset(object@datasetPL, object@datasetPL$Genotype %in% object@refGenotype)
	mutants  = subset(object@datasetPL, object@datasetPL$Genotype %in% object@testGenotype)
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
	if (is.null(model) || is.null(data) || is.null(term))
		return(FALSE)
	# if (!is.null(model$correctd))
	# 	model      = model$correctd
	r = termInTheModel(model = model ,
										 term = term,
										 message = FALSE) &&
		colLevelsSimple(data, all.vars0(model)[1]) >		threshold
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

message0 = function(...,
										breakLine = TRUE,
										capitalise = TRUE) {
	x = paste0(..., collapse = '')
	if (breakLine)
		nmessage = unlist(strsplit(x = x, split = '\n'))
	else
		nmessage = x
	if (capitalise)
		nmessage = capitalise(nmessage)
	message(paste(Sys.time(), nmessage, sep = '. ', collapse = '\n'))
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


all.vars0 = function(x, ...) {
	if (is.null(x))
		return(NULL)
	fif = all.vars(formula(x), ...)
	if (length(fif) > 0) {
		return(fif)
	} else{
		return(NULL)
	}
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
	cat    = all.vars(formula)[!sapply(data[, all.vars(formula)], is.numeric)]
	if (length(cat)) {
		FactsThatMustBeRemoved = cat[lapply(data[, cat, drop = FALSE], function(x) {
			length(unique(na.omit(x)))
		}) <= 1]
		if (length(FactsThatMustBeRemoved)) {
			message0(
				'The following terms from the model are removed because they only contain one level: ',
				pasteComma(FactsThatMustBeRemoved, replaceNull = FALSE)
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
																		 variable  = 'Genotype'                        ,
																		 anova     = TRUE                              ,
																		 what      = c('Pr(>|z|)', 'Pr(>Chi)', 'p-value'),
																		 debug     = TRUE) {
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
	return(as.vector(unlist(mSumFiltered)))
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
	depVar  = all.vars(formula)[1]
	# do not move me
	if (any(dim(x) == 0))
		return('empty dataset')
	
	cat    = all.vars(formula)[!sapply(x[, all.vars(formula)], is.numeric)]
	lvls   = interaction(x[, cat], sep = sep, drop = drop)
	isNumeric = is.numeric(x[, depVar])
	
	summaryT   = as.list(tapply(x[, depVar], INDEX = lvls, function(xx) {
		if (isNumeric) {
			c  = ifelse(length(na.omit(xx)) > 0, length(na.omit(xx)), 0)
			m  = ifelse(length(na.omit(xx)) > 0, mean(xx, na.rm = TRUE), NA)
			sd = ifelse(length(na.omit(xx)) > 1, sd(xx, na.rm = TRUE)  , NA)
			r = list(
				count = c                       ,
				mean = m                        ,
				sd = sd                         ,
				normality_test = ifelse(
					length(xx)           > 3    &&
						length(unique(xx)) > 3    &&
						length(xx)         < 5000 &&
						var(xx)            != 0,
					shapiro.test(xx)$p.value,
					'Not possible(Possible causes: <3 or >5000 unique data points'
				)
			)
		} else{
			c = ifelse(length(na.omit(xx)) > 0, length(na.omit(xx)), 0)
			r = list(count = c,
							 mean = NULL,
							 sd = NULL)
		}
		return(r)
	}, default = -999.991233210123))
	##
	fTmp = function(isNum) {
		if (isNum) {
			r = list(count = 0,
							 mean = NA,
							 sd = NA)
		} else{
			r = list(count = 0)
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
rndProce = function(procedure = NULL ) {
	if (length(procedure) < 1   ||
			is.null(procedure)      ||
			#!(procedure %in% lop()) ||
			length(procedure) > 1) {
		message0 (
			'Error in inputing the procedure symbol. Current input: ',
			ifelse(is.null(procedure),'NULL',procedure)              ,
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
	}else{
		random = reformulate  (' 1 |  Batch'      , response = NULL, intercept = TRUE)
	}
	return(random)
}

QuyalityTests = function(object,
												 levels    = c('Genotype', 'Sex', 'LifeStage'),
												 list      = TRUE,
												 noDataLab = 'No data',
												 sep       = '_',
												 collapse  = '_') {
	if (is.null(object))
		return(NULL)
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
					normality_test = if (!is.null(x) && !is.na(x) &&
															 length(x)           > 3  &&
															 length(unique(x)) > 3    &&
															 length(x)         < 5000 &&
															 var(x)            != 0) {
						shapiro.test(x)$p.value
					} else{
						'Not possible(Possible causes: <3 or >5000 unique data points'
					}
				})
				if (!is.null(result)) {
					flst[[counter]] = result
					names(flst)[counter] = paste(cmb[, j], collapse = collapse, sep = sep)
					counter = counter + 1
				}
			}
		}
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
	lst       = flst = olst = NULL
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
				object = object,
				depVariable = depVariable,
				effOfInd = eff,
				errorReturn = NULL
			)
			if (!is.null(lstTmp)) {
				lst[[counter1]] = lstTmp
				names(lst)[counter1] = eff
				counter1  = counter1  + 1
			}
			# 2. All subset interactions effect sizes
			if (length(cats) > 1) {
				for (j in 1:length(cats)) {
					cmbn = combn(x = cats, m = j)
					for (k in 1:ncol(cmbn)) {
						interact = interaction(data[, cmbn[, k] , drop = FALSE], sep = '.', drop = TRUE)
						for (lvl in levels(interact)) {
							message0('\tLevel:', pasteComma(unlist(lvl)))
							olstTmp = eff.size(
								object = object,
								data = droplevels(data[interact %in% lvl, ]),
								depVariable = depVariable,
								errorReturn = NULL,
								effOfInd = eff
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
			extractAIC(x, k = 0)[2L]
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
	bAIC <- extractAIC(fit, scale, k = k, ...)
	edf <- bAIC[1L]
	bAIC <- bAIC[2L]
	if (is.na(bAIC))
		stop("AIC is not defined for this model, so 'stepAIC' cannot proceed")
	if (bAIC == -Inf)
		stop("AIC is -infinity for this model, so 'stepAIC' cannot proceed")
	nm <- 1
	Terms <- terms(fit)
	if (trace) {
		cat("Start:  AIC=",
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
			ffac <- ffac[-st,]
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
			aod <- aod[nzdf,]
			if (is.null(aod) || ncol(aod) == 0)
				break
			nc <- match(c("Cp", "AIC"), names(aod))
			nc <- nc[!is.na(nc)][1L]
			o <- order(aod[, nc])
			if (trace) {
				print(aod[o,])
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
		bAIC <- extractAIC(fit, scale, k = k, ...)
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
		ans[1,] <- extractAIC(object, scale, k = k, ...)
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
			ans[i + 1,] <- extractAIC(nfit, scale, k = k, ...)
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
		aod <- aod[o,]
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
