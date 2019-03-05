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
pasteComma = function(..., replaceNull = TRUE) {
	if (replaceNull)
		paste(replaceNull(list(...), replaceBy = 'NULL'),
					sep = ', ',
					collapse = ', ')
	else
		paste(..., sep = ', ', collapse = ', ')
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
			'Some parts of the model does not included in the data. See: ',
			pasteComma(vars[!In], replaceNull = FALSE)
		)
		ft  = vars [!In]
		formula = update.formula(formula,
														 reformulate0(
														 	termlabels = c('.', ft, paste0(ft, ':.')),
														 	response   = NULL,
														 	intercept  = TRUE,
														 	sep        = '-'
														 ))
	}
	return(formula)
}

# Test type (decimal [Continues]/not decimal [Categorical])
IsCategorical = function(num = NA, max.consider = Inf) {
	if (length(num) < 1 ||
			is.null(num) || length(na.omit(num)) < 1)
		# do not change the order
		stop('\n ~> Null vector or all NAs \n')
	
	num = na.omit(num)
	if (!is.numeric(num)) {
		cat = TRUE
		method = 'FE'
	} else if (all(num %% 1 == 0) &&
						 length(unique(num)) < max.consider) {
		cat = TRUE
		method = 'GLM'
	} else{
		cat = FALSE
		method = 'MM'
	}
	return(list(cat = cat, method = method))
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

suppressMessagesANDErrors = function(exp, debug = TRUE) {
	if (!debug) {
		suppressMessages(exp)
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
	if (lvars > 0) {
		for (i in 1:lvars) {
			message0(i,
							 ' of ',
							 lvars,
							 '. Checking for the feasibility of the terms and interactions ...')
			cmb = combn(vars, i)
			for (j in 1:ncol(cmb)) {
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

variablesInData = function(df, names) {
	if (is.null(df) || is.null(names) || sum(names %in%  names(df)) < 1)
		return(NULL)
	newNames = names[names %in%  names(df)]
	message0('Variables that are found in data: ',pasteComma(newNames))
	return(newNames)
}

ComplementaryFeasibleTermsInContFormula = function(formula, data) {
	message0(
		'Checking for the feasibility of the terms and interactions. Formula: ',
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
				'The following term(s) removed due to No data  or no data in the interactions: ',
				pasteComma(fbm[fbm$min.freq <= 0, c('names')], replaceNull = FALSE)
			)
	}
	return(formula)
}
# check for missings
CheckMissing = function(data, formula) {
	org.data = data
	new.data = data[complete.cases(data[, all.vars(formula)]),]
	missings = ifelse(all(dim(org.data) == dim(new.data)),	0, dim(org.data)[1] -
											dim(new.data)[1])
	if (missings)
		message0('The data contain ',
						 missings,
						 ' missing(s).\nMissing data removed. \n ')
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
										effOfInd = 'Genotype',
										errorReturn = NULL,
										debug = FALSE) {
	if (all(is.null(data)))
		data = getData(object)
	# if (length(levels(data[, effOfInd])) != 2) {
	# 	message0('There must be two levels in ',
	# 					paste(effOfInd, collapse = ' and '))
	# 	return(NULL)
	# }
	f    = reformulate(termlabels = effOfInd, depVariable)
	agr  = aggregate(f, data = data, FUN = mean)
	if (debug) {
		print(agr)
	}
	if (any(dim(agr) < 2)) {
		message0(
			' Effect size estimation: No variation in "',
			paste(effOfInd, collapse = ' or '),
			'"; or less than two levels in "',
			paste0(effOfInd, collapse = ','),
			'"'
		)
		return(errorReturn)
	}
	MDiff        = max(dist(agr[, depVariable]))
	NModel       = update(object,  as.formula(paste('~', effOfInd)), data = data)
	CoefEffSizes = sapply(data[, effOfInd], FUN = is.numeric)
	if (!is.null(NModel)) {
		if (sum(CoefEffSizes)) {
			# For continues variables it is the coefficient
			efSi         = list(value = as.list(coef(NModel))[[effOfInd]][1], type =
														'coefficient')
		} else{
			# For categorical variables it is the mean difference
			r            = resid(NModel)
			efSi         = list(value = ifelse(sd(r) > 0, abs(MDiff) / sd(r), NA), type =
														'Mean difference')
		}
	} else{
		efSi = NULL
	}
	return(efSi)
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
								'Genotype',
								arg,
								collapse = ':',
								sep = ':'
							),
							argsCon
						),
						intercept = TRUE
					)
					
					l0    = tryCatch(
						update(F.Model,
									 newModel),
						error = function(e) {
							message0(e)
							return(NULL)
						} ,
						warning = function(w) {
							message0(w)
							return(NULL)
						}
					)
					message0(
						'Tested model: ',
						printformula(newModel),
						ifelse(!is.null(l0), ' [Successful]', ' [Failed]'),
						breakLine = FALSE
					)
					if (!is.null(l0)) {
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
	paste01(format(formula, trim = TRUE, width = 0), collapse = '')
}
# Categorical effect size
# https://github.com/mpi2/stats_working_group/raw/master/PhenStatUserGuide/PhenStatUsersGuide.pdf p122
cat.eff.size = function(xtb) {
	if (any(dim(xtb) < 1)) {
		r = NULL
	} else{
		r = max(apply(prop.table(xtb, margin = 2), 1, function(x) {
			max(dist(x, method = 'maximum', diag = TRUE), na.rm = TRUE)
		}), na.rm = TRUE)
	}
	return(r)
}
# Test engine
ctest = function(x,
								 formula,
								 asset = NULL,
								 rep = 1500,
								 ...) {
	message = c()
	#message0('\tForming the cross table ...')
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
		#message0('\tFisher exact test procedure in progress ...')
		r = tryCatch(
			fisher.test(
				xtb,
				simulate.p.value = rep > 0,
				conf.int = FALSE,
				B = rep,
				...
			),
			error = function(e) {
				message0(e)
				return(NULL)
			} ,
			warning = function(w) {
				message0(w)
				return(NULL)
			}
		)
		effect     = cat.eff.size(xtb)
		r$formula   = formula
		r$table     = xtb
	}
	#message0('\tGathering the results for the fisher exact test procedure ...')
	return(
		list(
			result      = r,
			effectSize  = list(value = effect, type = 'Proportion change'),
			note        = message,
			table       = xtb,
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
	message0('\tSplitting in progress...')
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
	
	message0('\tKeep variables of interest ...')
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
	message0('\tdichotomising the final tables ...')
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
	message0('\tfinalising the tables ....')
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

formulaTerms = function(formula) {
	r = attr(terms(as.formula(formula))    , which = 'term.labels')
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


ModelInReference = function(model, reference, responseIncluded = FALSE,veryLower = ~Genotype+1) {
	mo = formulaTerms(model)
	re = formulaTerms(reference)
	r  = re[re %in% mo]
	if (length(r) > ifelse(responseIncluded, 1, 0)) {
		message0('Some terms in the "lower" parameter ignored. See:\n\t',
						 pasteComma(mo[!(mo %in% re)]))
		if (responseIncluded)
			out = reformulate(termlabels = r[-1],
												response   = r[1],
												intercept  = TRUE)
		else
			out = reformulate(termlabels = r,
												response   = NULL,
												intercept  = TRUE)
	} else{
		out = veryLower
	}
	message0('The corrected lower: ', printformula(out))
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
	# SexLab = ifelse(
	# 	is.na(object$input$PhenListAgeing@dataset.colname.sex),
	# 	'Sex',
	# 	object$input$PhenListAgeing@dataset.colname.sex
	# )
	
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
	# GenotypeLab = ifelse(
	# 	is.na(object$input$PhenListAgeing@dataset.colname.genotype),
	# 	'Genotype',
	# 	object$input$PhenListAgeing@dataset.colname.genotype
	# )
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
	# BatchLab = ifelse(
	# 	is.na(object$input$PhenListAgeing@dataset.colname.batch),
	# 	'Batch',
	# 	object$input$PhenListAgeing@dataset.colname.batch
	# )
	BatchLab = 'Batch'
	##########
	# WeightLab = ifelse(
	# 	is.na(object$input$PhenListAgeing@dataset.colname.weight),
	# 	'Weight',
	# 	object$input$PhenListAgeing@dataset.colname.weight
	# )
	WeightLab = 'Weight'
	##########
	l = list(
		LifeStage = list(
			LifeStage    = 'LifeStage'   ,
			Early        = 'Early',
			Late         = 'Late'  ,
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
	controls = subset(object@datasetPL, Genotype %in% object@refGenotype)
	mutants  = subset(object@datasetPL, Genotype %in% object@testGenotype)
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
		qntl[duplicated(qntl)] = jitter(x = qntl[duplicated(qntl)], amount = 10 ^
																			-JitterPrecision)
		qntl = sort(qntl)
	}
	message0('quantile for cutting the data = ',
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
	controls = subset(object@datasetPL, Genotype %in% object@refGenotype)
	mutants  = subset(object@datasetPL, Genotype %in% object@testGenotype)
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
		#message0(dataPointsSummary[3])
		variabilityThreshold <- NoCombinations
		#if (NoCombinations==4)
		#variabilityThreshold <- 3
		
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
	message0('Extended formula (if special characters included): ',
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
#
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
