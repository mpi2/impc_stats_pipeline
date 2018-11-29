# Categorical core
crunner = function(object,
									 Sex           = TRUE,
									 LifeStage     = TRUE,
									 depVariable   = NULL,
									 rep           = 10 ^ 5,
									 method        = NULL,
									 categorical   = NULL) 
{
	if (method != 'FE' || categorical != TRUE)
		stop ('\n => Improper method for the type of data \n')
	orgobject = object
	object    = object@datasetPL
	Obj       =  CheckMissing(object, reformulate(termlabels = depVariable, response = NULL))
	object    = Obj$new.data
	sta.time = Sys.time()
	cat('\n => Fisher exact test with ', rep, ' repetition in progress ... \n')
	# Genotype
	G  = reformulate(termlabels = c(depVariable, 'Genotype'),
									 response = NULL)
	GT = ctest(x       = object,
						 formula = G,
						 rep     = rep)
	GTE  = GT$effect
	GT   = GT$result
	############
	# For all test 
	f = as.formula(paste(
		'~',
		depVariable,
		'+Genotype',
		ifelse(Sex, '+Sex', ''),
		ifelse(LifeStage, '+LifeStage', '')
	))
	alTst   = xtabs(formula = f, data = object)
	allGtst = AllCombinations = all.tables(
		as.data.frame(alTst),
		cl   = 3,
		cols = c(depVariable,'Genotype'),
		response.name = depVariable,
		shrink = TRUE
	)
	allGtst = lapply(allGtst, function(x) {
		ctest(x = x,
					formula = Freq ~ .,
					rep     = rep)
	})
	#add one ome for the life stage
	############
	if (Sex) {
		# Sex
		S  = reformulate(termlabels = c(depVariable, 'Sex'),
										 response = NULL)
		ST = ctest(x = object,
							 formula = S,
							 rep = rep)
		STE = ST$effect
		ST  = ST$result
		##### all for sex
		allStst = all.tables(
			as.data.frame(alTst),
			cl = 3,
			cols = c(depVariable,'Sex'),
			response.name = depVariable,
			shrink = TRUE
		)
		AllCombinations = append(AllCombinations,allStst)
		allStst = lapply(allStst, function(x) {
			ctest(x = x,
						formula = Freq ~ .,
						rep = rep)
		})
		allGtst = append(allGtst, allStst)
		###### end
	} else{
		ST = STE = NULL
	}
	if (LifeStage) {
		# LifeStage
		A  = reformulate(termlabels = c(depVariable, 'LifeStage'),
										 response = NULL)
		AT = ctest(x = object,
							 formula = A,
							 rep = rep)
		# Early/Late
		AEL  = reformulate(
			termlabels = c(depVariable, 'Genotype', 'LifeStage'),
			response = NULL
		)
		ATE = ctest(
			x = object,
			formula = AEL,
			asset = 1,
			rep = rep
		)
		ATL = ctest(
			x = object,
			formula = AEL,
			asset = 2,
			rep = rep
		)
		
		# 1
		ATGE = AT$effect
		AT   = AT$result
		# 2
		ATEE = ATE$effect
		ATE  = ATE$result
		# 3
		ATLE = ATL$effect
		ATL  = ATL$result
		
		##### all for LifeStage
		allLStst = all.tables(
			as.data.frame(alTst),
			cl = 3,
			cols = c(depVariable,'LifeStage'),
			response.name = depVariable,
			shrink = TRUE
		)
		AllCombinations = append(AllCombinations,allLStst)
		allLStst = lapply(allLStst, function(x) {
			ctest(x = x,
						formula = Freq ~ .,
						rep = rep)
		})
		allGtst = append(allGtst, allLStst)
		###### end
	} else{
		AT = ATE = ATL = ATGE = ATEE = ATLE = NULL
	}
	
	
	cat('\n => Function executed in (',
			round(Sys.time() - sta.time, 4),
			'seconds )\n')
	return(
		list(
			Genotype.eff = GTE       ,
			LifeStage.eff = ATGE     ,
			Sex.eff       = STE      ,
			LateGenotype.eff = ATLE  ,
			EarlyGenotype.eff = ATEE ,
			test.Genotype = GT       ,
			test.Sex = ST            ,
			test.LifeStage = AT      ,
			test.late = ATL          ,
			test.early = ATE         ,
			AllOtherTests = allGtst  ,
			data  = orgobject        ,
			Sex = Sex                ,
			LifeStage = LifeStage    ,
			depVariable = depVariable,
			rep = rep                ,
			AllCombinations=AllCombinations  ,
			missings        = Obj$missings   ,
			OriginalDataSet = Obj$org.data   ,
			method = method                  ,
			categorical = categorical
		)
	)
}


# Count summary core
summary.PhenlistCategoricalAgeingModel = function(object, ...) {
	if (is.null(object))
		stop('\n => NULL object\n')
	
	cat('\n => Fisher exact test results based on ',
			object$rep,
			' iterations \n')
	
	l = l2 = list(
		GenPval    = object$test.Genotype$p.value,
		SexPval    = object$test.Sex$p.value,
		LifeStagePval      = object$test.LifeStage$p.value,
		GenForLatePval     = object$test.late$p.value,
		GenForEarlyPval    = object$test.early$p.value,
		GenEsize           = object$Genotype.eff,
		SexEsize           = object$Sex.eff,
		LifestageEsize     = object$LifeStage.eff,
		LateStageEsize     = object$LateGenotype.eff,
		EarlyStageEsize    = object$EarlyGenotype.eff,
		otherPvals         = paste(
			names(object$AllOtherTests),
			lapply(object$AllOtherTests, function(x) {
				paste(x$result$p.value, collapse = '')
			}),
			sep = ': ',
			collapse = ';  '
		)
	)
	l = lapply(l, function(x) {
		if (is.null(x)) {
			x = '-'
		} else{
			x = x
		}
	})
	cat(
		'\nGenotype pval                       :',
		l$GenPval,
		'\nSex pval                            :',
		l$SexPval,
		'\nLifeStage pval                      :',
		l$LifeStagePval,
		'\nGenotype for only late stage  pval  :',
		l$GenForLatePval,
		'\nGenotype for only early stage pval  :',
		l$GenForEarlyPval,
		'\nGenotype effect size                :',
		l$GenEsize,
		'\nSex effect size                     :',
		l$SexEsize,
		'\nLifestage effect size               :',
		l$LifestageEsize,
		'\nEarly Stage Genotype effect size    :',
		l$EarlyStageEsize,
		'\nLate Stage Genotype effect size     :',
		l$LateStageEsize,
		'\nOther p-values                      :',
		l$otherPvals,
		'\n'
	)
	
	outp = list(object = object     ,
							summary = l2)
	outp$JSON = toJSONI(outp)
	return(invisible(outp))
}


# Plot
plot.PhenlistCategoricalAgeingModel = function(x, ...) {
	if (is.null(x))
		stop('\n => NULL object\n')
	
	cat('\n => There is no plot available for the categorical data \n')
	
}
# Test engine
ctest = function(x,
								 formula,
								 asset = NULL,
								 rep = 10 ^ 5,
								 ...) {
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
			stop ('\n => There are some empty levels in data\n')
		}
	}
	if (any(dim(xtb) == 1)) {
		r = list(p.value = 1)
		effect = 1
	} else{
		r = fisher.test(
			xtb,
			simulate.p.value = TRUE,
			conf.int = FALSE,
			B = rep,
			...
		)
		#https://stats.stackexchange.com/questions/22508/effect-size-for-fishers-exact-test
		effect     = abs(xtb[2, 2] / (xtb[1, 2] + xtb[2, 2]) - xtb[2, 1] / (xtb[1, 1] +
																																					xtb[2, 1]))
		r$formula   = formula
		r$table     = xtb
	}
	return(list(result = r, effect = effect))
}

# Change with super extra care
# This is a complicated function for modeling all variation of the variables
all.tables = function(dframe        = NULL,  # dataframe
											cl            = 0   ,  # lock the number of columns: works only if shrinke = TRUE
											response.name = NULL,  # response name in the beginning of each variable name [only]
											cols          = NULL,  # filter on the columns of interest : works only if shrinke = TRUE
											shrink        =  FALSE # remove no variation levels (columns)
) {
	
	if(is.null(dframe))
		stop('\n => Null data frame \n')
	
	l2 = list()
	lnames = c()
	counter = 0
	
	allCategorical = sapply(dframe, function(x) {
		is.factor(x)
	})
	
	cat = which(allCategorical == TRUE)
	con = which(allCategorical == FALSE)
	
	for (i in seq(cat)) {
		cb = combn(cat, i)
		for (j in 1:ncol(cb)) {
			out = split(dframe, apply(dframe, 1, function(x)
				paste(x[-c(cb[, j], con)], collapse = ".")))
			l2[(counter + 1):(counter + length(out))] = out[1:length(out)]
			lnames = c(lnames, names(out[1:length(out)]))
			counter = counter + length(out)
		}
	}
	
	lnames[lnames == ''] = 'All'
	names(l2) = lnames
	
	if (shrink)
		l2 = lapply(l2, function(x) {
			# remove fixed value columns
			x[, !c(apply(x, 2, function(xx) {
				all(na.omit(xx) == na.omit(xx)[1])
			}) & allCategorical)]
		})
	
	if (!is.null(cols))
		l2 = l2[as.logical(lapply(
			# keep certain columns
			l2,
			FUN = function(x) {
				all(cols %in% names(x))
			}
		))]
	
	if (all(cl > 0))
		# Fixed number of columns in each element of the list
		l2 = l2[which(lapply(
			l2,
			FUN = function(x) {
				ncol(x) %in% cl
			}
		) == TRUE)]
	
	if(!is.null(response.name))
		names(l2) = paste(response.name, names(l2), sep = ': ')
	
	return(l2)
}
