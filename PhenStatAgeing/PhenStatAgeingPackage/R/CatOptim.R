# Categorical core
crunner = function(object              ,
									 formula       = ~ category + Genotype + Sex + LifeStage,
									 rep           = 2000,
									 method        = NULL,
									 categorical   = NULL)
{
	if (method != 'FE' || categorical != TRUE) {
		message0 ('Improper method for the type of data')
		return(NULL)
	}
	####
	depVariable = all.vars(formula)[1]
	vars        = all.vars(formula)[-1]
	l           = lcomb = names = alTbls = NULL
	CmbiVars    = (length(vars) > 1)
	####
	orgobject   = object
	object      = object@datasetPL
	Obj         =  CheckMissing(object, reformulate(termlabels = depVariable, response = NULL))
	object      = Obj$new.data
	sta.time    = Sys.time()
	message0('Fisher exact test with ',
					 ifelse(rep > 0, rep, 'No'),
					 ' repetition in progress ...')
	####
	counter  = 1
	for (j in 1:length(vars)) {
		l[[counter]] = ctest(x = object,
												 formula = reformulate(
												 	termlabels = c(depVariable, vars[j]),
												 	response = NULL
												 ))
		names = c(names, vars[j])
		counter = counter + 1
	}
	names(l)       = names
	if (CmbiVars) {
		alTbls = AllTables(
			dframe = as.data.frame(xtabs(formula = formula, data = object))  ,
			vars   = vars                                                    ,
			cl     = 3:(length(vars) + 2)                                    ,
			cols   = NULL                                                    ,
			response.name = depVariable                                      ,
			shrink = TRUE
		)
		lcomb = c(lcomb, lapply(alTbls, function(x) {
			ctest(x = x,
						formula = Freq ~ .,
						rep     = rep)
		}))
		names(lcomb) = paste(lapply(
			names(lcomb),
			FUN = function(nam) {
				n1 = names(dimnames(lcomb[[nam]]$table))
				n2 = paste0(n1[!n1 %in% depVariable], collapse = '.')
				return(n2)
			}
		), names(lcomb),
		sep = '_')
		l = c(l, lcomb)
	}
	message0('Function executed in ', (round(Sys.time() - sta.time, 2)), ' seconds')
	return(list(
		output = l[!duplicated(l)]         ,
		input  = list(
			PhenListAgeing  = orgobject      ,
			depVariable     = depVariable    ,
			rep             = rep            ,
			method          = method         ,
			categorical     = categorical
		),
		extra  = list(
			missings        = Obj$missings   ,
			UsedData        = Obj$new.data   ,
			OrgData         = Obj$org.data   ,
			table           = alTbls
		)
	))
}
# Count summary core
summary.PhenlistCategoricalAgeingModel = function(object, ...) {
	if (is.null(object))
		message0('NULL object\n')
	
	message0('Fisher exact test results based on ',
			object$input$rep,
			' iterations \n')
	
	l = l2 = list(
		GenPval    = object$output$Genotype$result$p.value,
		SexPval    = object$output$Sex$result$p.value,
		LifeStagePval      = object$output$LifeStage$result$p.value,
		GenForLatePval     = object$output,
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
		stop('\n ~> NULL object\n')
	
	cat('\n ~> There is no plot available for the categorical data \n')
	
}

