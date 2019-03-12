# Categorical OPT core
crunner = function(object              ,
									 formula          = ~ category + Genotype + Sex + LifeStage,
									 rep              = 1500,
									 method           = NULL,
									 fullComparisions = TRUE)
{
	if (sum(method != c('FE', 'RR')) == 0 ||
			is.null(object)                   ||
			is.null(all.vars0(formula))) {
		message0 ('Improper method (',
							method,
							') for the type of data, or the "formula" is left blank')
		return(NULL)
	}
	message0('FE framework in progress ...')
	sta.time    = Sys.time()
	message0('Top framework: ', method)
	message0('Fisher exact test with ',
					 ifelse(rep > 0, rep, 'No'),
					 ' repetition(s) in progress ...')
	newFormula    = checkModelTermsInData(
		formula = formula,
		data = object@datasetPL,
		responseIsTheFirst = TRUE
	)
	####
	allTerms    = all.vars(newFormula)
	if (length(allTerms) < 1) {
		message0('No variable to test.')
		return(NULL)
	}
	####
	lComplete   = NULL
	for (indx in 1:ifelse(fullComparisions, pmax(1, length(allTerms) - 1), 1)) {
		depVariable = allTerms[indx]
		vars        = allTerms[-c(1:indx)]
		l           = lcomb = names = alTbls = NULL
		CmbiVars    = (length(vars) > 1)
		####
		message0('Step ', indx, '. Testing "', depVariable,'"')
		####
		newObject   = object@datasetPL
		Obj         = CheckMissing(newObject,
															 reformulate(termlabels = depVariable, response = NULL))
		newObject   = Obj$new.data
		####
		counter  = 1
		for (j in 1:length(vars)) {
			message0('Testing for the main effect: ',
							 pasteComma(vars[j], replaceNull = FALSE))
			l[[counter]] = ctest(x = newObject,
													 formula = reformulate(
													 	termlabels = c(depVariable, vars[j]),
													 	response = NULL
													 ))
			names = c(names, vars[j])
			counter = counter + 1
		}
		names(l)       = names
		if (CmbiVars) {
			message0('Combined effects in progress ...')
			alTbls = AllTables(
				dframe = as.data.frame(xtabs(
					formula = newFormula, data = newObject
				))  ,
				vars   = vars                                                    ,
				cl     = 3:(length(vars) + 2)                                    ,
				cols   = NULL                                                    ,
				response.name = depVariable                                      ,
				shrink = TRUE
			)
			message0('Testing for the combined effects ... ')
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
		# fix naming issue (make them short)
		names(l) = gsub("\\.{4}.*", "", names(l))
		lComplete[[depVariable]] = l
	}
	
	message0('Total tested categories = ',
					 length(lComplete),
					 ': ',
					 pasteComma(names(lComplete)))
	message0('\tTotal tests =  ', sum(sapply(lComplete, length)))
	message0('FE framework executed in ', round(difftime(Sys.time() , sta.time, units = 'sec'), 2), ' seconds')
	OutR = list(
		output = list(SplitModels = lComplete) ,
		input  = list(
			PhenListAgeing  = object   ,
			data            = object@datasetPL   ,
			depVariable     = allTerms[1]        ,
			rep             = rep                ,
			method          = method             ,
			formula         = formula
		),
		extra  = list(
			missings        = Obj$missings   ,
			UsedData        = Obj$new.data   ,
			AllTable        = alTbls         ,
			Cleanedformula  = newFormula
		)
	)
	class(OutR) <- 'PhenStatAgeingFE'
	return(OutR)
}




# 
# 
# 
# 
# # Count summary core
# summary.PhenlistCategoricalAgeingModel = function(object, ...) {
# 	if (is.null(object))
# 		message0('NULL object\n')
# 	
# 	message0('Fisher exact test results based on ',
# 					 object$input$rep,
# 					 ' iterations \n')
# 	
# 	l = l2 = list(
# 		GenPval    = object$output$Genotype$result$p.value,
# 		SexPval    = object$output$Sex$result$p.value,
# 		LifeStagePval      = object$output$LifeStage$result$p.value,
# 		GenForLatePval     = object$output,
# 		GenForEarlyPval    = object$test.early$p.value,
# 		GenEsize           = object$Genotype.eff,
# 		SexEsize           = object$Sex.eff,
# 		LifestageEsize     = object$LifeStage.eff,
# 		LateStageEsize     = object$LateGenotype.eff,
# 		EarlyStageEsize    = object$EarlyGenotype.eff,
# 		otherPvals         = paste(
# 			names(object$AllOtherTests),
# 			lapply(object$AllOtherTests, function(x) {
# 				paste(x$result$p.value, collapse = '')
# 			}),
# 			sep = ': ',
# 			collapse = ';  '
# 		)
# 	)
# 	l = lapply(l, function(x) {
# 		if (is.null(x)) {
# 			x = '-'
# 		} else{
# 			x = x
# 		}
# 	})
# 	cat(
# 		'\nGenotype pval                       :',
# 		l$GenPval,
# 		'\nSex pval                            :',
# 		l$SexPval,
# 		'\nLifeStage pval                      :',
# 		l$LifeStagePval,
# 		'\nGenotype for only late stage  pval  :',
# 		l$GenForLatePval,
# 		'\nGenotype for only early stage pval  :',
# 		l$GenForEarlyPval,
# 		'\nGenotype effect size                :',
# 		l$GenEsize,
# 		'\nSex effect size                     :',
# 		l$SexEsize,
# 		'\nLifestage effect size               :',
# 		l$LifestageEsize,
# 		'\nEarly Stage Genotype effect size    :',
# 		l$EarlyStageEsize,
# 		'\nLate Stage Genotype effect size     :',
# 		l$LateStageEsize,
# 		'\nOther p-values                      :',
# 		l$otherPvals,
# 		'\n'
# 	)
# 	
# 	outp = list(object = object     ,
# 							summary = l2)
# 	outp$JSON = toJSONI(outp)
# 	return(invisible(outp))
# }
# 
# 
