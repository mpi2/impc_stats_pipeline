# Reference range OPT core
RRrunner = function(object              ,
										formula       = ~ category + Genotype + Sex + LifeStage,
										rep           = 1500,
										method        = NULL,
										RRprop        = .05)
{
	sta.time    = Sys.time()
	allTerms    = all.vars0(x = formula)
	if (!method %in% c('RR')        ||
			is.null(allTerms)           ||
			is.null(object)) {
		#####
		message0 ('Improper method (',
							method,
							') for the type of data, or the "formula" is left empty.')
		return(NULL)
	}
	
	message0('Discritizing the continuous data into discrite data. The quantile = ',
					 RRprop)
	message0('Stp 1. Low versus Normal/High')
	RRobject_low = RRNewObjectAndFormula(
		object = object,
		RRprop = RRprop,
		formula = formula ,
		labels = c('Low', 'NormalHigh')
	)
	RRresult_low = crunner(
		object = RRobject_low$newobject,
		formula = RRobject_low$newFormula,
		rep = rep,
		method = 'RR'
	)
	
	message0('Stp 2. Low/Normal versus High')
	RRobject_high = RRNewObjectAndFormula(
		object = object,
		RRprop = 1 - RRprop,
		formula = formula ,
		labels = c('LowNormal', 'High')
	)
	RRresult_high = crunner(
		object = RRobject_high$newobject,
		formula = RRobject_high$newFormula,
		rep = rep,
		method = 'RR'
	)
	message0('RR framework executed in ', round(difftime(Sys.time() , sta.time, units = 'sec')), ' seconds')
	OutR = list(
		output = list(
			SplitModels = list(
				Low  = RRresult_low$output$SplitModels  ,
				High = RRresult_high$output$SplitModels
			)
		) ,
		input  = list(
			PhenListAgeing  = object           ,
			data            = object@datasetPL ,
			depVariable     = allTerms[1]      ,
			rep             = rep              ,
			method          = method           ,
			formula         = formula          ,
			RRprop          = RRprop
		),
		extra  = list(
			missings        = RRresult_low$extra$missings                                                               ,
			UsedData        = RRresult_low$input$data[!names(RRresult_low$input$data) %in% RRobject_low$newDepVariable] ,
			AllTable        = list(
				Low  = RRresult_low$extra$AllTable                                                                        ,
				High = RRresult_high$extra$AllTable
			)                                                                                                           ,
			Cleanedformula  = RRresult_low$extra$Cleanedformula
		)
	)
	class(unclass(OutR)) = 'PhenStatAgeingRR'
	return(OutR)
}
