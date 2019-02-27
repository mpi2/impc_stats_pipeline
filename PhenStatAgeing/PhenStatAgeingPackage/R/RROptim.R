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
							') for the type of data, or the "formula" is left blank')
		return(NULL)
	}
	message0('RR+ framework in progress ...')
	cleanFormulaForOutput   = checkModelTermsInData(
		formula = formula,
		data = object@datasetPL,
		responseIsTheFirst = TRUE
	)
	message0('Discritizing the continuous data into discrite data. The quantile = ',
					 RRprop)
	message0('Stp 1. Low versus Normal/High')
	RRobject_low = RRNewObjectAndFormula(
		object = object                   ,
		RRprop = RRprop                   ,
		formula = formula                 ,
		labels = c('Low', 'NormalHigh')   ,
		depVarPrefix = 'Low'
	)
	RRresult_low = crunner(
		object = RRobject_low$newobject   ,
		formula = RRobject_low$newFormula ,
		rep = rep                         ,
		method = 'RR'                     ,
		fullComparisions = TRUE
	)
	
	message0('Stp 2. Low/Normal versus High')
	RRobject_high = RRNewObjectAndFormula(
		object = object                    ,
		RRprop = 1 - RRprop                ,
		formula = formula                  ,
		labels = c('LowNormal', 'High')    ,
		depVarPrefix = 'High'
	)
	RRresult_high = crunner(
		object = RRobject_high$newobject   ,
		formula = RRobject_high$newFormula ,
		rep = rep                          ,
		method = 'RR'                      ,
		fullComparisions = TRUE
	)
	message0('RR framework executed in ', round(difftime(Sys.time() , sta.time, units = 'sec'), 2), ' seconds')
	#####
	SpltResult = c(
		renameVariableNameInList(
			list = RRresult_low$output$SplitModels ,
			name = RRobject_low$newDepVariable,
			prefix  = 'Low',
			not     = TRUE
		)  ,
		renameVariableNameInList(
			list = RRresult_high$output$SplitModels,
			name = RRobject_high$newDepVariable,
			prefix  = 'High',
			not = TRUE
		)
	)
	#names(RRresult_low$output$SplitModels)  = paste('Low' , names(RRresult_low$output$SplitModels ), sep =  '_')
	#names(RRresult_high$output$SplitModels) = paste('High', names(RRresult_high$output$SplitModels), sep =  '_')
	OutR = list(
		output = list(SplitModels = SpltResult)                                  ,
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
			missings        = list(
				Low  = RRresult_low$extra$missings  ,
				High = RRresult_high$extra$missings
			),
			UsedData        = list(Low  = RRresult_low$input$data,
														 High = RRresult_high$input$data),
			#				RRresult_low$input$data[!names(RRresult_low$input$data) %in% RRobject_low$newDepVariable] ,
			AllTable        = list(
				Low  = RRresult_low$extra$AllTable                                                                        ,
				High = RRresult_high$extra$AllTable
			)                                                                                                           ,
			Cleanedformula           = cleanFormulaForOutput           ,
			DiscritiseCleanedformula = list(
				Low  = RRresult_low$extra$Cleanedformula ,
				High = RRresult_high$extra$Cleanedformula
			)
		)
	)
	class(OutR) <- 'PhenStatAgeingRR'
	return(OutR)
}
