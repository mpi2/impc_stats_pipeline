# Reference range OPT core
RRrunner = function(object              ,
										formula       = ~ category + Genotype + Sex + LifeStage,
										rep           = 1500,
										method        = NULL,
										RRprop        = .05,
										...)
{
	sta.time    = Sys.time()
	allTerms    = all_vars0(x = formula)
	if (!method %in% c('RR')        ||
			is.null(allTerms)           ||
			is.null(object)             ||
			sum(allTerms %in% names(object@datasetPL)) < 2) {
		#####
		message0 (
			'Improper method (',
			method,
			') for the type of data, or the `formula` is not properly specified/left blank. \n\tFormula: ',
			printformula(formula)
		)
		return(NULL)
	}
	message0('RR+ framework in progress ...')
	if (is.null(RRprop)     ||
			!is.numeric(RRprop) ||
			RRprop <= 0.5       ||
			RRprop >= 1.0) {
		message0('"RRprop" must be a value greater than 0.5 and less than 1')
		warnings('Improper value for "RRprop"')
		return(NULL)
	}
	cleanFormulaForOutput   = checkModelTermsInData(
		formula = formula,
		data = object@datasetPL,
		responseIsTheFirst = TRUE
	)
	message0('Discritizing the continuous data into discrete data. The quantile = ',
					 RRprop)
	message0('Stp 1. Low versus Normal/High')
	RRobject_low = RRNewObjectAndFormula(
		object = object                   ,
		RRprop = 1 - RRprop                 ,
		formula = formula                 ,
		labels = c('Low', 'NormalHigh')   ,
		depVarPrefix = 'Low'
	)
	RRresult_low = crunner(
		object = RRobject_low$newobject   ,
		formula = RRobject_low$newFormula ,
		rep = rep                         ,
		method = 'RR'                     ,
		fullComparisions = TRUE           ,
		noteToFinish = 'in step 1'        ,
		...
	)
	
	message0('Stp 2. Low/Normal versus High')
	RRobject_high = RRNewObjectAndFormula(
		object = object                    ,
		RRprop = RRprop                    ,
		formula = formula                  ,
		labels = c('LowNormal', 'High')    ,
		depVarPrefix = 'High'
	)
	RRresult_high = crunner(
		object = RRobject_high$newobject   ,
		formula = RRobject_high$newFormula ,
		rep = rep                          ,
		method = 'RR'                      ,
		fullComparisions = TRUE            ,
		noteToFinish = 'in step 2'         ,
		...
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
	OutR = list(
		output = list(SplitModels = SpltResult),
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
			UsedData        = list(Low  = RRresult_low$input$data      ,
														 High = RRresult_high$input$data)    ,
			AllTable        = list(
				Low  = RRresult_low$extra$AllTable                       ,
				High = RRresult_high$extra$AllTable
			)                                                          ,
			Cleanedformula           = cleanFormulaForOutput           ,
			DiscritiseCleanedformula = list(
				Low  = RRresult_low$extra$Cleanedformula                 ,
				High = RRresult_high$extra$Cleanedformula
			)
		)
	)
	class(OutR) <- 'PhenStatAgeingRR'
	return(OutR)
}
