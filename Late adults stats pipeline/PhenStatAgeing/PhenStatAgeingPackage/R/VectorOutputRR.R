vectorOutputRR =	function(object)
{
	if (!is.null(object$messages))
		return (NULL)
	#####################################################################
	Labels         = PhenListAgeingLevels(object = object)
	Fmodel         = object$extra$Cleanedformula
	frm            = formula(Fmodel)
	depVariable    = all_vars0(frm)[1]
	#	equation       = NULL
	formula        = printformula(frm)
	framework      = paste0('Reference Range Plus Test framework; quantile = ',
													object$input$RRprop)
	#	fittingMethod  = NULL
	#####################################################################
	Vsplit         = object$output$SplitModels
	low            = pasteUnderscore('Low' , depVariable, 'discretised')
	high           = pasteUnderscore('High', depVariable, 'discretised')
	VsplitLow      = Vsplit[[low]]
	VsplitHig      = Vsplit[[high]]
	#####################################################################
	x                = object$input$data
	columnOfInterest = x[, c(depVariable)]
	#####################################################################
	variability      = length(unique(columnOfInterest)) / length(columnOfInterest)
	#####################################################################
	DSsize            = SummaryStats(
		x = x,
		formula = object$input$formula,
		#label = 'Summary statistics',
		lower = TRUE,
		drop = TRUE,
		sep = ' '
	)
	MultiBatch = ifelse(multiBatch(x),
											'Dataset contains multi batches',
											'Dataset contains single batch')
	addInfo           = list(
		Data = list(
			'Data signature'          = dataSignature(formula = frm,
																								data   = x)  ,
			'Variability'            = variability                 ,
			'Summary statistics'      = DSsize
		),
		Analysis = list(
			# 'Formula'                = list(
			# 	input   = printformula(object$input$formula),
			# 	final   = printformula(formula)
			# ),
			'Model setting'          = extractFERRTerms(object),
			'Is model optimised'     = NULL                    , 
			'Multibatch in analysis' = MultiBatch,
			'Gender included in analysis' = ifelse(
				nlevels(x$Sex) > 1,
				'Both sexes included',
				paste0('Only one sex included in the analysis; ', levels(x$Sex))
			),
			'Further models' = if (!is.null(Vsplit)) {
				setNames(sapply(Vsplit, function(v) {
					lapply(
						v,
						FUN = function(v2) {
							list('p-val'       = v2$result$p.value,
									 'effect size' = v2$effectSize)
						}
					)
				}), nm = names(Vsplit))
			} else{
				NULL
			},
			'Effect sizes'                   = 'Look at the individual models',
			'Other residual normality tests' = NULL
		)
	)
	#####################################################################
	percentageChanges = NA
	#####################################################################
	vectorOutput      = list(
		'Method'                               = 	framework,
		'Dependent variable'                   =	pasteComma(depVariable,
																												pasteUnderscore(depVariable, 'discretised'),
																												truncate = FALSE), 
		'Batch included'                       =	 NULL,
		'Batch p-val'                          =  'Legacy',
		'Residual variances homogeneity'       =   NULL,
		'Residual variances homogeneity p-val' =  'Legacy',
		#####################################################################
		'Genotype contribution' =	list(
			Overal = lowHighList(
				VsplitLow$Genotype$result$p.value,
				VsplitHig$Genotype$result$p.value
			),
			'Sex FvKO p-val'   =	lowHighList(
				VsplitLow$Genotype_Female$result$p.value,
				VsplitHig$Genotype_Female$result$p.value
			),
			'Sex MvKO p-val'   =  lowHighList(
				VsplitLow$Genotype_Male$result$p.value,
				VsplitHig$Genotype_Male$result$p.value
			),
			'Sexual dimorphism detected' = 'Sex specific results for Low/High tables are always reported'
		),
		'Genotype estimate'           = list(
			low  = CatEstimateAndCI(VsplitLow$Genotype$result),
			high = CatEstimateAndCI(VsplitHig$Genotype$result)
		),
		'Genotype standard error'     = NULL, 
		'Genotype p-val'              = lowHighList(
			VsplitLow$Genotype$result$p.value  ,
			VsplitHig$Genotype$result$p.value  
		),
		'Genotype percentage change'  =	percentageChanges                                 ,
		'Genotype effect size'        = list(
			low  = VsplitLow$Genotype$effectSize     ,
			high = VsplitHig$Genotype$effectSize
		)   ,
		#####################################################################
		'Sex estimate'                =	list(
			low  = CatEstimateAndCI(VsplitLow$Sex$result),
			high = CatEstimateAndCI(VsplitHig$Sex$result)
		), 
		'Sex standard error'          = NULL,
		'Sex p-val'                   =	lowHighList(
			VsplitLow$Sex$result$p.value      ,
			VsplitHig$Sex$result$p.value      
		)   ,
		'Sex effect size'             =	list(
			low  = VsplitLow$Sex$effectSize   ,
			high = VsplitHig$Sex$effectSize
		)  ,
		#####################################################################
		'LifeStage estimate'          =	list(
			low  = CatEstimateAndCI(VsplitLow$LifeStage$result),
			high = CatEstimateAndCI(VsplitHig$LifeStage$result)
		), 
		'LifeStage standard error'    =	NULL,
		'LifeStage p-val'             =	lowHighList(
			VsplitLow$LifeStage$result$p.value ,
			VsplitHig$LifeStage$result$p.value 
		),
		'LifeStage effect size'       = list(
			low  = VsplitLow$LifeStage$effectSize     ,
			high = VsplitHig$LifeStage$effectSize
		)    ,
		#####################################################################
		'Weight estimate'             =	NULL,
		'Weight standard error'       =	NULL,
		'Weight p-val'                =	NULL,
		'Weight effect size'          = NULL,
		#####################################################################
		'Gp1 genotype'                     =	Labels$Genotype$Control		,
		'Gp1 Residuals normality test'     =	NULL                      ,
		'Gp2 genotype'                     =	Labels$Genotype$Mutant		,
		'Gp2 Residuals normality test'     =	NULL                      ,
		#####################################################################
		'Blups test'                       =  'Legacy',
		'Rotated residuals normality test' =  'Legacy',
		#####################################################################
		'Intercept estimate'               =	NULL,
		'Intercept standard error'         =	NULL,
		'Intercept p-val'                  =	NULL,
		#####################################################################
		'Interactions included'          =	list(
			'Genotype Sex'                   =  NULL,
			'Genotype LifeStage'             =  NULL,
			'Sex LifeStage'                  =  NULL,
			'Genotype Sex LifeStage'         =  NULL
		),
		#####################################################################
		################ interaction
		'Interactions p-val'      =	list(
			'Genotype Sex'            = NULL,
			'Genotype LifeStage'      = NULL,
			'Sex LifeStage'           = NULL,
			'Genotype Sex LifeStage'  = NULL
		),
		#####################################################################
		################ Sex interactions
		'Sex FvKO estimate'         = list(
			low  = CatEstimateAndCI(VsplitLow$Genotype_Female$result),
			high = CatEstimateAndCI(VsplitHig$Genotype_Female$result)
		), 
		'Sex FvKO standard error'   = NULL,
		'Sex FvKO p-val'            = lowHighList(
			VsplitLow$Genotype_Female$result$p.value ,
			VsplitHig$Genotype_Female$result$p.value 
		),
		'Sex FvKO effect size'      = list(
			low  = VsplitLow$Genotype_Female$effectSize     ,
			high = VsplitHig$Genotype_Female$effectSize
		)    ,
		#####################################################################
		'Sex MvKO estimate'         = list(
			low  = CatEstimateAndCI(VsplitLow$Genotype_Male$result),
			high = CatEstimateAndCI(VsplitHig$Genotype_Male$result)
		), 
		'Sex MvKO standard error'   = NULL,
		'Sex MvKO p-val'            = lowHighList(
			VsplitLow$Genotype_Male$result$p.value   ,
			VsplitHig$Genotype_Male$result$p.value  
		)  ,
		'Sex MvKO effect size'      = list(
			low  = VsplitLow$Genotype_Male$effectSize       ,
			high = VsplitHig$Genotype_Male$effectSize
		)      ,
		#####################################################################
		################ LifeStage interaction
		'LifeStage EvKO estimate'         =	list(
			low  = CatEstimateAndCI(VsplitLow$Genotype_Early$result),
			high = CatEstimateAndCI(VsplitHig$Genotype_Early$result)
		), 
		'LifeStage EvKO standard error'   =	NULL,
		'LifeStage EvKO p-val'            =	lowHighList(
			VsplitLow$Genotype_Early$result$p.value  ,
			VsplitHig$Genotype_Early$result$p.value 
		) ,
		'LifeStage EvKO effect size'      = list(
			low  = VsplitLow$Genotype_Early$effectSize      ,
			high = VsplitHig$Genotype_Early$effectSize
		)     ,
		#####################################################################
		'LifeStage LvKO estimate'         =	list(
			low  = CatEstimateAndCI(VsplitLow$Genotype_Late$result),
			high = CatEstimateAndCI(VsplitHig$Genotype_Late$result)
		), 
		'LifeStage LvKO standard error'   =	NULL,
		'LifeStage LvKO p-val'            =	lowHighList(
			VsplitLow$Genotype_Late$result$p.value ,
			VsplitHig$Genotype_Late$result$p.value 
		),
		'LifeStage LvKO effect size'      = list(
			low  = VsplitLow$Genotype_Late$effectSize     ,
			high = VsplitHig$Genotype_Late$effectSize
		)    ,
		#####################################################################
		################ Sex LifeStage Genotype interactions
		# 1.
		'LifeStageSexGenotype FvEvKO estimate'        =	list(
			low  = CatEstimateAndCI(VsplitLow$Genotype_Female.Early$result),
			high = CatEstimateAndCI(VsplitHig$Genotype_Female.Early$result)
		), 
		'LifeStageSexGenotype FvEvKO standard error'  =	NULL,
		'LifeStageSexGenotype FvEvKO p-val'           =	lowHighList(
			VsplitLow$Genotype_Female.Early$result$p.value ,
			VsplitHig$Genotype_Female.Early$result$p.value 
		),
		'LifeStageSexGenotype FvEvKO effect size'     = list(
			low  = VsplitLow$Genotype_Female.Early$effectSize    ,
			high = VsplitHig$Genotype_Female.Early$effectSize
		)   ,
		# 2.
		'LifeStageSexGenotype MvEvKO estimate'        =	list(
			low  = CatEstimateAndCI(VsplitLow$Genotype_Male.Early$result),
			high = CatEstimateAndCI(VsplitHig$Genotype_Male.Early$result)
		), 
		'LifeStageSexGenotype MvEvKO standard error'  =	NULL,
		'LifeStageSexGenotype MvEvKO p-val'           =	lowHighList(
			VsplitLow$Genotype_Male.Early$result$p.value ,
			VsplitHig$Genotype_Male.Early$result$p.value 
		),
		'LifeStageSexGenotype MvEvKO effect size'     = list(
			low  = VsplitLow$Genotype_Male.Early$effectSize    ,
			high = VsplitHig$Genotype_Male.Early$effectSize
		)   ,
		# 3.
		'LifeStageSexGenotype FvLvKO estimate'        = list(
			low  = CatEstimateAndCI(VsplitLow$Genotype_Female.Late$result),
			high = CatEstimateAndCI(VsplitHig$Genotype_Female.Late$result)
		), 
		'LifeStageSexGenotype FvLvKO standard error'  = NULL,
		'LifeStageSexGenotype FvLvKO p-val'           = lowHighList(
			VsplitLow$Genotype_Female.Late$result$p.value ,
			VsplitHig$Genotype_Female.Late$result$p.value 
		),
		'LifeStageSexGenotype FvLvKO effect size'     = list(
			low  = VsplitLow$Genotype_Female.Late$effectSize    ,
			high = VsplitHig$Genotype_Female.Late$effectSize
		)   ,
		
		'LifeStageSexGenotype MvLvKO estimate'        = list(
			low  = CatEstimateAndCI(VsplitLow$Genotype_Male.Late$result),
			high = CatEstimateAndCI(VsplitHig$Genotype_Male.Late$result)
		),
		'LifeStageSexGenotype MvLvKO standard error'  = NULL,
		'LifeStageSexGenotype MvLvKO p-val'           =	lowHighList(
			VsplitLow$Genotype_Male.Late$result$p.value,
			VsplitHig$Genotype_Male.Late$result$p.value
		),
		'LifeStageSexGenotype MvLvKO effect size'     = list(
			low  = VsplitLow$Genotype_Male.Late$effectSize ,
			high = VsplitHig$Genotype_Male.Late$effectSize 
		) ,
		################
		'Classification tag'                          =	NA      ,
		'Transformation'                              =	'Legacy',
		'Additional information'                      =	addInfo
	)
	
	return(vectorOutput)
}
