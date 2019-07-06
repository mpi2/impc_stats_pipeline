vectorOutputCat =	function(object)
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
	framework      = 'Fisher Exact Test framework'
	#	fittingMethod  = NULL
	#####################################################################
	x                = object$input$data
	columnOfInterest = x[, c(depVariable)]
	#####################################################################
	variability      =  list('Value' = length(unique(columnOfInterest)) / max(length(columnOfInterest), 1), 
													 'Type'  = 'Length of unique response divided by total number of response')
	#####################################################################
	DSsize            = SummaryStats(
		x = x,
		formula = object$input$formula,
		#label = 'Summary statistics',
		lower = TRUE,
		drop = TRUE,
		sep = '_'
	)
	MultiBatch = ifelse(multiBatch(x),
											'Dataset contains multi batches',
											'Dataset contains single batch')
	addInfo           = list(
		Data = list(
			'Data signature'         = dataSignature(formula = frm,
																							 data   = x),
			'Variability'            = variability      ,
			'Summary statistics'      = DSsize
		),
		Analysis = list(
			# 'Formula'                = list(
			# 	input   = printformula(object$input$formula),
			# 	final   = printformula(formula)
			# ),
			'Model setting' =  extractFERRTerms(object),
			'Is model optimised'     = NULL                    , 
			'Multibatch in analysis' = MultiBatch,
			'Gender included in analysis' = ifelse(
				nlevels(x$Sex) > 1,
				'Both sexes included',
				paste0('Only one sex included in the analysis; ', levels(x$Sex))
			),
			'Further models' = if (!is.null(object$output$SplitModels)) {
				setNames(sapply(object$output$SplitModels, function(v) {
					lapply(
						v,
						FUN = function(v2) {
							list('P-value'     = v2$result$p.value,
									 'Effect size' = v2$effectSize)
						}
					)
				}),
				nm = names(object$output$SplitModels))
			} else{
				NULL
			},
			'Effect sizes'                   = 'Look at the individual models',
			'Other residual normality tests' = NULL
		)
	)
	#####################################################################
	percentageChanges = NULL
	#####################################################################
	vectorOutput      = list(
		'Applied method'                       = 	framework,
		'Dependent variable'                   =	depVariable,
		'Batch included'                       =	 NULL   ,
		'Batch p-value'                        =   NULL,
		'Residual variances homogeneity'       =   NULL   ,
		'Residual variances homogeneity p-value' =   NULL,
		#####################################################################
		'Genotype contribution' =	list(
			Overal = object$output$SplitModels[[depVariable]]$Genotype$result$p.value,
			'Sex FvKO p-value'   =	object$output$SplitModels[[depVariable]]$Genotype_Female$result$p.value,
			'Sex MvKO p-value'   =  object$output$SplitModels[[depVariable]]$Genotype_Male$result$p.value	,
			'Sexual dimorphism detected' = 'Sex specific results are always reported'
		),
		'Genotype estimate'           = CatEstimateAndCI(object$output$SplitModels[[depVariable]]$Genotype$result),
		'Genotype standard error'     = NULL,
		'Genotype p-value'            = object$output$SplitModels[[depVariable]]$Genotype$result$p.value ,
		'Genotype percentage change'  =	percentageChanges                                 ,
		'Genotype effect size'        = object$output$SplitModels[[depVariable]]$Genotype$effectSize     ,
		#####################################################################
		'Sex estimate'                =	CatEstimateAndCI(object$output$SplitModels[[depVariable]]$Sex$result),
		'Sex standard error'          = NULL,
		'Sex p-value'                 =	object$output$SplitModels[[depVariable]]$Sex$result$p.value,
		'Sex effect size'             =	object$output$SplitModels[[depVariable]]$Sex$effectSize    ,
		#####################################################################
		'LifeStage estimate'          =	CatEstimateAndCI(object$output$SplitModels[[depVariable]]$LifeStage$result),
		'LifeStage standard error'    =	NULL,
		'LifeStage p-value'           =	object$output$SplitModels[[depVariable]]$LifeStage$result$p.value,
		'LifeStage effect size'       = object$output$SplitModels[[depVariable]]$LifeStage$effectSize    ,
		#####################################################################
		'Weight estimate'             =	NULL,
		'Weight standard error'       =	NULL,
		'Weight p-value'              =	NULL,
		'Weight effect size'          = NULL,
		#####################################################################
		'Gp1 genotype'                     =	Labels$Genotype$Control		,
		'Gp1 Residuals normality test'     =	NULL                      ,
		'Gp2 genotype'                     =	Labels$Genotype$Mutant		,
		'Gp2 Residuals normality test'     =	NULL                      ,
		#####################################################################
		'Blups test'                       =  NULL,
		'Rotated residuals normality test' =  NULL,
		#####################################################################
		'Intercept estimate'               =	NULL,
		'Intercept standard error'         =	NULL,
		'Intercept p-value'                =	NULL,
		#####################################################################
		'Interactions included'          =	list(
			'Genotype Sex'                   =  NULL,
			'Genotype LifeStage'             =  NULL,
			'Sex LifeStage'                  =  NULL,
			'Genotype Sex LifeStage'         =  NULL
		),
		#####################################################################
		################ interaction
		'Interactions p-value'      =	list(
			'Genotype Sex'            = NULL,
			'Genotype LifeStage'      = NULL,
			'Sex LifeStage'           = NULL,
			'Genotype Sex LifeStage'  = NULL
		),
		#####################################################################
		################ Sex interactions
		'Sex FvKO estimate'         = CatEstimateAndCI(object$output$SplitModels[[depVariable]]$Genotype_Female$result),
		'Sex FvKO standard error'   = NULL,
		'Sex FvKO p-value'            = object$output$SplitModels[[depVariable]]$Genotype_Female$result$p.value,
		'Sex FvKO effect size'      = object$output$SplitModels[[depVariable]]$Genotype_Female$effectSize,
		#####################################################################
		'Sex MvKO estimate'         = CatEstimateAndCI(object$output$SplitModels[[depVariable]]$Genotype_Male$result),
		'Sex MvKO standard error'   = NULL,
		'Sex MvKO p-value'          = object$output$SplitModels[[depVariable]]$Genotype_Male$result$p.value,
		'Sex MvKO effect size'      = object$output$SplitModels[[depVariable]]$Genotype_Male$effectSize,
		#####################################################################
		################ LifeStage interaction
		'LifeStage EvKO estimate'         =	CatEstimateAndCI(object$output$SplitModels[[depVariable]]$Genotype_Early$result),
		'LifeStage EvKO standard error'   =	NULL,
		'LifeStage EvKO p-value'          =	object$output$SplitModels[[depVariable]]$Genotype_Early$result$p.value,
		'LifeStage EvKO effect size'      = object$output$SplitModels[[depVariable]]$Genotype_Early$effectSize,
		#####################################################################
		'LifeStage LvKO estimate'         =	CatEstimateAndCI(object$output$SplitModels[[depVariable]]$Genotype_Late$result),
		'LifeStage LvKO standard error'   =	NULL,
		'LifeStage LvKO p-value'          =	object$output$SplitModels[[depVariable]]$Genotype_Late$result$p.value,
		'LifeStage LvKO effect size'      = object$output$SplitModels[[depVariable]]$Genotype_Late$effectSize,
		#####################################################################
		################ Sex LifeStage Genotype interactions
		# 1.
		'LifeStageSexGenotype FvEvKO estimate'        =	CatEstimateAndCI(object$output$SplitModels[[depVariable]]$Genotype_Female.Early$result),
		'LifeStageSexGenotype FvEvKO standard error'  =	NULL,
		'LifeStageSexGenotype FvEvKO p-value'         =	object$output$SplitModels[[depVariable]]$Genotype_Female.Early$result$p.value,
		'LifeStageSexGenotype FvEvKO effect size'     = object$output$SplitModels[[depVariable]]$Genotype_Female.Early$effectSize,
		# 2.
		'LifeStageSexGenotype MvEvKO estimate'        =	CatEstimateAndCI(object$output$SplitModels[[depVariable]]$Genotype_Male.Early$result),
		'LifeStageSexGenotype MvEvKO standard error'  =	NULL,
		'LifeStageSexGenotype MvEvKO p-value'         =	object$output$SplitModels[[depVariable]]$Genotype_Male.Early$result$p.value,
		'LifeStageSexGenotype MvEvKO effect size'     = object$output$SplitModels[[depVariable]]$Genotype_Male.Early$effectSize,
		# 3.
		'LifeStageSexGenotype FvLvKO estimate'        = CatEstimateAndCI(object$output$SplitModels[[depVariable]]$Genotype_Female.Late$result),
		'LifeStageSexGenotype FvLvKO standard error'  = NULL,
		'LifeStageSexGenotype FvLvKO p-value'         = object$output$SplitModels[[depVariable]]$Genotype_Female.Late$result$p.value,
		'LifeStageSexGenotype FvLvKO effect size'     = object$output$SplitModels[[depVariable]]$Genotype_Female.Late$effectSize,
		
		'LifeStageSexGenotype MvLvKO estimate'        = CatEstimateAndCI(object$output$SplitModels[[depVariable]]$Genotype_Male.Late$result),
		'LifeStageSexGenotype MvLvKO standard error'  = NULL,
		'LifeStageSexGenotype MvLvKO p-value'         =	object$output$SplitModels[[depVariable]]$Genotype_Male.Late$result$p.value,
		'LifeStageSexGenotype MvLvKO effect size'     = object$output$SplitModels[[depVariable]]$Genotype_Male.Late$effectSize ,
		################
		'Classification tag'                          =	NULL,
		'Transformation'                              =	NULL,
		'Additional information'                      =	addInfo
	)
	return(vectorOutput)
}
