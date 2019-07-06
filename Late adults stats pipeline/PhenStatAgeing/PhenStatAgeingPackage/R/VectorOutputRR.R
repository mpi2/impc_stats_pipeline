vectorOutputRR =	function(object)
{
	if (!is.null(object$messages))
		return (NULL)
	#####################################################################
	Labels         = PhenListAgeingLevels(object = object)
	Fmodel         = object$extra$Cleanedformula
	frm            = formula(Fmodel)
	depVariable    = all_vars0(frm)[1]
	refVariable    = all_vars0(frm)[2]
	#	equation     = NULL
	formula        = printformula(frm)
	framework      = paste0('Reference Range Plus Test framework; quantile = ',
													object$input$prop)
	#	fittingMethod  = NULL
	#####################################################################
	Vsplit         = object$output$SplitModels
	low            = pastedot('Low' , depVariable, refVariable)
	high           = pastedot('High', depVariable, refVariable)
	VsplitLow      = Vsplit[[low ]]
	VsplitHig      = Vsplit[[high]]
	#### Sex
	VsplitLowFemale= Vsplit[[pastedot(low ,'Female')]]
	VsplitHigFemale= Vsplit[[pastedot(high,'Female')]]
	VsplitLowMale  = Vsplit[[pastedot(low ,'Male'  )]]
	VsplitHigMale  = Vsplit[[pastedot(high,'Male'  )]]
	### Lifestage
	VsplitLowEarly = Vsplit[[pastedot(low ,'Early'  )]]
	VsplitHigEarly = Vsplit[[pastedot(high,'Early'  )]]
	VsplitLowLate  = Vsplit[[pastedot(low ,'Late'   )]]
	VsplitHigLate  = Vsplit[[pastedot(high,'Late'   )]]
	### LifeStage x Sex
	# Female
	VsplitLowEarlyFemale = Vsplit[[pastedot(low ,'Female','Early'  )]]
	VsplitHigEarlyFemale = Vsplit[[pastedot(high,'Female','Early'  )]]
	VsplitLowLateFemale  = Vsplit[[pastedot(low ,'Female','Late'   )]]
	VsplitHigLateFemale  = Vsplit[[pastedot(high,'Female','Late'   )]]
	# Male
	VsplitLowEarlyMale = Vsplit[[pastedot(low ,'Male','Early'  )]]
	VsplitHigEarlyMale = Vsplit[[pastedot(high,'Male','Early'  )]]
	VsplitLowLateMale  = Vsplit[[pastedot(low ,'Male','Late'   )]]
	VsplitHigLateMale  = Vsplit[[pastedot(high,'Male','Late'   )]]
	###
	GenotypeDiscLabel  = 'Data is discritised by Genotype levels'
	SexDiscLabel       = 'Data is first split by Sex levels then discritised by Genotype levels'
	LifeStageDiscLabel = 'Data is first split by LifeStage levels then discritised by Genotype levels'
	LifeStageSexDiscLabel = 'Data is first split by LifeStage-Sex levels then discritised by Genotype levels'
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
		lower  = TRUE,
		drop   = TRUE,
		sep    = '_'
	)
	MultiBatch = ifelse(multiBatch(x),
											'Dataset contains multi batches',
											'Dataset contains single batch')
	addInfo           = list(
		Data = list(
			'Data signature'          = dataSignature(formula = frm,
																								data   = x)  ,
			'Variability'             = variability                 ,
			'Summary statistics'      = DSsize
		),
		Analysis = list(
			'Model setting'          = extractFERRTerms(object),
			'Is model optimised'     = NULL                    , 
			'Multibatch in analysis' = MultiBatch,
			'Gender included in analysis' = ifelse(
				nlevels(x$Sex) > 1,
				'Both sexes included',
				paste0('Only one sex included in the analysis; ', levels(x$Sex))
			),
			'Further models' = if (!is.null(Vsplit)) {
				setNames(lapply(Vsplit, function(v) {
					lapply(
						v,
						FUN = function(v2) {
							list('P-value'     = v2$result$p.value,
									 'Effect size' = v2$effectSize    ,
									 'Extra'       = v2$RRextra)
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
	percentageChanges = NULL
	#####################################################################
	vectorOutput      = list(
		'Applied method'                       = 	framework   ,
		'Dependent variable'                   =	depVariable , 
		'Batch included'                       =	 NULL       ,
		'Batch p-value'                        =   NULL       ,
		'Residual variances homogeneity'       =   NULL       ,
		'Residual variances homogeneity p-value' =   NULL       ,
		#####################################################################
		'Genotype contribution' =	list(
			Overal = lowHighList(
				VsplitLow$Genotype$result$p.value,
				VsplitHig$Genotype$result$p.value,
				'Details'  = GenotypeDiscLabel
			),
			'Sex FvKO p-value'   =	lowHighList(
				VsplitLowFemale$Genotype$result$p.value,
				VsplitHigFemale$Genotype$result$p.value,
				'Details' = SexDiscLabel
			),
			'Sex MvKO p-value'   =  lowHighList(
				VsplitLowMale$Genotype$result$p.value,
				VsplitHigMale$Genotype$result$p.value,
				'Details' = SexDiscLabel
			),
			'Sexual dimorphism detected' = 'Sex specific results for Low/High tables are always reported'
		),
		'Genotype estimate'           = lowHighList(
			CatEstimateAndCI(VsplitLow$Genotype$result),
			CatEstimateAndCI(VsplitHig$Genotype$result),
			'Details' = GenotypeDiscLabel
		), 
		'Genotype standard error'     = NULL, 
		'Genotype p-value'            = lowHighList(
			VsplitLow$Genotype$result$p.value,
			VsplitHig$Genotype$result$p.value,
			'Details' = GenotypeDiscLabel
		),
		'Genotype percentage change'  =	percentageChanges                            ,
		'Genotype effect size'        = lowHighList(VsplitLow$Genotype$effectSize    ,
																								VsplitHig$Genotype$effectSize    ,
																								'Details' = GenotypeDiscLabel)     ,
		#####################################################################
		'Sex estimate'                =	lowHighList(
			CatEstimateAndCI(VsplitLow$Sex$result),
			CatEstimateAndCI(VsplitHig$Sex$result),
			'Details' = GenotypeDiscLabel
		), 
		'Sex standard error'          = NULL,
		'Sex p-value'                 =	lowHighList(
			VsplitLow$Sex$result$p.value      ,
			VsplitHig$Sex$result$p.value      ,
			'Details' = GenotypeDiscLabel
		)   ,
		'Sex effect size'             =	lowHighList(VsplitLow$Sex$effectSize    ,
																								VsplitHig$Sex$effectSize    ,
																								'Details' = GenotypeDiscLabel), 
		#####################################################################
		'LifeStage estimate'          =	lowHighList(
			CatEstimateAndCI(VsplitLow$LifeStage$result),
			CatEstimateAndCI(VsplitHig$LifeStage$result),
			'Details' = GenotypeDiscLabel
		), 
		'LifeStage standard error'    =	NULL,
		'LifeStage p-value'           =	lowHighList(
			VsplitLow$LifeStage$result$p.value,
			VsplitHig$LifeStage$result$p.value,
			'Details' = GenotypeDiscLabel
		),
		'LifeStage effect size'       = lowHighList(
			VsplitLow$LifeStage$effectSize     ,
			VsplitHig$LifeStage$effectSize     ,
			'Details' = GenotypeDiscLabel
		)    ,
		#####################################################################
		'Weight estimate'                  =	NULL,
		'Weight standard error'            =	NULL,
		'Weight p-value'                   =	NULL,
		'Weight effect size'               =  NULL,
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
		'Interactions p-value'            =	list(
			'Genotype Sex'                  = NULL,
			'Genotype LifeStage'            = NULL,
			'Sex LifeStage'                 = NULL,
			'Genotype Sex LifeStage'        = NULL
		),
		#####################################################################
		################ Sex interactions
		'Sex FvKO estimate'         = lowHighList(
			CatEstimateAndCI(VsplitLowFemale$Genotype$result),
			CatEstimateAndCI(VsplitHigFemale$Genotype$result),
			'Details' = SexDiscLabel
		), 
		'Sex FvKO standard error'   = NULL,
		'Sex FvKO p-value'          = lowHighList(
			VsplitLowFemale$Genotype$result$p.value,
			VsplitHigFemale$Genotype$result$p.value,
			'Details' = SexDiscLabel
		),
		'Sex FvKO effect size'      = lowHighList(
			VsplitLowFemale$Genotype$effectSize,
			VsplitHigFemale$Genotype$effectSize,
			'Details' = SexDiscLabel
		),
		#####################################################################
		'Sex MvKO estimate'         = lowHighList(
			CatEstimateAndCI(VsplitLowMale$Genotype$result),
			CatEstimateAndCI(VsplitHigMale$Genotype$result),
			'Details' = SexDiscLabel
		), 
		'Sex MvKO standard error'   = NULL,
		'Sex MvKO p-value'          = lowHighList(
			VsplitLowMale$Genotype$result$p.value,
			VsplitHigMale$Genotype$result$p.value,
			'Details' = SexDiscLabel
		)  ,
		'Sex MvKO effect size'      = lowHighList(
			VsplitLowMale$Genotype$effectSize,
			VsplitHigMale$Genotype$effectSize,
			'Details' = SexDiscLabel
		) ,
		#####################################################################
		################ LifeStage interaction
		'LifeStage EvKO estimate'         =	lowHighList(
			CatEstimateAndCI(VsplitLowEarly$Genotype$result),
			CatEstimateAndCI(VsplitHigEarly$Genotype$result),
			'Details' = LifeStageDiscLabel
		), 
		'LifeStage EvKO standard error'   =	NULL,
		'LifeStage EvKO p-value'          =	lowHighList(
			VsplitLowEarly$Genotype$result$p.value,
			VsplitHigEarly$Genotype$result$p.value,
			'Details' = LifeStageDiscLabel
		) ,
		'LifeStage EvKO effect size'      = lowHighList(
			VsplitLowEarly$Genotype$effectSize,
			VsplitHigEarly$Genotype$effectSize,
			'Details' = LifeStageDiscLabel
		)  ,
		#####################################################################
		'LifeStage LvKO estimate'         =	lowHighList(
			CatEstimateAndCI(VsplitLowLate$Genotype$result),
			CatEstimateAndCI(VsplitHigLate$Genotype$result),
			'Details' = LifeStageDiscLabel
		), 
		'LifeStage LvKO standard error'   =	NULL,
		'LifeStage LvKO p-value'          =	lowHighList(
			VsplitLowLate$Genotype$result$p.value ,
			VsplitHigLate$Genotype$result$p.value ,
			'Details' = LifeStageDiscLabel
		),
		'LifeStage LvKO effect size'      = lowHighList(
			VsplitLowLate$Genotype$effectSize,
			VsplitHigLate$Genotype$effectSize,
			'Details' = LifeStageDiscLabel
		),
		#####################################################################
		################ Sex LifeStage Genotype interactions
		#####################################################################
		# 1.
		'LifeStageSexGenotype FvEvKO estimate'        =	lowHighList(
			CatEstimateAndCI(VsplitLowEarlyFemale$Genotype$result),
			CatEstimateAndCI(VsplitHigEarlyFemale$Genotype$result),
			'Details' = LifeStageSexDiscLabel
		), 
		'LifeStageSexGenotype FvEvKO standard error'  =	NULL,
		'LifeStageSexGenotype FvEvKO p-value'         =	lowHighList(
			VsplitLowEarlyFemale$Genotype$result$p.value,
			VsplitHigEarlyFemale$Genotype$result$p.value,
			'Details' = LifeStageSexDiscLabel
		),
		'LifeStageSexGenotype FvEvKO effect size'     = lowHighList(
			VsplitLowEarlyFemale$Genotype$effectSize,
			VsplitHigEarlyFemale$Genotype$effectSize,
			'Details' = LifeStageSexDiscLabel
		)   ,
		# 2.
		'LifeStageSexGenotype MvEvKO estimate'        =	lowHighList(
			CatEstimateAndCI(VsplitLowEarlyMale$Genotype$result),
			CatEstimateAndCI(VsplitHigEarlyMale$Genotype$result),
			'Details' = LifeStageSexDiscLabel
		), 
		'LifeStageSexGenotype MvEvKO standard error'  =	NULL,
		'LifeStageSexGenotype MvEvKO p-value'         =	lowHighList(
			VsplitLowEarlyMale$Genotype$result$p.value,
			VsplitHigEarlyMale$Genotype$result$p.value,
			'Details' = LifeStageSexDiscLabel
		),
		'LifeStageSexGenotype MvEvKO effect size'     = lowHighList(
			VsplitLowEarlyMale$Genotype$effectSize,
			VsplitHigEarlyMale$Genotype$effectSize,
			'Details' = LifeStageSexDiscLabel
		),
		# 3.
		'LifeStageSexGenotype FvLvKO estimate'        = lowHighList(
			CatEstimateAndCI(VsplitLowLateFemale$Genotype$result),
			CatEstimateAndCI(VsplitHigLateFemale$Genotype$result),
			'Details' = LifeStageSexDiscLabel
		), 
		'LifeStageSexGenotype FvLvKO standard error'  = NULL,
		'LifeStageSexGenotype FvLvKO p-value'         = lowHighList(
			VsplitLowLateFemale$Genotype$result$p.value ,
			VsplitHigLateFemale$Genotype$result$p.value ,
			'Details' = LifeStageSexDiscLabel
		),
		'LifeStageSexGenotype FvLvKO effect size'     = lowHighList(
			VsplitLowLateFemale$Genotype$effectSize,
			VsplitHigLateFemale$Genotype$effectSize,
			'Details' = LifeStageSexDiscLabel
		)   ,
		'LifeStageSexGenotype MvLvKO estimate'        = lowHighList(
			CatEstimateAndCI(VsplitLowLateMale$Genotype$result),
			CatEstimateAndCI(VsplitHigLateMale$Genotype$result),
			'Details' = LifeStageSexDiscLabel
		),
		'LifeStageSexGenotype MvLvKO standard error'  = NULL,
		'LifeStageSexGenotype MvLvKO p-value'         =	lowHighList(
			VsplitLowLateMale$Genotype$result$p.value,
			VsplitHigLateMale$Genotype$result$p.value,
			'Details' = LifeStageSexDiscLabel
		),
		'LifeStageSexGenotype MvLvKO effect size'     = lowHighList(
			VsplitLowLateMale$Genotype$effectSize,
			VsplitHigLateMale$Genotype$effectSize,
			'Details' = LifeStageSexDiscLabel
		),
		################
		'Classification tag'                          =	NULL,
		'Transformation'                              =	NULL,
		'Additional information'                      =	addInfo
	)
	
	return(vectorOutput)
}
