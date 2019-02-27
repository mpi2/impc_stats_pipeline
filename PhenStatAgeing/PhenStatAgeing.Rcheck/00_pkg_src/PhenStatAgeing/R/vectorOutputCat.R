vectorOutputCat =	function(object, json = FALSE)
{
	if (is.null(object))
		return (NULL)
	#####################################################################
	Labels         = PhenListAgeingLevels(object = object)
	Fmodel         = object$extra$Cleanedformula
	frm            = formula(Fmodel)
	depVariable    = all.vars(frm)[1]
	#	equation       = NULL
	formula        = printformula(frm)
	framework      = 'Fisher Exact Test framework'
	#	fittingMethod  = NULL
	#####################################################################
	x                = object$input$data
	columnOfInterest = x[, c(depVariable)]
	#####################################################################
	variability      = length(unique(columnOfInterest)) / length(columnOfInterest)
	#####################################################################
	DSsize            = SummaryStats(
		x = x,
		formula = object$extra$Cleanedformula,
		label = 'Summary statistics',
		lower = TRUE,
		drop = TRUE,
		sep = ' '
	)
	MultiBatch = ifelse(multiBatch(x),
											'Dataset contains multi batches',
											'Dataset contains single batch')
	
	OtherModels        = list(
		'Further model(s)' = setNames(sapply(object$output$SplitModels, function(v) {
			lapply(
				v,
				FUN = function(v2) {
					list('p-value' = v2$result$p.value,
							 'effect size' = v2$effectSize)
				}
			)
		}),nm = names(object$output$SplitModels)),
		'Other residual normality test(s)' = NULL
	)
	addInfo           = c(
		'Variability'            = variability,
		'Multibatch in analysis' = MultiBatch,
		'Formula'                = formula,
		"Gender included in analysis" = ifelse(
			nlevels(x$Sex) > 1,
			'Both sexes included',
			'Only one sex included in the analysis'
		),
		DSsize,
		OtherModels
	)
	#####################################################################
	percentageChanges = NA
	#####################################################################
	vectorOutput      = list(
		'Method'                               = 	framework,
		'Dependent variable'                   =	depVariable,
		'Batch included'                       =	 NULL,
		'Batch p-val'                          =  'Legacy',
		'Residual variances homogeneity'       =   NULL,
		'Residual variances homogeneity p-val' =  'Legacy',
		#####################################################################
		'Genotype contribution' =	list(
			Overal = object$output$SplitModels$category$Genotype$result$p.value,
			'Sex FvKO p-val'   =	object$output$SplitModels$category$Genotype_Female$result$p.value,
			'Sex MvKO p-val'   =  object$output$SplitModels$category$Genotype_Male$result$p.value	,
			'Sexual dimorphism detected' = 'Sex specific results are always reported'
		),
		'Genotype estimate'           = NULL,
		'Genotype standard error'     = NULL,
		'Genotype p-val'              = object$output$SplitModels$category$Genotype$result$p.value ,
		'Genotype percentage change'  =	percentageChanges                                 ,
		'Genotype effect size'        = object$output$SplitModels$category$Genotype$effectSize     ,
		#####################################################################
		'Sex estimate'                =	NULL,
		'Sex standard error'          = NULL,
		'Sex p-val'                   =	object$output$SplitModels$category$Sex$result$p.value,
		'Sex effect size'             =	object$output$SplitModels$category$Sex$effectSize    ,
		#####################################################################
		'LifeStage estimate'          =	NULL,
		'LifeStage standard error'    =	NULL,
		'LifeStage p-val'             =	object$output$SplitModels$category$LifeStage$result$p.value,
		'LifeStage effect size'       = object$output$SplitModels$category$LifeStage$effectSize    ,
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
		'Interaction(s) included'          =	list(
			'Genotype Sex'                   =  NULL,
			'Genotype LifeStage'             =  NULL,
			'Sex LifeStage'                  =  NULL,
			'Genotype Sex LifeStage'         =  NULL
		),
		#####################################################################
		################ interaction
		'Interaction(s) p-val'      =	list(
			'Genotype Sex'            = NULL,
			'Genotype LifeStage'      = NULL,
			'Sex LifeStage'           = NULL,
			'Genotype Sex LifeStage'  = NULL
		),
		#####################################################################
		################ Sex interactions
		'Sex FvKO estimate'         = NULL,
		'Sex FvKO standard error'   = NULL,
		'Sex FvKO p-val'            = object$output$SplitModels$category$Genotype_Female$result$p.value,
		'Sex FvKO effect size'      = object$output$SplitModels$category$Genotype_Female$effectSize,
		#####################################################################
		'Sex MvKO estimate'         = NULL,
		'Sex MvKO standard error'   = NULL,
		'Sex MvKO p-val'            = object$output$SplitModels$category$Genotype_Male$result$p.value,
		'Sex MvKO effect size'      = object$output$SplitModels$category$Genotype_Male$effectSize,
		#####################################################################
		################ LifeStage interaction
		'LifeStage EvKO estimate'         =	NULL,
		'LifeStage EvKO standard error'   =	NULL,
		'LifeStage EvKO p-val'            =	object$output$SplitModels$category$Genotype_Early$result$p.value,
		'LifeStage EvKO effect size'      = object$output$SplitModels$category$Genotype_Early$effectSize,
		#####################################################################
		'LifeStage LvKO estimate'         =	NULL,
		'LifeStage LvKO standard error'   =	NULL,
		'LifeStage LvKO p-val'            =	object$output$SplitModels$category$Genotype_Late$result$p.value,
		'LifeStage LvKO effect size'      = object$output$SplitModels$category$Genotype_Late$effectSize,
		#####################################################################
		################ Sex LifeStage Genotype interactions
		# 1.
		'LifeStageSexGenotype FvEvKO estimate'        =	NULL,
		'LifeStageSexGenotype FvEvKO standard error'  =	NULL,
		'LifeStageSexGenotype FvEvKO p-val'           =	object$output$SplitModels$category$Genotype_Female.Early$result$p.value,
		'LifeStageSexGenotype FvEvKO effect size'     = object$output$SplitModels$category$Genotype_Female.Early$effectSize,
		# 2.
		'LifeStageSexGenotype MvEvKO estimate'        =	NULL,
		'LifeStageSexGenotype MvEvKO standard error'  =	NULL,
		'LifeStageSexGenotype MvEvKO p-val'           =	object$output$SplitModels$category$Genotype_Male.Early$result$p.value,
		'LifeStageSexGenotype MvEvKO effect size'     = object$output$SplitModels$category$Genotype_Male.Early$effectSize,
		# 3.
		'LifeStageSexGenotype FvLvKO estimate'        = NULL,
		'LifeStageSexGenotype FvLvKO standard error'  = NULL,
		'LifeStageSexGenotype FvLvKO p-val'           = object$output$SplitModels$category$Genotype_Female.Late$result$p.value,
		'LifeStageSexGenotype FvLvKO effect size'     = object$output$SplitModels$category$Genotype_Female.Late$effectSize,
		
		'LifeStageSexGenotype MvLvKO estimate'        = NULL,
		'LifeStageSexGenotype MvLvKO standard error'  = NULL,
		'LifeStageSexGenotype MvLvKO p-val'           =	object$output$SplitModels$category$Genotype_Male.Late$result$p.value,
		'LifeStageSexGenotype MvLvKO effect size'     = object$output$SplitModels$category$Genotype_Male.Late$effectSize ,
		################
		'Classification tag'                          =	NA      ,
		'Transformation'                              =	'Legacy',
		'Additional information'                      =	addInfo
	)
	
	if (json) {
		for (i in 1:10) {
			vectorOutput = toJSON(
				vectorOutput,
				auto_unbox = TRUE,
				null = 'null',
				digits = 500
			)
			if (i != 10)
				vectorOutput = fromJSON(txt = vectorOutput)
		}
		write(
			vectorOutput,
			file = 'd:/jsontestCat.txt',
			ncolumns = 10 ^ 3,
			append = FALSE
		)
	}
	
	
	return(vectorOutput)
}
