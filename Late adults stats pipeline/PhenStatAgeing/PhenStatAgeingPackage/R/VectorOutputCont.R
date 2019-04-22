vectorOutputCont =	function(object,
														debug = FALSE)
{
	if (!is.null(object$messages))
		return (NULL)
	#####################################################################
	Labels         = PhenListAgeingLevels(object = object)
	Fmodel         = object$output$Final.Model
	frm            = formula(Fmodel)
	depVariable    = all.vars(frm)[1]
	equation       = ifelse(
		Labels$Weight %in% all.vars(frm),
		paste0('equation with'   , Labels$Weight),
		paste0('equation without', Labels$Weight)
	)
	formula        = printformula(frm)
	#modelContrast  = modelContrasts(formula = frm,data = object$input$data)
	framework      = switch(
		object$output$Final.Model.Tag                                  ,
		lme  = "Linear Mixed Model framework"                          ,
		gls  = "Linear Model Using Generalized Least Squares framework",
		glm  = "Generalized Linear Model framework"
	)
	fittingMethod    = toupper(object$output$Final.Model.Tag)
	#####################################################################
	x                = object$input$PhenListAgeing@datasetPL
	columnOfInterest = x[, c(depVariable)]
	#####################################################################
	variability      = length(unique(columnOfInterest)) / length(columnOfInterest)
	#####################################################################
	DSsize            = SummaryStats(
		x = object$input$data         ,
		formula = checkModelTermsInData(
			formula = object$input$fixed,
			data = x,
			responseIsTheFirst = TRUE
		),
		#label = 'Summary statistics',
		lower = TRUE,
		drop = TRUE,
		sep = ' '
	)
	MultiBatch = ifelse(multiBatch(x),
											'Dataset contains multi batches',
											'Dataset contains single batch')
	addInfo           = list(
		'Formula'                = list(
			input   = printformula(object$input$fixed),
			final   =	printformula(formula)
		),
		'Variability'            = variability             ,
		'Multibatch in analysis' = MultiBatch              ,
		'Gender included in analysis' = ifelse(
			nlevels(x$Sex) > 1                               ,
			'Both sexes included'                            ,
			'Only one sex included in the analysis'
		)                                                  ,
		'Is model optimised'  = object$output$optimised    ,
		'data code'           = dataCode(formula = object$input$fixed, 
																		 data    = object$input$data),
		'Summary statistics'  = DSsize                     ,
		'Further models' = if (!is.null(object$output$SplitModels)) {
			lapply(object$output$SplitModels, function(v) {
				if (class(v) %in% c('lme', 'gls', 'glm')) {
					as.list(unmatrix0(summary(v)$tTable))
				} else{
					v
				}
			})
		} else{
			NULL
		},
		'Effect sizes'                   = object$output$EffectSizes,
		'Other residual normality tests' = object$output$ResidualNormalityTests
	)
	#####################################################################
	pcS = object$output$EffectSizes$CombinedEffectSizes.Genotype_Sex$percentageChange
	pcO = object$output$EffectSizes$Genotype$percentageChange
	percentageChanges = if (!is.null(pcS)) {
		pcS
	} else{
		pcO
	}
	#####################################################################
	vectorOutput      = list(
		'Method'                               = 	paste0(framework, ", ", fittingMethod, ', ', format(equation)),
		'Dependent variable'                   =	depVariable,
		'Batch included'                       =	object$output$BatchIn,
		'Batch p-val'                          =  'Legacy',
		'Residual variances homogeneity'       =	object$output$VarHomoIn,
		'Residual variances homogeneity p-val' =  'Legacy',
		#####################################################################
		'Genotype contribution' =	list(
			Overal = TermInFormulaReturn(
				active = TRUE,
				formula = frm,
				term = CombineLevels(Labels$Genotype$Genotype,Labels$Sex$Sex,debug = debug),
				return = NULL,
				not = modelSummaryPvalueExtract(
					x = Fmodel,
					variable = Labels$Genotype$Genotype,
					anova = TRUE,
					debug = debug
				),
				debug = debug
			),
			'Sex FvKO p-val'   =	TermInFormulaReturn(
				active = TRUE,
				formula = frm,
				term = Labels$Sex$Sex,
				not = NULL,
				return = modelSummaryPvalueExtract(
					x = object$output$SplitModels$Genotype_Sex,
					# SexFemale:Genotypeexperimental
					variable = CombineLevels(
						paste0('Sex', Labels$Sex$Female),
						Labels$Genotype$Levels,
						debug = debug
					),
					anova = FALSE,
					debug = debug
				),
				debug = debug
			),
			'Sex MvKO p-val'  =	TermInFormulaReturn(
				active = TRUE,
				formula = frm,
				term = Labels$Sex$Sex,
				not = NULL,
				return = modelSummaryPvalueExtract(
					x = object$output$SplitModels$Genotype_Sex,
					variable = CombineLevels(
						paste0('Sex', Labels$Sex$Male),
						Labels$Genotype$Levels,
						debug = debug
					),
					anova = FALSE,
					debug = debug
				),
				debug = debug
			),
			'Sexual dimorphism detected' = TermInFormulaReturn(
				active = TRUE,
				formula = frm,
				term = CombineLevels(Labels$Genotype$Genotype,Labels$Sex$Sex,debug = debug),
				return = TRUE,
				not = FALSE,
				debug = debug
			)
		),
		'Genotype estimate' =
			modelSummaryPvalueExtract(
				x = Fmodel,
				variable = unlist(Labels$Genotype$Levels),
				anova = FALSE,
				what = 'Value',
				debug = debug
			),
		'Genotype standard error'  =
			modelSummaryPvalueExtract(
				x = Fmodel,
				variable = unlist(Labels$Genotype$Levels),
				anova = FALSE,
				what = 'Std.Error',
				debug = debug
			),
		'Genotype p-val'           =
			modelSummaryPvalueExtract(
				x = Fmodel,
				variable = Labels$Genotype$Genotype,
				anova = TRUE,
				debug = debug
			),
		'Genotype percentage change'           =	percentageChanges,
		'Genotype effect size'                 = object$output$EffectSizes[[Labels$Genotype$Genotype]],
		#####################################################################
		'Sex estimate'                         =	modelSummaryPvalueExtract(
			x = Fmodel,
			variable = unlist(Labels$Sex$Levels),
			anova = FALSE,
			what = 'Value',
			debug = debug
		),
		'Sex standard error'                   = modelSummaryPvalueExtract(
			x = Fmodel,
			variable = unlist(Labels$Sex$Levels),
			anova = FALSE,
			what = 'Std.Error',
			debug = debug
		),
		'Sex p-val'                            =	modelSummaryPvalueExtract(
			x = Fmodel,
			variable = Labels$Sex$Sex,
			anova = TRUE,
			debug = debug
		),
		'Sex effect size'                         =	object$output$EffectSizes[[Labels$Sex$Sex]],
		#####################################################################
		'LifeStage estimate'                      =	modelSummaryPvalueExtract(
			x = Fmodel,
			variable = unlist(Labels$LifeStage$Levels),
			anova = FALSE,
			what = 'Value',
			debug = debug
		),
		'LifeStage standard error'                =	modelSummaryPvalueExtract(
			x = Fmodel,
			variable = unlist(Labels$LifeStage$Levels),
			anova = FALSE,
			what = 'Std.Error',
			debug = debug
		),
		'LifeStage p-val'                         =	modelSummaryPvalueExtract(
			x = Fmodel,
			variable = Labels$LifeStage$LifeStage,
			anova = TRUE,
			debug = debug
		),
		'LifeStage effect size'                = object$output$EffectSizes[[Labels$LifeStage$LifeStage]] ,
		#####################################################################
		'Weight estimate'                      =	modelSummaryPvalueExtract(
			x = Fmodel,
			variable = Labels$Weight,
			anova = FALSE,
			what = 'Value',
			debug = debug
		),
		'Weight standard error'                =	modelSummaryPvalueExtract(
			x = Fmodel,
			variable = Labels$Weight,
			anova = FALSE,
			what = 'Std.Error',
			debug = debug
		),
		'Weight p-val'                         =	modelSummaryPvalueExtract(
			x = Fmodel,
			variable = Labels$Weight,
			anova = TRUE,
			debug = debug
		),
		'Weight effect size'                   =  object$output$EffectSizes[[Labels$Weight]],
		#####################################################################
		'Gp1 genotype'                         =	Labels$Genotype$Control		,
		'Gp1 Residuals normality test'         =	object$output$ResidualNormalityTests$Genotype[Labels$Genotype$Control],
		'Gp2 genotype'                         =	Labels$Genotype$Mutant				,
		'Gp2 Residuals normality test'         =	object$output$ResidualNormalityTests$Genotype[Labels$Genotype$Mutant],
		#####################################################################
		'Blups test'                           =  'Legacy',
		'Rotated residuals normality test'     =  'Legacy',
		#####################################################################
		'Intercept estimate'                   =	modelSummaryPvalueExtract(
			x = Fmodel,
			variable = '(Intercept)',
			anova = FALSE,
			what = 'Value',
			debug = debug
		),
		'Intercept standard error'             =	modelSummaryPvalueExtract(
			x = Fmodel,
			variable = '(Intercept)',
			anova = FALSE,
			what = 'Std.Error',
			debug = debug
		),
		'Intercept p-val'                      =	modelSummaryPvalueExtract(
			x = Fmodel,
			variable = '(Intercept)',
			anova = TRUE,
			debug = debug
		),
		#####################################################################
		'Interactions included'              =	list(
			'Genotype Sex'  =  !is.null(
				modelSummaryPvalueExtract(
					x = Fmodel,
					variable = CombineLevels(Labels$Genotype$Genotype, Labels$Sex$Sex,debug = debug),
					anova = TRUE,
					debug = debug
				)
			),
			'Genotype LifeStage'  =  !is.null(
				modelSummaryPvalueExtract(
					x = Fmodel,
					variable = CombineLevels(Labels$Genotype$Genotype, Labels$LifeStage$LifeStage,debug = debug),
					anova = TRUE,
					debug = debug
				)
			),
			'Sex LifeStage'  =  !is.null(
				modelSummaryPvalueExtract(
					x = Fmodel,
					variable = CombineLevels(Labels$Sex$Sex, Labels$LifeStage$LifeStage,debug = debug),
					anova = TRUE,
					debug = debug
				)
			),
			'Genotype Sex LifeStage'  =  !is.null(
				modelSummaryPvalueExtract(
					x = Fmodel,
					variable = CombineLevels(
						Labels$Genotype$Genotype,
						Labels$Sex$Sex,
						Labels$LifeStage$LifeStage,
						debug = debug
					),
					anova = TRUE,
					debug = debug
				)
			)
		),
		#####################################################################
		################ interaction
		'Interactions p-val'      =	list(
			'Genotype Sex'  =
				modelSummaryPvalueExtract(
					x = Fmodel,
					variable =  CombineLevels(Labels$Genotype$Genotype, Labels$Sex$Sex,debug = debug),
					anova = TRUE,
					debug = debug
				),
			'Genotype LifeStage'  =
				modelSummaryPvalueExtract(
					x = Fmodel,
					variable =  CombineLevels(Labels$Genotype$Genotype, Labels$LifeStage$LifeStage,debug = debug),
					anova = TRUE,
					debug = debug
				),
			'Sex LifeStage'  =
				modelSummaryPvalueExtract(
					x = Fmodel,
					variable =  CombineLevels(Labels$Sex$Sex, Labels$LifeStage$LifeStage,debug = debug),
					anova = TRUE,
					debug = debug
				),
			'Genotype Sex LifeStage'  =
				modelSummaryPvalueExtract(
					x = Fmodel,
					variable =  CombineLevels(
						Labels$Genotype$Genotype,
						Labels$Sex$Sex,
						Labels$LifeStage$LifeStage,
						debug = debug
					),
					anova = TRUE,
					debug = debug
				)
		),
		#####################################################################
		################ Sex interactions
		'Sex FvKO estimate'                    =
			modelSummaryPvalueExtract(
				x = object$output$SplitModels$Genotype_Sex,
				variable = CombineLevels(
					paste0('Sex', Labels$Sex$Female),
					Labels$Genotype$Levels,
					debug = debug
				),
				anova = FALSE,
				what = 'Value',
				debug = debug
			),
		'Sex FvKO standard error'              =	modelSummaryPvalueExtract(
			x = object$output$SplitModels$Genotype_Sex,
			variable = CombineLevels(
				paste0('Sex', Labels$Sex$Female),
				Labels$Genotype$Levels,
				debug = debug
			),
			anova = FALSE,
			what = 'Std.Error',
			debug = debug
		),
		'Sex FvKO p-val'                       = modelSummaryPvalueExtract(
			x = object$output$SplitModels$Genotype_Sex,
			variable = CombineLevels(
				paste0('Sex', Labels$Sex$Female),
				Labels$Genotype$Levels,
				debug = debug
			),
			anova = FALSE,
			debug = debug
		),
		'Sex FvKO effect size'                 = object$output$EffectSizes[[paste(Labels$Genotype$Genotype, Labels$Sex$Female, sep = '_')]],
		#####################################################################
		'Sex MvKO estimate'                    =	modelSummaryPvalueExtract(
			x = object$output$SplitModels$Genotype_Sex,
			variable = CombineLevels(
				paste0('Sex', Labels$Sex$Male),
				Labels$Genotype$Levels,
				debug = debug
			),
			anova = FALSE,
			what = 'Value',
			debug = debug
		),
		'Sex MvKO standard error'              =	 modelSummaryPvalueExtract(
			x = object$output$SplitModels$Genotype_Sex,
			variable = CombineLevels(
				paste0('Sex', Labels$Sex$Male),
				Labels$Genotype$Levels,
				debug = debug
			),
			anova = FALSE,
			what = 'Std.Error',
			debug = debug
		),
		'Sex MvKO p-val'                       =	modelSummaryPvalueExtract(
			x = object$output$SplitModels$Genotype_Sex,
			variable = CombineLevels(
				paste0('Sex', Labels$Sex$Male),
				Labels$Genotype$Levels,
				debug = debug
			),
			anova = FALSE,
			debug = debug
		),
		'Sex MvKO effect size'                 = object$output$EffectSizes[[paste(Labels$Genotype$Genotype, Labels$Sex$Male, sep = '_')]],
		#####################################################################
		################ LifeStage interaction
		'LifeStage EvKO estimate'                    =	modelSummaryPvalueExtract(
			x = object$output$SplitModels$Genotype_LifeStage,
			variable = CombineLevels(
				paste0(Labels$LifeStage$LifeStage, Labels$LifeStage$Early),
				Labels$Genotype$Levels,
				debug = debug
			),
			anova = FALSE,
			what = 'Value',
			debug = debug
		),
		'LifeStage EvKO standard error'              =	modelSummaryPvalueExtract(
			x = object$output$SplitModels$Genotype_LifeStage,
			variable =  CombineLevels(
				paste0(Labels$LifeStage$LifeStage, Labels$LifeStage$Early),
				Labels$Genotype$Levels,
				debug = debug
			),
			anova = FALSE,
			what = 'Std.Error',
			debug = debug
		),
		'LifeStage EvKO p-val'                       =	modelSummaryPvalueExtract(
			x = object$output$SplitModels$Genotype_LifeStage,
			variable =  CombineLevels(
				paste0(Labels$LifeStage$LifeStage, Labels$LifeStage$Early),
				Labels$Genotype$Levels,
				debug = debug
			),
			anova = FALSE,
			debug = debug
		),
		'LifeStage EvKO effect size'                 = object$output$EffectSizes$Genotype_Early ,
		#####################################################################
		'LifeStage LvKO estimate'                    =	modelSummaryPvalueExtract(
			x = object$output$SplitModels$Genotype_LifeStage,
			variable =  CombineLevels(
				paste0(Labels$LifeStage$LifeStage, Labels$LifeStage$Late),
				Labels$Genotype$Levels,
				debug = debug
			),
			anova = FALSE,
			what = 'Value',
			debug = debug
		),
		'LifeStage LvKO standard error'              =	modelSummaryPvalueExtract(
			x = object$output$SplitModels$Genotype_LifeStage,
			variable = CombineLevels(
				paste0(Labels$LifeStage$LifeStage, Labels$LifeStage$Late),
				Labels$Genotype$Levels,
				debug = debug
			),
			anova = FALSE,
			what = 'Std.Error',
			debug = debug
		),
		'LifeStage LvKO p-val'                       =	modelSummaryPvalueExtract(
			x = object$output$SplitModels$Genotype_LifeStage,
			variable = CombineLevels(
				paste0(Labels$LifeStage$LifeStage, Labels$LifeStage$Late),
				Labels$Genotype$Levels,
				debug = debug
			),
			anova = FALSE,
			debug = debug
		),
		'LifeStage LvKO effect size'                 = object$output$EffectSizes$Genotype_Late ,
		#####################################################################
		################ Sex LifeStage Genotype interactions
		# 1.
		'LifeStageSexGenotype FvEvKO estimate'           =	modelSummaryPvalueExtract(
			x = object$output$SplitModels$Genotype_Sex.LifeStage,
			variable = 		CombineLevels(
				paste0(Labels$Sex$Sex, Labels$Sex$Female),
				paste0(Labels$LifeStage$LifeStage, Labels$LifeStage$Early),
				Labels$Genotype$Levels,
				debug = debug
			),
			anova = FALSE,
			what = 'Value',
			debug = debug
		),
		'LifeStageSexGenotype FvEvKO standard error'           =	modelSummaryPvalueExtract(
			x = object$output$SplitModels$Genotype_Sex.LifeStage,
			variable = 		CombineLevels(
				paste0(Labels$Sex$Sex, Labels$Sex$Female),
				paste0(Labels$LifeStage$LifeStage, Labels$LifeStage$Early),
				Labels$Genotype$Levels,
				debug = debug
			),
			anova = FALSE,
			what = 'Std.Error',
			debug = debug
		),
		'LifeStageSexGenotype FvEvKO p-val'           =	modelSummaryPvalueExtract(
			x = object$output$SplitModels$Genotype_Sex.LifeStage,
			variable = 		CombineLevels(
				paste0(Labels$Sex$Sex, Labels$Sex$Female),
				paste0(Labels$LifeStage$LifeStage, Labels$LifeStage$Early),
				Labels$Genotype$Levels,
				debug = debug
			),
			anova = FALSE,
			debug = debug
		),
		'LifeStageSexGenotype FvEvKO effect size'        = object$output$EffectSizes$Genotype_Female.Early ,
		# 2.
		'LifeStageSexGenotype MvEvKO estimate'           =	modelSummaryPvalueExtract(
			x = object$output$SplitModels$Genotype_Sex.LifeStage,
			variable = 		CombineLevels(
				paste0(Labels$Sex$Sex, Labels$Sex$Male),
				paste0(Labels$LifeStage$LifeStage, Labels$LifeStage$Early),
				Labels$Genotype$Levels,
				debug = debug
			),
			anova = FALSE,
			what = 'Value',
			debug = debug
		),
		'LifeStageSexGenotype MvEvKO standard error'           =	modelSummaryPvalueExtract(
			x = object$output$SplitModels$Genotype_Sex.LifeStage,
			variable = CombineLevels(
				paste0(Labels$Sex$Sex, Labels$Sex$Male),
				paste0(Labels$LifeStage$LifeStage, Labels$LifeStage$Early),
				Labels$Genotype$Levels,
				debug = debug
			),
			anova = FALSE,
			what = 'Std.Error',
			debug = debug
		),
		'LifeStageSexGenotype MvEvKO p-val'           =	modelSummaryPvalueExtract(
			x = object$output$SplitModels$Genotype_Sex.LifeStage,
			variable = CombineLevels(
				paste0(Labels$Sex$Sex, Labels$Sex$Male),
				paste0(Labels$LifeStage$LifeStage, Labels$LifeStage$Early),
				Labels$Genotype$Levels,
				debug = debug
			),
			anova = FALSE,
			debug = debug
		),
		'LifeStageSexGenotype MvEvKO effect size'        = object$output$EffectSizes$Genotype_Male.Early ,
		# 3.
		'LifeStageSexGenotype FvLvKO estimate'           =	modelSummaryPvalueExtract(
			x = object$output$SplitModels$Genotype_Sex.LifeStage,
			variable = CombineLevels(
				paste0(Labels$Sex$Sex, Labels$Sex$Female),
				paste0(Labels$LifeStage$LifeStage, Labels$LifeStage$Late),
				Labels$Genotype$Levels,
				debug = debug
			),
			anova = FALSE,
			what = 'Value',
			debug = debug
		),
		'LifeStageSexGenotype FvLvKO standard error'           =	modelSummaryPvalueExtract(
			x = object$output$SplitModels$Genotype_Sex.LifeStage,
			variable =  CombineLevels(
				paste0(Labels$Sex$Sex, Labels$Sex$Female),
				paste0(Labels$LifeStage$LifeStage, Labels$LifeStage$Late),
				Labels$Genotype$Levels,
				debug = debug
			),
			anova = FALSE,
			what = 'Std.Error',
			debug = debug
		),
		'LifeStageSexGenotype FvLvKO p-val'           =	modelSummaryPvalueExtract(
			x = object$output$SplitModels$Genotype_Sex.LifeStage,
			variable =  CombineLevels(
				paste0(Labels$Sex$Sex, Labels$Sex$Female),
				paste0(Labels$LifeStage$LifeStage, Labels$LifeStage$Late),
				Labels$Genotype$Levels,
				debug = debug
			),
			anova = FALSE,
			debug = debug
		),
		'LifeStageSexGenotype FvLvKO effect size'        = object$output$EffectSizes$Genotype_Female.Late ,
		#4.
		'LifeStageSexGenotype MvLvKO estimate'           =	modelSummaryPvalueExtract(
			x = object$output$SplitModels$Genotype_Sex.LifeStage,
			variable =  CombineLevels(
				paste0(Labels$Sex$Sex, Labels$Sex$Male),
				paste0(Labels$LifeStage$LifeStage, Labels$LifeStage$Late),
				Labels$Genotype$Levels,
				debug = debug
			),
			anova = FALSE,
			what = 'Value',
			debug = debug
		),
		'LifeStageSexGenotype MvLvKO standard error'           =	modelSummaryPvalueExtract(
			x = object$output$SplitModels$Genotype_Sex.LifeStage,
			variable = CombineLevels(
				paste0(Labels$Sex$Sex, Labels$Sex$Male),
				paste0(Labels$LifeStage$LifeStage, Labels$LifeStage$Late),
				Labels$Genotype$Levels,
				debug = debug
			),
			anova = FALSE,
			what = 'Std.Error',
			debug = debug
		),
		'LifeStageSexGenotype MvLvKO p-val'           =	modelSummaryPvalueExtract(
			x = object$output$SplitModels$Genotype_Sex.LifeStage,
			variable = CombineLevels(
				paste0(Labels$Sex$Sex, Labels$Sex$Male),
				paste0(Labels$LifeStage$LifeStage, Labels$LifeStage$Late),
				Labels$Genotype$Levels,
				debug = debug
			),
			anova = FALSE,
			debug = debug
		),
		'LifeStageSexGenotype MvLvKO effect size'        = object$output$EffectSizes$Genotype_Male.Late ,
		################
		'Classification tag'                   =	NULL,
		'Transformation'                       =	'Legacy',
		'Additional information'               =	addInfo
	)
	
	return(vectorOutput)
}
