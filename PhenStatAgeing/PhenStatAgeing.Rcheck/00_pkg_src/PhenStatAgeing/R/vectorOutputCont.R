vectorOutputCont =	function(object, json = FALSE)
{
	if (is.null(object))
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
	framework      = switch(
		object$output$Final.Model.Tag                                  ,
		lme  = "Linear Mixed Model framework"                          ,
		gls  = "Linear Model Using Generalized Least Squares framework",
		glm  = "Generalized Linear Model framework"
	)
	fittingMethod    = toupper(object$output$Final.Model.Tag)
	#####################################################################
	x                = object$input$data
	columnOfInterest = x[, c(depVariable)]
	#####################################################################
	variability      = length(unique(columnOfInterest)) / length(columnOfInterest)
	#####################################################################
	DSsize            = SummaryStats(
		x = x,
		formula = object$input$Fullfixed$initial,
		label = 'Summary statistics',
		lower = TRUE,
		drop = TRUE,
		sep = ' '
	)
	MultiBatch = ifelse(multiBatch(x),
											'Dataset contains multi batches',
											'Dataset contains single batch')
	
	OtherModels        = list(
		'Further model(s)' = lapply(object$output$SplitModels, function(v) {
			if (class(v) %in% c('lme', 'gls', 'glm')) {
				as.list(unmatrix0(summary(v)$tTable))
			} else{
				v
			}
		}),
		'Effect size(s)'                   = object$output$EffectSizes,
		'Other residual normality test(s)' = object$output$ResidualNormalityTests
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
				term = Labels$Sex$Sex,
				return = NA,
				not = modelSummaryPvalueExtract(
					x = Fmodel,
					variable = Labels$Genotype$Genotype,
					anova = TRUE
				)
			),
			'Sex FvKO p-val'   =	TermInFormulaReturn(
				active = TRUE,
				formula = frm,
				term = Labels$Sex$Sex,
				not = NA,
				return = modelSummaryPvalueExtract(
					x = object$output$SplitModels$Genotype_Sex,
					# SexFemale:Genotypeexperimental
					variable = CombineLevels(paste0('Sex', Labels$Sex$Female), Labels$Genotype$Levels),
					anova = FALSE
				)
			),
			'Sex MvKO p-val'  =	TermInFormulaReturn(
				active = TRUE,
				formula = frm,
				term = Labels$Sex$Sex,
				not = NA,
				return = modelSummaryPvalueExtract(
					x = object$output$SplitModels$Genotype_Sex,
					variable = CombineLevels(paste0('Sex', Labels$Sex$Male), Labels$Genotype$Levels),
					anova = FALSE
				)
			),
			'Sexual dimorphism detected' = TermInFormulaReturn(
				active = TRUE,
				formula = frm,
				term = Labels$Sex$Sex,
				return = TRUE,
				not = FALSE
			)
		),
		'Genotype estimate' =
			modelSummaryPvalueExtract(
				x = Fmodel,
				variable = unlist(Labels$Genotype$Levels),
				anova = FALSE,
				what = 'Value'
			),
		'Genotype standard error'  =
			modelSummaryPvalueExtract(
				x = Fmodel,
				variable = unlist(Labels$Genotype$Levels),
				anova = FALSE,
				what = 'Std.Error'
			),
		'Genotype p-val'           =
			modelSummaryPvalueExtract(
				x = Fmodel,
				variable = Labels$Genotype$Genotype,
				anova = TRUE
			),
		'Genotype percentage change'           =	percentageChanges,
		'Genotype effect size'                 = object$output$EffectSizes[[Labels$Genotype$Genotype]],
		#####################################################################
		'Sex estimate'                         =	modelSummaryPvalueExtract(
			x = Fmodel,
			variable = unlist(Labels$Sex$Levels),
			anova = FALSE,
			what = 'Value'
		),
		'Sex standard error'                   = modelSummaryPvalueExtract(
			x = Fmodel,
			variable = unlist(Labels$Sex$Levels),
			anova = FALSE,
			what = 'Std.Error'
		),
		'Sex p-val'                            =	modelSummaryPvalueExtract(
			x = Fmodel,
			variable = Labels$Sex$Sex,
			anova = TRUE
		),
		'Sex effect size'                         =	object$output$EffectSizes[[Labels$Sex$Sex]],
		#####################################################################
		'LifeStage estimate'                      =	modelSummaryPvalueExtract(
			x = Fmodel,
			variable = unlist(Labels$LifeStage$Levels),
			anova = FALSE,
			what = 'Value'
		),
		'LifeStage standard error'                =	modelSummaryPvalueExtract(
			x = Fmodel,
			variable = unlist(Labels$LifeStage$Levels),
			anova = FALSE,
			what = 'Std.Error'
		),
		'LifeStage p-val'                         =	modelSummaryPvalueExtract(
			x = Fmodel,
			variable = Labels$LifeStage$LifeStage,
			anova = TRUE
		),
		'LifeStage effect size'                = object$output$EffectSizes[[Labels$LifeStage$LifeStage]] ,
		#####################################################################
		'Weight estimate'                      =	modelSummaryPvalueExtract(
			x = Fmodel,
			variable = Labels$Weight,
			anova = FALSE,
			what = 'Value'
		),
		'Weight standard error'                =	modelSummaryPvalueExtract(
			x = Fmodel,
			variable = Labels$Weight,
			anova = FALSE,
			what = 'Std.Error'
		),
		'Weight p-val'                         =	modelSummaryPvalueExtract(
			x = Fmodel,
			variable = Labels$Weight,
			anova = TRUE
		),
		'Weight effect size'                   =  object$output$EffectSizes[[Labels$Weight]],
		#####################################################################
		'Gp1 genotype'                         =	Labels$Genotype$Control		,
		'Gp1 Residuals normality test'         =	as.numeric01(object$output$ResidualNormalityTests$Genotype[Labels$Genotype$Control]),
		'Gp2 genotype'                         =	Labels$Genotype$Mutant				,
		'Gp2 Residuals normality test'         =	as.numeric01(
			object$output$ResidualNormalityTests$Genotype[Labels$Genotype$Mutant]
		),
		#####################################################################
		'Blups test'                           =  'Legacy',
		'Rotated residuals normality test'     =  'Legacy',
		#####################################################################
		'Intercept estimate'                   =	modelSummaryPvalueExtract(
			x = Fmodel,
			variable = '(Intercept)',
			anova = FALSE,
			what = 'Value'
		),
		'Intercept standard error'             =	modelSummaryPvalueExtract(
			x = Fmodel,
			variable = '(Intercept)',
			anova = FALSE,
			what = 'Std.Error'
		),
		'Intercept p-val'                      =	modelSummaryPvalueExtract(
			x = Fmodel,
			variable = '(Intercept)',
			anova = TRUE
		),
		#####################################################################
		'Interaction(s) included'              =	list(
			'Genotype Sex'  =  !is.null(
				modelSummaryPvalueExtract(
					x = Fmodel,
					variable = paste(Labels$Genotype$Genotype, Labels$Sex$Sex, sep = ':'),
					anova = TRUE
				)
			),
			'Genotype LifeStage'  =  !is.null(
				modelSummaryPvalueExtract(
					x = Fmodel,
					variable = paste(Labels$Genotype$Genotype, Labels$LifeStage$LifeStage, sep = ':'),
					anova = TRUE
				)
			),
			'Sex LifeStage'  =  !is.null(
				modelSummaryPvalueExtract(
					x = Fmodel,
					variable = paste(Labels$Sex$Sex, Labels$LifeStage$LifeStage, sep = ':'),
					anova = TRUE
				)
			),
			'Genotype Sex LifeStage'  =  !is.null(
				modelSummaryPvalueExtract(
					x = Fmodel,
					variable = paste(
						Labels$Genotype$Genotype,
						Labels$Sex$Sex,
						Labels$LifeStage$LifeStage,
						sep = ':'
					),
					anova = TRUE
				)
			)
		),
		#####################################################################
		################ interaction
		'Interaction(s) p-val'      =	list(
			'Genotype Sex'  =
				modelSummaryPvalueExtract(
					x = Fmodel,
					variable =  paste(Labels$Genotype$Genotype, Labels$Sex$Sex, sep = ':'),
					anova = TRUE
				),
			'Genotype LifeStage'  =
				modelSummaryPvalueExtract(
					x = Fmodel,
					variable =  paste(Labels$Genotype$Genotype, Labels$LifeStage$LifeStage, sep = ':'),
					anova = TRUE
				),
			'Sex LifeStage'  =
				modelSummaryPvalueExtract(
					x = Fmodel,
					variable =  paste(Labels$Sex$Sex, Labels$LifeStage$LifeStage, sep = ':'),
					anova = TRUE
				),
			'Genotype Sex LifeStage'  =
				modelSummaryPvalueExtract(
					x = Fmodel,
					variable =  paste(
						Labels$Genotype$Genotype,
						Labels$Sex$Sex,
						Labels$LifeStage$LifeStage,
						sep = ':'
					),
					anova = TRUE
				)
		),
		#####################################################################
		################ Sex interactions
		'Sex FvKO estimate'                    =
			modelSummaryPvalueExtract(
				x = object$output$SplitModels$Genotype_Sex,
				variable = CombineLevels(paste0('Sex', Labels$Sex$Female), Labels$Genotype$Levels),
				anova = FALSE,
				what = 'Value'
			),
		'Sex FvKO standard error'              =	modelSummaryPvalueExtract(
			x = object$output$SplitModels$Genotype_Sex,
			variable = CombineLevels(paste0('Sex', Labels$Sex$Female), Labels$Genotype$Levels),
			anova = FALSE,
			what = 'Std.Error'
		),
		'Sex FvKO p-val'                       = modelSummaryPvalueExtract(
			x = object$output$SplitModels$Genotype_Sex,
			variable = CombineLevels(paste0('Sex', Labels$Sex$Female), Labels$Genotype$Levels),
			anova = FALSE
		),
		'Sex FvKO effect size'                 = object$output$EffectSizes[[paste(Labels$Genotype$Genotype, Labels$Sex$Female, sep = '_')]],
		#####################################################################
		'Sex MvKO estimate'                    =	modelSummaryPvalueExtract(
			x = object$output$SplitModels$Genotype_Sex,
			variable = CombineLevels(paste0('Sex', Labels$Sex$Male), Labels$Genotype$Levels),
			anova = FALSE,
			what = 'Value'
		),
		'Sex MvKO standard error'              =	 modelSummaryPvalueExtract(
			x = object$output$SplitModels$Genotype_Sex,
			variable = CombineLevels(paste0('Sex', Labels$Sex$Male), Labels$Genotype$Levels),
			anova = FALSE,
			what = 'Std.Error'
		),
		'Sex MvKO p-val'                       =	modelSummaryPvalueExtract(
			x = object$output$SplitModels$Genotype_Sex,
			variable = CombineLevels(paste0('Sex', Labels$Sex$Male), Labels$Genotype$Levels),
			anova = FALSE
		),
		'Sex MvKO effect size'                 = object$output$EffectSizes[[paste(Labels$Genotype$Genotype, Labels$Sex$Male, sep = '_')]],
		#####################################################################
		################ LifeStage interaction
		'LifeStage EvKO estimate'                    =	modelSummaryPvalueExtract(
			x = object$output$SplitModels$Genotype_LifeStage,
			variable = CombineLevels(
				paste0(Labels$LifeStage$LifeStage, Labels$LifeStage$Early),
				Labels$Genotype$Levels
			),
			anova = FALSE,
			what = 'Value'
		),
		'LifeStage EvKO standard error'              =	modelSummaryPvalueExtract(
			x = object$output$SplitModels$Genotype_LifeStage,
			variable =  CombineLevels(
				paste0(Labels$LifeStage$LifeStage, Labels$LifeStage$Early),
				Labels$Genotype$Levels
			),
			anova = FALSE,
			what = 'Std.Error'
		),
		'LifeStage EvKO p-val'                       =	modelSummaryPvalueExtract(
			x = object$output$SplitModels$Genotype_LifeStage,
			variable =  CombineLevels(
				paste0(Labels$LifeStage$LifeStage, Labels$LifeStage$Early),
				Labels$Genotype$Levels
			),
			anova = FALSE
		),
		'LifeStage EvKO effect size'                 = object$output$EffectSizes$Genotype_Early ,
		#####################################################################
		'LifeStage LvKO estimate'                    =	modelSummaryPvalueExtract(
			x = object$output$SplitModels$Genotype_LifeStage,
			variable =  CombineLevels(
				paste0(Labels$LifeStage$LifeStage, Labels$LifeStage$Late),
				Labels$Genotype$Levels
			),
			anova = FALSE,
			what = 'Value'
		),
		'LifeStage LvKO standard error'              =	modelSummaryPvalueExtract(
			x = object$output$SplitModels$Genotype_LifeStage,
			variable = CombineLevels(
				paste0(Labels$LifeStage$LifeStage, Labels$LifeStage$Late),
				Labels$Genotype$Levels
			),
			anova = FALSE,
			what = 'Std.Error'
		),
		'LifeStage LvKO p-val'                       =	modelSummaryPvalueExtract(
			x = object$output$SplitModels$Genotype_LifeStage,
			variable = CombineLevels(
				paste0(Labels$LifeStage$LifeStage, Labels$LifeStage$Late),
				Labels$Genotype$Levels
			),
			anova = FALSE
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
				Labels$Genotype$Levels
			),
			anova = FALSE,
			what = 'Value'
		),
		'LifeStageSexGenotype FvEvKO standard error'           =	modelSummaryPvalueExtract(
			x = object$output$SplitModels$Genotype_Sex.LifeStage,
			variable = 		CombineLevels(
				paste0(Labels$Sex$Sex, Labels$Sex$Female),
				paste0(Labels$LifeStage$LifeStage, Labels$LifeStage$Early),
				Labels$Genotype$Levels
			),
			anova = FALSE,
			what = 'Std.Error'
		),
		'LifeStageSexGenotype FvEvKO p-val'           =	modelSummaryPvalueExtract(
			x = object$output$SplitModels$Genotype_Sex.LifeStage,
			variable = 		CombineLevels(
				paste0(Labels$Sex$Sex, Labels$Sex$Female),
				paste0(Labels$LifeStage$LifeStage, Labels$LifeStage$Early),
				Labels$Genotype$Levels
			),
			anova = FALSE
		),
		'LifeStageSexGenotype FvEvKO effect size'        = object$output$EffectSizes$Genotype_Female.Early ,
		# 2.
		'LifeStageSexGenotype MvEvKO estimate'           =	modelSummaryPvalueExtract(
			x = object$output$SplitModels$Genotype_Sex.LifeStage,
			variable = 		CombineLevels(
				paste0(Labels$Sex$Sex, Labels$Sex$Male),
				paste0(Labels$LifeStage$LifeStage, Labels$LifeStage$Early),
				Labels$Genotype$Levels
			),
			anova = FALSE,
			what = 'Value'
		),
		'LifeStageSexGenotype MvEvKO standard error'           =	modelSummaryPvalueExtract(
			x = object$output$SplitModels$Genotype_Sex.LifeStage,
			variable = CombineLevels(
				paste0(Labels$Sex$Sex, Labels$Sex$Male),
				paste0(Labels$LifeStage$LifeStage, Labels$LifeStage$Early),
				Labels$Genotype$Levels
			),
			anova = FALSE,
			what = 'Std.Error'
		),
		'LifeStageSexGenotype MvEvKO p-val'           =	modelSummaryPvalueExtract(
			x = object$output$SplitModels$Genotype_Sex.LifeStage,
			variable = CombineLevels(
				paste0(Labels$Sex$Sex, Labels$Sex$Male),
				paste0(Labels$LifeStage$LifeStage, Labels$LifeStage$Early),
				Labels$Genotype$Levels
			),
			anova = FALSE
		),
		'LifeStageSexGenotype MvEvKO effect size'        = object$output$EffectSizes$Genotype_Male.Early ,
		# 3.
		'LifeStageSexGenotype FvLvKO estimate'           =	modelSummaryPvalueExtract(
			x = object$output$SplitModels$Genotype_Sex.LifeStage,
			variable = CombineLevels(
				paste0(Labels$Sex$Sex, Labels$Sex$Female),
				paste0(Labels$LifeStage$LifeStage, Labels$LifeStage$Late),
				Labels$Genotype$Levels
			),
			anova = FALSE,
			what = 'Value'
		),
		'LifeStageSexGenotype FvLvKO standard error'           =	modelSummaryPvalueExtract(
			x = object$output$SplitModels$Genotype_Sex.LifeStage,
			variable =  CombineLevels(
				paste0(Labels$Sex$Sex, Labels$Sex$Female),
				paste0(Labels$LifeStage$LifeStage, Labels$LifeStage$Late),
				Labels$Genotype$Levels
			),
			anova = FALSE,
			what = 'Std.Error'
		),
		'LifeStageSexGenotype FvLvKO p-val'           =	modelSummaryPvalueExtract(
			x = object$output$SplitModels$Genotype_Sex.LifeStage,
			variable =  CombineLevels(
				paste0(Labels$Sex$Sex, Labels$Sex$Female),
				paste0(Labels$LifeStage$LifeStage, Labels$LifeStage$Late),
				Labels$Genotype$Levels
			),
			anova = FALSE
		),
		'LifeStageSexGenotype FvLvKO effect size'        = object$output$EffectSizes$Genotype_Female.Late ,
		#4.
		'LifeStageSexGenotype MvLvKO estimate'           =	modelSummaryPvalueExtract(
			x = object$output$SplitModels$Genotype_Sex.LifeStage,
			variable =  CombineLevels(
				paste0(Labels$Sex$Sex, Labels$Sex$Male),
				paste0(Labels$LifeStage$LifeStage, Labels$LifeStage$Late),
				Labels$Genotype$Levels
			),
			anova = FALSE,
			what = 'Value'
		),
		'LifeStageSexGenotype MvLvKO standard error'           =	modelSummaryPvalueExtract(
			x = object$output$SplitModels$Genotype_Sex.LifeStage,
			variable = CombineLevels(
				paste0(Labels$Sex$Sex, Labels$Sex$Male),
				paste0(Labels$LifeStage$LifeStage, Labels$LifeStage$Late),
				Labels$Genotype$Levels
			),
			anova = FALSE,
			what = 'Std.Error'
		),
		'LifeStageSexGenotype MvLvKO p-val'           =	modelSummaryPvalueExtract(
			x = object$output$SplitModels$Genotype_Sex.LifeStage,
			variable = CombineLevels(
				paste0(Labels$Sex$Sex, Labels$Sex$Male),
				paste0(Labels$LifeStage$LifeStage, Labels$LifeStage$Late),
				Labels$Genotype$Levels
			),
			anova = FALSE
		),
		'LifeStageSexGenotype MvLvKO effect size'        = object$output$EffectSizes$Genotype_Male.Late ,
		################
		'Classification tag'                   =	NA,
		'Transformation'                       =	'Legacy',
		'Additional information'               =	addInfo
	)
	
	if (json) {
		for (i in 1:10) {
			vectorOutput = toJSON(vectorOutput, auto_unbox = TRUE, null = 'null',digits = 500)
			if (i != 10)
				vectorOutput = fromJSON(txt = vectorOutput)
		}
		write(vectorOutput, file = 'd:/jsontestCon.txt', ncolumns = 10 ^ 3,append = FALSE)
	}
	return(vectorOutput)
}
