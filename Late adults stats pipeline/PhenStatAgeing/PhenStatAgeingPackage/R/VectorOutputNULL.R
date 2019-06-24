vectorOutputNULL = function(object) {
	VO      = list(
		'Method'                               = 	NULL,
		'Dependent variable'                   =	NULL,
		'Batch included'                       =	NULL,
		'Batch p-val'                          =  NULL,
		'Residual variances homogeneity'       =  NULL,
		'Residual variances homogeneity p-val' =  NULL,
		#####################################################################
		'Genotype contribution' =	list(
			Overal             = NULL,
			'Sex FvKO p-val'   = NULL,
			'Sex MvKO p-val'   = NULL	,
			'Sexual dimorphism detected' = NULL
		),
		'Genotype estimate'           = NULL ,
		'Genotype standard error'     = NULL ,
		'Genotype p-val'              = NULL ,
		'Genotype percentage change'  =	NULL ,
		'Genotype effect size'        = NULL ,
		#####################################################################
		'Sex estimate'                =	NULL ,
		'Sex standard error'          = NULL ,
		'Sex p-val'                   =	NULL ,
		'Sex effect size'             =	NULL ,
		#####################################################################
		'LifeStage estimate'          =	NULL ,
		'LifeStage standard error'    =	NULL ,
		'LifeStage p-val'             =	NULL ,
		'LifeStage effect size'       = NULL ,
		#####################################################################
		'Weight estimate'             =	NULL ,
		'Weight standard error'       =	NULL ,
		'Weight p-val'                =	NULL ,
		'Weight effect size'          = NULL ,
		#####################################################################
		'Gp1 genotype'                     =	NULL ,
		'Gp1 Residuals normality test'     =	NULL ,
		'Gp2 genotype'                     =	NULL ,
		'Gp2 Residuals normality test'     =	NULL ,
		#####################################################################
		'Blups test'                       =  NULL,
		'Rotated residuals normality test' =  NULL,
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
		'Sex FvKO estimate'         = NULL ,
		'Sex FvKO standard error'   = NULL ,
		'Sex FvKO p-val'            = NULL ,
		'Sex FvKO effect size'      = NULL ,
		#####################################################################
		'Sex MvKO estimate'         = NULL ,
		'Sex MvKO standard error'   = NULL ,
		'Sex MvKO p-val'            = NULL ,
		'Sex MvKO effect size'      = NULL ,
		#####################################################################
		################ LifeStage interaction
		'LifeStage EvKO estimate'         =	NULL ,
		'LifeStage EvKO standard error'   =	NULL ,
		'LifeStage EvKO p-val'            =	NULL ,
		'LifeStage EvKO effect size'      = NULL ,
		#####################################################################
		'LifeStage LvKO estimate'         =	NULL ,
		'LifeStage LvKO standard error'   =	NULL ,
		'LifeStage LvKO p-val'            =	NULL ,
		'LifeStage LvKO effect size'      = NULL ,
		#####################################################################
		################ Sex LifeStage Genotype interactions
		# 1.
		'LifeStageSexGenotype FvEvKO estimate'        =	NULL ,
		'LifeStageSexGenotype FvEvKO standard error'  =	NULL ,
		'LifeStageSexGenotype FvEvKO p-val'           =	NULL ,
		'LifeStageSexGenotype FvEvKO effect size'     = NULL ,
		# 2.
		'LifeStageSexGenotype MvEvKO estimate'        =	NULL ,
		'LifeStageSexGenotype MvEvKO standard error'  =	NULL ,
		'LifeStageSexGenotype MvEvKO p-val'           =	NULL ,
		'LifeStageSexGenotype MvEvKO effect size'     = NULL ,
		# 3.
		'LifeStageSexGenotype FvLvKO estimate'        = NULL ,
		'LifeStageSexGenotype FvLvKO standard error'  = NULL ,
		'LifeStageSexGenotype FvLvKO p-val'           = NULL ,
		'LifeStageSexGenotype FvLvKO effect size'     = NULL ,
		
		'LifeStageSexGenotype MvLvKO estimate'        = NULL ,
		'LifeStageSexGenotype MvLvKO standard error'  = NULL ,
		'LifeStageSexGenotype MvLvKO p-val'           =	NULL ,
		'LifeStageSexGenotype MvLvKO effect size'     = NULL ,
		################
		'Classification tag'                          =	NULL ,
		'Transformation'                              =	NULL ,
		'Additional information'                      =	object$messages
	)
	return(VO)
}
